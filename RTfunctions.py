import os, sys, numpy as np, pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.vectors import FloatVector
from pyteomics import mzxml


class progressBar:
    def __init__(self, total):
        self.total = total
        self.barLength = 20
        self.count = 0
        self.progress = 0
        self.block = 0
        self.status = ""

    def increment(self, nIncrement=None):
        if nIncrement == None:
            self.count += 1
        else:
            self.count = nIncrement
        self.progress = self.count / self.total
        self.block = int(round(self.barLength * self.progress))
        if self.progress == 1:
            self.status = "Done...\r\n"
        else:
            self.status = ""
        #         self.status = str(self.count) + "/" + str(self.total)
        text = "\r  Progress: [{0}] {1}% {2}".format("#" * self.block + "-" * (self.barLength - self.block),
                                                     int(self.progress * 100), self.status)
        sys.stdout.write(text)
        sys.stdout.flush()


def loess():
    rstring = """
    loess.as = function(x, y, degree = 1, criterion="aicc", family="gaussian", user.span=NULL, plot=FALSE, ...) {

        criterion <- match.arg(criterion)
        family <- match.arg(family)
        x <- as.matrix(x)

        if ((ncol(x) != 1) & (ncol(x) != 2)) stop("The predictor 'x' should be one or two dimensional!!")
        if (!is.numeric(x)) stop("argument 'x' must be numeric!")
        if (!is.numeric(y)) stop("argument 'y' must be numeric!")
        if (any(is.na(x))) stop("'x' contains missing values!")
        if (any(is.na(y))) stop("'y' contains missing values!")
        if (!is.null(user.span) && (length(user.span) != 1 || !is.numeric(user.span))) 
            stop("argument 'user.span' must be a numerical number!")
        if(nrow(x) != length(y)) stop("'x' and 'y' have different lengths!")
        if(length(y) < 3) stop("not enough observations!")

        data.bind <- data.frame(x=x, y=y)
        if (ncol(x) == 1) {
            names(data.bind) <- c("x", "y")
        } else { names(data.bind) <- c("x1", "x2", "y") }

        opt.span <- function(model, criterion=c("aicc", "gcv"), span.range=c(.05, .95)){	
            as.crit <- function (x) {
                span <- x$pars$span
                traceL <- x$trace.hat
                sigma2 <- sum(x$residuals^2 ) / (x$n-1)
                aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (x$n-traceL-2)
                gcv <- x$n*sigma2 / (x$n-traceL)^2
                result <- list(span=span, aicc=aicc, gcv=gcv)
                return(result)
            }
            criterion <- match.arg(criterion)
            fn <- function(span) {
                mod <- update(model, span=span)
                as.crit(mod)[[criterion]]
            }
            result <- optimize(fn, span.range)
            return(list(span=result$minimum, criterion=result$objective))
        }

        control = loess.control(surface = "direct")
        if (ncol(x)==1) {
            if (is.null(user.span)) {
                fit0 <- loess(y ~ x, degree=degree, family=family, data=data.bind, control=control, ...)
                span1 <- opt.span(fit0, criterion=criterion)$span
            } else {
                span1 <- user.span
            }		
            fit <- loess(y ~ x, degree=degree, span=span1, family=family, data=data.bind, control=control, ...)
        } else {
            if (is.null(user.span)) {
                fit0 <- loess(y ~ x1 + x2, degree=degree,family=family, data.bind, control=control, ...)
                span1 <- opt.span(fit0, criterion=criterion)$span
            } else {
                span1 <- user.span
            }		
            fit <- loess(y ~ x1 + x2, degree=degree, span=span1, family=family, data=data.bind, control=control...)
        }
        return(fit)
    }
    """
    return ro.r(rstring)


def getMs2ToSurvey(reader):
    # Read mzXML file and make a relationship between every MS2 and its survey scan
    res = {}
    for spec in reader:
        if spec["msLevel"] == 1:
            survey = int(spec["num"])
        elif spec["msLevel"] == 2:
            res[int(spec["num"])] = survey
    return res


def getPrecursorPeak(reader, psmScanNumber, surveyScanNumber, isoWindow):
    # Find the precursor peak of the PSM (the strongest peak within the isolation window)
    nominalPrecMz = float(reader[str(psmScanNumber)]["filterLine"].split("@")[0].split(" ")[-1])
    spec = reader[str(surveyScanNumber)]
    mzArray = spec["m/z array"]
    intArray = spec["intensity array"]
    ind = (mzArray >= nominalPrecMz - 0.5 * isoWindow) & (mzArray <= nominalPrecMz + 0.5 * isoWindow)
    if sum(ind > 0):
        subMzArray = mzArray[ind]
        subIntArray = intArray[ind]
        ind2 = np.argmax(subIntArray)
        precMz = subMzArray[ind2]
        precIntensity = subIntArray[ind2]
    else:
        precMz = -1
        precIntensity = -1

    precRt = spec["retentionTime"]   # Unit of minute
    return precMz, precIntensity, precRt


def formatRtTable(df):
    # Input dataframe is from "getRT" function with columns ["key", "run", "RT"]
    # Output is the dataframe with columns ["key", <run#1>, <run#2>, ... <run#n>]
    # Each column, <run#i>, contains RT of the key in the run#i
    # Re-organization of the data frame
    keys = list(df["key"].drop_duplicates())
    runs = list(df["run"].drop_duplicates())
    res = pd.DataFrame({"key": keys})
    for run in runs:
        res[run] = np.nan
        res[run + "_nPSMs"] = np.nan

    for i in range(res.shape[0]):
        key = res.loc[i]["key"]
        idx = df["key"] == key
        runs = df[idx]["run"].values
        rts = df[idx]["RT"].values
        npsms = df[idx]["nPSMs"].values
        for j in range(len(runs)):
            res.loc[i, runs[j]] = rts[j]
            res.loc[i, runs[j] + "_nPSMs"] = npsms[j]

    return res


# # Suresh's revision for speeding up
# def formatRtTable(df):
#     df_nPSMs = df.set_index(['key', 'run']).nPSMs.unstack().reset_index()
#     df_RT = df.set_index(['key', 'run']).RT.unstack().reset_index()
#     df_nPSMs2 = df_nPSMs.set_index("key")
#     df_RT2 = df_RT.set_index("key")
#     col_keys_nPSM ={}
#     for val in df_nPSMs2.columns:
#         new_val = val + "_nPSMs"
#         col_keys_nPSM[val] = new_val
#
#     df_nPSMs3 = df_nPSMs2.rename(columns=col_keys_nPSM)
#     res = pd.concat([df_RT2, df_nPSMs3], axis=1)
#     res.reset_index(inplace=True)
#     return res



def inferRT(idTxt, runs, isolationWindow):
    # Input
    # 1. mzXML files
    # 2. ID.txt file containing all identified PSMs
    print("  Extraction and assignment of RTs to the identified PSMs")
    print("  =======================================================")

    # Read ID.txt files to extract PSM information
    print("  Read ID.txt file: to extract PSM information")
    psms = pd.read_csv(idTxt, skiprows=1, sep=";")  # Note that ID.txt file is delimited by semicolon
    psms = psms[["Peptide", "Outfile", "XCorr"]].drop_duplicates()
    psms["charge"] = [outfile.split("/")[-1].split(".")[-2] for outfile in psms["Outfile"]]
    psms["key"] = psms["Peptide"] + "_" + psms["charge"]
    print("  Done ...")

    # RT extraction/assignment for each mzXML file
    res = []
    for run in runs:
        runName = os.path.basename(run).split(".")[0]
        print("  RT of every identified peptide in {} is being inferred and assigned".format(runName))

        # Read a mzXML file and extract PSMs corresponding to the mzXML file
        reader = mzxml.MzXML(run)
        ms2ToSurvey = getMs2ToSurvey(reader)
        subPsms = psms[psms["Outfile"].str.contains(runName + "/")]

        # Unique key is peptide-charge pair
        keys = subPsms["key"]
        keys = list(set(keys))
        progress = progressBar(len(keys))
        for key in keys:
            progress.increment()
            rtArray = np.array([])
            intArray = np.array([])
            for _, psm in subPsms[subPsms["key"] == key].iterrows():
                [_, psmScanNum, _, _, _] = os.path.basename(psm["Outfile"]).split(".")
                psmScanNum = int(psmScanNum)
                surveyScanNum = ms2ToSurvey[psmScanNum]
                _, precIntensity, precRt = getPrecursorPeak(reader, int(psmScanNum), surveyScanNum, isolationWindow)
                rtArray = np.append(rtArray, precRt)
                intArray = np.append(intArray, precIntensity)
            rt = sum(rtArray * intArray) / sum(intArray)   # Unit of minute
            res.append([key, runName, rt, len(rtArray)])

    print()
    res = pd.DataFrame(res, columns=["key", "run", "RT", "nPSMs"])
    res = formatRtTable(res)
    return res


def alignRT(df, refRunIdx=0):
    # Input: df = a pandas dataframe containing keys (peptide-charge pairs) and RTs over the runs (i.e., mzXML files)
    #        refRunIdx = the index of the reference run (by default, 0, the first run)

    print("  Alignment and calibration of RTs over multiple runs")
    print("  ===================================================")
    # Initialization
    res = df.copy()
    rLoess = loess()
    rPredict = ro.r("predict")
    colN = [c for c in df.columns if c.endswith("nPSMs")]
    colRt = df.columns.drop(colN)
    colRt = colRt.drop("key")

    # RT-alignment using a specific run as a reference (sequential alignment)
    # This approach is composed of three steps
    # 1. Set a specific run as a reference
    # 2. Align and calibrate one of remaining run against the reference using a LOESS model
    # 3. Update the reference by merging the old reference and the aligned/calibrated run (using weighted average)

    # Set the initial reference
    refIdx = refRunIdx
    ref = df[colRt[refIdx]].copy()
    refN = df[colN[refIdx]].copy()
    print("  {} is set to the initial reference for RT-alignment".format(colRt[refIdx]))

    # Alignment and calibration of RT
    for i in range(len(colRt)):
        if i == refIdx:  # Skip for the reference run
            continue
        print("  The RTs in {} are being aligned and calibrated".format(colRt[i]))

        # Prepare the modeling: find the shared peptides between the reference and the current run
        idx = (~ref.isna()) & (~df[colRt[i]].isna())
        x = ref[idx]
        y = df[idx][colRt[i]]

        # Build a LOESS model and calibrate RTs
        mod = rLoess(FloatVector(x), FloatVector(y))  # LOESS model based on the shared peptides
        res[colRt[i]] = rPredict(mod, FloatVector(df[colRt[i]]))  # Calibration is applied to whole peptides of the current run

        # Update the reference by merging (taking the union of) the current reference and the calibrated run
        # There are three types of peptides (i.e., row indices)
        # 1. Peptides shared between the current reference and the calibrated run -> update the reference by taking the weighted average
        fullRt = pd.concat([ref[idx], res[colRt[i]][idx]], axis=1)
        fullN = pd.concat([refN[idx], res[colN[i]][idx]], axis=1)
        rt = (fullRt * fullN.values).sum(axis=1) / fullN.sum(axis=1)  # Reference as the weighted mean of RTs
        ref[idx] = rt
        refN[idx] = refN[idx] + res[colN[i]][idx]
        # 2. Peptides only found in the calibrated run -> replace the reference with the calibrated RTs
        idx = (ref.isna()) & (~df[colRt[i]].isna())
        ref[idx] = res[colRt[i]][idx]
        refN[idx] = res[colN[i]][idx]
        # 3. Peptides only found  in the current reference -> no action

    # Organization of the output dataframe
    # Calculation of the weighted standard deviation of RTs (https://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf)
    M = (~res[colN].isna()).sum(axis=1)
    den = ((M - 1) / M) * res[colN].sum(axis=1)
    num = ((res[colRt].sub(ref, axis=0) ** 2) * res[colN].values).sum(axis=1)
    sdRt = np.sqrt(num / den)
    sdRt[den == 0] = 0
    res["SdRT"] = sdRt
    # Calculation of the weighted average RTs
    res["AvgRT"] = ref  # In fact, the final reference is equivalent to the weighted average of the aligned/calibrated RTs

    print("  Done ...\n")
    return res


def getRT(idTxt, runs):
    # Input: idTxt = the path of "ID.txt" file containing the information of PSMs
    #        runs = an array containing the paths of mzXML files

    # Only parameter (users should specify) is the "isolation window"
    # Here, it is arbitrarily set to 1 -> +/- 0.5 m/z will be examined to search for the precursor peak
    isolationWindow = 1

    # Extraction of RTs for the precursor peaks of PSMs, and assignment/calculation of RTs to PSMs
    df = inferRT(idTxt, runs, isolationWindow)

    # Alignment and calibration of RTs over the "runs"
    res = alignRT(df)

    return res