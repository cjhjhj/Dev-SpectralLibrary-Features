import sys, os, re, numpy as np, pandas as pd
import matplotlib.pyplot as plt
import rpy2.robjects as ro
from rpy2.robjects.vectors import FloatVector
from RTassignment import getRT


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

def getAlignedRT(mode, refRunIdx=0):
    mzXML = [r"../Data/FTLD_Batch2_F50.mzXML", r"../Data/FTLD_Batch2_F51.mzXML", r"../Data/FTLD_Batch2_F52.mzXML"]
    idTxt = "../Data/ID.txt"
    # df = getRT(mzXML, idTxt)
    df = pd.read_csv("rtAlignment.txt", sep="\t")

    # Find the keys fully identified across the runs
        # To-do: almost fully identified keys may be used

    # Initialization
    res = df.copy()
    rLoess = loess()
    rPredict = ro.r("predict")
    colN = [c for c in df.columns if c.endswith("nPSMs")]
    colRt = df.columns.drop(colN)
    colRt = colRt.drop("key")

    if mode == 1:   # RT-alignment based on the fully-identified peptides/PSMs
        # 1. RT-alignment using the reference derived from weighted average RTs
        fullRt = df[df[colRt].isna().sum(axis=1) == 0][colRt]
        fullN = df[df[colRt].isna().sum(axis=1) == 0][colN]
        y = (fullRt * fullN.values).sum(axis=1) / fullN.sum(axis=1) # Reference as the weighted mean of RTs
        for col in colRt:
            x = fullRt[col]
            mod = rLoess(FloatVector(x), FloatVector(y))
            res[col] = rPredict(mod, FloatVector(df[col]))
    elif mode == 2:
        # 2. RT-alignment using a specific run as a reference (sequential alignment)
        # Junmin's preference

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
            if i == refIdx: # Skip for the reference run
                continue

            print("  The RTs in {} are being aligned and calibrated".format(colRt[i]))
            # Prepare the modeling: find the shared peptides between the reference and the current run
            idx = (~ref.isna()) & (~df[colRt[i]].isna())
            x = ref[idx]
            y = df[idx][colRt[i]]

            # Build a LOESS model and calibrate RTs
            mod = rLoess(FloatVector(x), FloatVector(y))    # LOESS model based on the shared peptides
            res[colRt[i]] = rPredict(mod, FloatVector(df[colRt[i]]))    # Calibration is applied to whole peptides of the current run

            # Update the reference by merging the current reference and the calibrated run
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

    print("  RT-alignment is finished\n")
    return res
