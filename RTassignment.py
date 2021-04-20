import sys, os, re, numpy as np, pandas as pd
from pyteomics import mzxml


###############
# Subroutines #
###############
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


def getMs2ToSurvey(reader):
    # Read mzXML file and make a relationship between every MS2 and its survey scan
    res = {}
    for spec in reader:
        if spec["msLevel"] == 1:
            survey = int(spec["num"])
        elif spec["msLevel"] == 2:
            res[int(spec["num"])] = survey
    return res


def getPrecursorPeak(reader, psmScanNumber, surveyScanNumber, params):
    isoWindow = float(params["isolation_window"])

    # Find the precursor peak of the PSM (the strongest peak within the isolation window)
    nominalPrecMz = float(reader[str(psmScanNumber)]["filterLine"].split("@")[0].split(" ")[-1])
    # nominalPrecMz = float(filterLine.search(reader[str(psmScanNumber)]["filterLine"]).group(1))
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


def reformatRtTable(df):
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


#################
# Main function #
#################
def getRT(runs, idTxt):
    # Input
    # 1. mzXML files
    # 2. ID.txt file containing all identified PSMs

    # Parameter(s)
    params = {"isolation_window": 1}  # isolation window size 1= +/-0.5
    # filterLine = re.compile("([0-9.]+)\\@")

    # Read ID.txt files to extract PSM information
    print("  Read ID.txt file: to extract PSM information")
    psms = pd.read_csv(idTxt, skiprows=1, sep=";")  # Note that ID.txt file is delimited by semicolon
    psms = psms[["Peptide", "Outfile", "XCorr"]].drop_duplicates()
    psms["charge"] = [outfile.split("/")[-1].split(".")[-2] for outfile in psms["Outfile"]]
    psms["key"] = psms["Peptide"] + "_" + psms["charge"]
    print("  Done ...\n")

    # RT extraction/assignment for each mzXML file
    res = []
    for run in runs:
        runName = os.path.basename(run).split(".")[0]

        # Read a mzXML file and extract PSMs corresponding to the mzXML file
        reader = mzxml.MzXML(run)
        ms2ToSurvey = getMs2ToSurvey(reader)
        subPsms = psms[psms["Outfile"].str.contains(runName)]

        # Unique key is peptide-charge pair
        print("  RT of every identified peptide in {} is being inferred and assigned".format(runName))
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
                _, precIntensity, precRt = getPrecursorPeak(reader, int(psmScanNum), surveyScanNum, params)
                rtArray = np.append(rtArray, precRt)
                intArray = np.append(intArray, precIntensity)
            rt = sum(rtArray * intArray) / sum(intArray)   # Unit of minute
            res.append([key, runName, rt, len(rtArray)])
        print("  Done ...\n")

    res = pd.DataFrame(res, columns=["key", "run", "RT", "nPSMs"])
    res = reformatRtTable(res)
    return res


    # df = {}
    # for mzXML in mzXMLs:
    #     mzXMLBaseName = os.path.basename(mzXML).split(".")[0]
    #
    #     # Read mzXML file
    #     reader = mzxml.MzXML(mzXML)
    #     ms2ToSurvey = getMs2ToSurvey(reader)
    #
    #     # Assign RT to each PSM
    #     res = []
    #     progress = progressBar(len(keys))
    #     for key in keys:
    #         progress.increment()
    #         rtArray = np.array([])
    #         intArray = np.array([])
    #         outfiles = psms[(psms["key"] == key) & (psms["Outfile"].str.contains(mzXMLBaseName))]["Outfile"].values
    #         if len(outfiles) > 0:
    #             psmScanNums = [int(os.path.basename(outfile).split(".")[1]) for outfile in outfiles]
    #             for psmScanNum in psmScanNums:
    #                 surveyScanNum = ms2ToSurvey[psmScanNum]
    #                 _, precIntensity, precRt = getPrecursorPeak(reader, int(psmScanNum), surveyScanNum, params)
    #                 rtArray = np.append(rtArray, precRt)
    #                 intArray = np.append(intArray, precIntensity)
    #             rt = sum(rtArray * intArray) / sum(intArray)  # Unit of minute
    #             res.append([key, rt])
    #     print("  Done ...\n")
    #     res = pd.DataFrame(res, columns=["key", "RT"])
    #     df[mzXMLBaseName] = res
    #     return df
