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
    print("  Read a mzXML file: to find survey scans of MS2 scans")
    res = {}
    for spec in reader:
        if spec["msLevel"] == 1:
            survey = int(spec["num"])
        elif spec["msLevel"] == 2:
            res[int(spec["num"])] = survey
    print("  Done ...\n")
    return res


def getPrecursorPeak(reader, psmScanNumber, surveyScanNumber, params):
    isoWindow = float(params["isolation_window"])

    # Find the precursor peak of the PSM (the strongest peak within the isolation window)
    nominalPrecMz = float(re.search("([0-9.]+)\\@", reader[str(psmScanNumber)]["filterLine"]).group(1))
    mzArray = reader[str(surveyScanNumber)]["m/z array"]
    intArray = reader[str(surveyScanNumber)]["intensity array"]
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

    precRt = reader[str(surveyScanNumber)]["retentionTime"]   # Unit of minute
    return precMz, precIntensity, precRt


#################
# Main function #
#################
def getRT(mzXML, idTxt):
    # Input
    # 1. mzXML file
    # 2. ID.txt file containing all identified PSMs
    mzXMLBaseName = os.path.basename(mzXML).split(".")[0]

    # Parameter(s)
    params = {"isolation_window": 1}    # isolation window size 1= +/-0.5

    # Read mzXML file
    reader = mzxml.MzXML(mzXML)
    ms2ToSurvey = getMs2ToSurvey(reader)

    # Read ID.txt files to extract PSM information
    print("  Read ID.txt file: to extract PSM information")
    psms = pd.read_csv(idTxt, skiprows=1, sep=";")  # Note that ID.txt file is delimited by semicolon
    psms = psms[["Peptide", "Outfile", "XCorr"]]
    psms = psms.loc[psms["Outfile"].str.contains(mzXMLBaseName)]    # Extract PSMs from FTLD_Batch2_F50.mzXML
    psms["charge"] = [outfile.split("/")[-1].split(".")[-2] for outfile in psms["Outfile"]]
    psms = psms.drop_duplicates()
    print("  Done ...\n")

    # Unique key is peptide-charge pair
    print("  RT of every identified peptide is being inferred and assigned")
    keys = psms["Peptide"] + "_" + psms["charge"]
    keys = list(set(keys))
    res = []
    progress = progressBar(len(keys))
    for key in keys:
        progress.increment()
        pep, z = key.split("_")
        rtArray = np.array([])
        intArray = np.array([])
        for _, psm in psms[(psms["Peptide"] == pep) & (psms["charge"] == z)].iterrows():
            [_, psmScanNum, _, _, _] = os.path.basename(psm["Outfile"]).split(".")
            psmScanNum = int(psmScanNum)
            surveyScanNum = ms2ToSurvey[psmScanNum]
            _, precIntensity, precRt = getPrecursorPeak(reader, int(psmScanNum), surveyScanNum, params)
            rtArray = np.append(rtArray, precRt)
            intArray = np.append(intArray, precIntensity)
        rt = sum(rtArray * intArray) / sum(intArray)   # Unit of minute
        # rtStd = np.std(rtArray)
        res.append([pep, z, rt])
    print("  Done ...\n")

    res = pd.DataFrame(res, columns=["peptide", "charge", "RT"])
    return res