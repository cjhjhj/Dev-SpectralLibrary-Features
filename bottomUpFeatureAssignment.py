import sys, os, re, utils, numpy as np, pandas as pd
from pyteomics import mzxml
from datetime import datetime


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
    mzArray = reader[str(surveyScanNum)]["m/z array"]
    intArray = reader[str(surveyScanNum)]["intensity array"]
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

    precRt = reader[str(surveyScanNum)]["retentionTime"] * 60   # Unit of second
    return precMz, precIntensity, precRt


# Input
# 1. mzXML file
# 2. ID.txt file containing all identified PSMs
# args = sys.argv
# del args[0]
# mzXML = args[0]
# idTxt = args[1:]
t1 = datetime.now()
mzXML = "../Data/FTLD_Batch2_F50.mzXML"
idTxt = "../Data/ID.txt"
mzXMLBaseName = os.path.basename(mzXML).split(".")[0]

# Parameters
params = {}
params["data_acquisition_mode"] = "1"        # 1 = centroid, 2 = profile for full scan and centroid for MS/MS scan
params["first_scan_extraction"] = "1"        # the first scan used for search
params["last_scan_extraction"] = "100000"    # the last scan used for search
params["isolation_window"] = "1"             # isolation window size 1= +/-0.5
params["mass_correction"] = "0"              # 0 = no correction, 1 = MS1-based
params["signal_noise_ratio"] = "0"           # fold of the minimum signal noise ratio
params["max_percentage_RT_range"] = "100"    # threshold maximum percentage of the range of retention time of a peak
params["min_peak_intensity"] = "10000"       # threshold of a peak intensity
params["skipping_scans"] = "3"               # number of skipping scans during 3D formation
params["mass_tolerance_peak_matching"] = "3" # mass tolerance for peak matching during 3D formation

# Read mzXML file
reader = mzxml.MzXML(mzXML)
ms2ToSurvey = getMs2ToSurvey(reader)

# Read ID.txt files to extract PSM information
print("  Read ID.txt file and feature file")
psms = pd.read_csv(idTxt, skiprows=1, sep=";")  # Note that ID.txt file is delimited by semicolon
psms = psms[["Peptide", "Outfile", "measuredMH", "XCorr"]]
psms = psms.loc[psms["Outfile"].str.contains(mzXMLBaseName)]    # Extract PSMs from FTLD_Batch2_F50.mzXML
psms["charge"] = [outfile.split("/")[-1].split(".")[-2] for outfile in psms["Outfile"]]
psms = psms.drop_duplicates()
print("  PSM information has been parsed\n")

# Unique key is peptide-charge pair
keys = psms["Peptide"] + "_" + psms["charge"]
keys = list(set(keys))
res = []
progress = utils.progressBar(len(keys))
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
    rt = sum(rtArray * intArray) / sum(intArray) / 60
    rtStd = np.std(rtArray / 60)
    res.append([key, rt, rtStd])

res = pd.DataFrame(res, columns=["peptide_charge", "RT", "stdRT"])
outputFile = mzXMLBaseName + "_BottomUp_RT_assignment.txt"
res.to_csv(outputFile, sep="\t", index=False)
t2 = (datetime.now() - t1).total_seconds()
print(t2)