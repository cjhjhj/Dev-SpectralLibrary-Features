import sys, os, re, utils, numpy as np, pandas as pd
from pyteomics import mzxml
from featureDetection import detectFeatures


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
args = sys.argv
del args[0]
mzXML = args[0]
idTxt = args[1:]
# mzXML = "../Data/FTLD_Batch2_F50.mzXML"
# idTxt = "../Data/ID.txt"

# Feature detection
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
features, ms1ToFeatures = detectFeatures(mzXML, params)
# features = pd.read_pickle("FTLD_Batch2_F50_Features_SN0_Gap3_3ppm.pickle")
# ms1ToFeatures = pd.read_pickle("FTLD_Batch2_F50_MS1_to_Features.pickle")

# Read mzXML file
reader = mzxml.MzXML(mzXML)
ms2ToSurvey = getMs2ToSurvey(reader)

# Read ID.txt files to extract PSM information
print("  Read ID.txt file and feature file")
psms = pd.read_csv(idTxt, skiprows=1, sep=";")  # Note that ID.txt file is delimited by semicolon
psms = psms[["Peptide", "Outfile", "measuredMH", "XCorr"]]
psms = psms.loc[psms["Outfile"].str.contains("FTLD_Batch2_F50")]    # Extract PSMs from FTLD_Batch2_F50.mzXML
psms["precMz"] = np.nan
psms["charge"] = np.nan
psms["featureIndex"] = np.nan
psms["category"] = ""
psms = psms.drop_duplicates()
print("  PSM information has been parsed\n")

# Find the match between features and PSMs
n1, n2 = 0, 0   # n1 = number of PSMs mapped to feature(s), n2 = number of PSMs not mapped to any feature
isoWindow = 1   # Isolation window size for a precursor peak
proton = 1.007276466812
featureName = os.path.basename(mzXML).split(".")[0]
progress = utils.progressBar(psms.shape[0])
for idx, psm in psms.iterrows():
    progress.increment()
    [psmRunName, psmScanNum, _, psmZ, _] = os.path.basename(psm["Outfile"]).split(".")
    psms.loc[idx, "charge"] = psmZ
    if psmRunName == featureName:
        # Extract the precursor m/z
        psmScanNum = int(psmScanNum)
        surveyScanNum = ms2ToSurvey[psmScanNum]
        precMz, _, _ = getPrecursorPeak(reader, psmScanNum, surveyScanNum, params)

        # Assign the feature corresponding to the PSM, based on the precursor m/z and feature's m/z
        if precMz != -1:
            psms.loc[idx, "precMz"] = precMz
            fInd = ms1ToFeatures[((ms1ToFeatures["ms1ScanNumber"] == str(surveyScanNum)) & (ms1ToFeatures["peakMz"] == precMz))]["featureIndex"].values
        else:
            fInd = []
        if len(fInd) > 0:   # Using
            psms.loc[idx, "featureIndex"] = fInd
            n1 += 1
        else:
            n2 += 1
    else:
        continue
# print("There are {} unique PSMs".format(len(psms)))
# print("{} PSMs match to unique features".format(n1))
# print("{} PSMs do not match to any feature".format(n2))

# psms = pd.read_csv("psms_test.txt", sep="\t")
# features = pd.read_pickle("FTLD_Batch2_F50_Features_SN0_Gap3_3ppm.pickle")
# ms1ToFeatures = pd.read_pickle("FTLD_Batch2_F50_MS1_to_Features.pickle")

# Assign RT values to PSMs
res = []
keys = psms["Peptide"] + "_" + psms["charge"].astype(str)
keys = list(set(keys))  # Keep unique "key"s
for key in keys:
    pep, z = key.split("_")
    fInd = psms[(psms["Peptide"] == pep) & (psms["charge"] == int(z))]["featureIndex"].drop_duplicates().dropna()
    if len(fInd) > 0:   # There is/are feature(s) corresponding to the key (i.e., peptide-charge pair)
        # Top-down RT (RT derived from features)
        rtArray = features.loc[fInd]["RT"].values
        intArray = features.loc[fInd]["intensity"].values
    else:
        # Bottom-up RT (RT derived from PSM(s) of the key)
        df = psms[(psms["Peptide"] == pep) & (psms["charge"] == int(z))]
        rtArray = np.array([])
        intArray = np.array([])
        for ind, psm in df.iterrows():
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
res.to_csv("FTLD_Batch2_F50_RT_assignment.txt", sep="\t", index=False)