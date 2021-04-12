import sys, os, re, utils, numpy as np, pandas as pd
from pyteomics import mzxml
from featureDetection import detectFeatures

# Input
# 1. mzXML file
# 2. ID.txt file containing all identified PSMs
# args = sys.argv
# del args[0]
# featureFile = args[0]
# idTxt = args[1:]
mzXML = "../Data/FTLD_Batch2_F50.mzXML"
idTxt = "../Data/ID.txt"

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
# features, ms1ToFeatures = detectFeatures(mzXML, params)
features = pd.read_pickle("FTLD_Batch2_F50_Features_SN0_Gap3_3ppm.pickle")
ms1ToFeatures = pd.read_pickle("FTLD_Batch2_F50_MS1_to_Features.pickle")
print()

# Read mzXML file
print("  Read a mzXML file: to find survey scans of MS2 scans")
ms2ToSurvey = {}
reader = mzxml.MzXML(mzXML)
for spec in reader:
    if spec["msLevel"] == 1:
        survey = int(spec["num"])
    elif spec["msLevel"] == 2:
        ms2ToSurvey[int(spec["num"])] = survey
print("  Done ...\n")

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
proton = 1.007276466812
featureName = "FTLD_Batch2_F50"
# featureName = os.path.basename(featureFile).split(".")[0]
n1, n2, n3 = 0, 0, 0
progress = utils.progressBar(psms.shape[0])
for idx, psm in psms.iterrows():
    progress.increment()
    [psmRunName, psmScanNum, _, psmZ, _] = os.path.basename(psm["Outfile"]).split(".")
    psms.loc[idx, "charge"] = psmZ
    if psmRunName == featureName:
        # Find the precursor peak of the PSM
        precMz = float(re.search("([0-9.]+)\\@", reader[psmScanNum]["filterLine"]).group(1))
        surveyScanNum = ms2ToSurvey[int(psmScanNum)]
        mzArray = reader[str(surveyScanNum)]["m/z array"]
        intArray = reader[str(surveyScanNum)]["intensity array"]
        ii = np.where((mzArray >= precMz - 0.5) & (mzArray <= precMz + 0.5))
        if len(ii[0] > 0):
            subMzArray = mzArray[ii]
            subIntArray = intArray[ii]
            jj = np.argmax(subIntArray)
            precMz = subMzArray[jj]
            psms.loc[idx, "precMz"] = precMz
            fInd = ms1ToFeatures[((ms1ToFeatures["ms1ScanNumber"] == str(surveyScanNum)) & (ms1ToFeatures["peakMz"] == precMz))]["featureIndex"].values
        else:
            fInd = []
            print(psm)

        if len(fInd) == 1:
            fInd = fInd[0]
            psms.loc[idx, "featureIndex"] = fInd
            n1 += 1
        elif len(fInd) > 1:
            n2 += 1
        else:
            n3 += 1
    else:
        continue

print("There are {} unique PSMs".format(len(psms)))
print("{} PSMs match to unique features".format(n1))
print("{} PSMs match to multiple features".format(n2))
print("{} PSMs do not match to any feature".format(n3))


# For the peptides mapped to features, their mapping pattern/category is determined
subDf = psms[psms["featureIndex"] >=0][["Peptide", "charge", "featureIndex"]] # Sub-DataFrame for the peptides mapped to features
subDf = subDf.drop_duplicates()
pepDict = {}
for i, row in subDf.iterrows():
    pep = row["Peptide"] + "_" + row["charge"]
    fInd = row["featureIndex"]
    if pep in pepDict:
        pepDict[pep].append(fInd)
    else:
        pepDict[pep] = [fInd]

n1, n2, n3 = 0, 0, 0
p1, p2, p3 = [], [], []
f1, f2, f3 = [], [], []
for k, v in pepDict.items():
    if len(v) == 1:
        v = v[0]
        ll = np.where(subDf["featureIndex"] == v)[0]
        if len(ll) == 1:    # One-to-one (peptide-to-feature)
            n1 += 1
            p1.append(k)
            f1.append(v)
        else:   # Many-to-one
            n2 += 1
            p2.append(k)
            f2.append(v)
    else:   # One-to-many
        n3 += 1
        p3.append(k)
        f3.extend(v)

print("One-to-one: peptide-to-feature: {}".format(n1))
print("Many-to-one: peptide-to-feature: {}".format(n2))
print("One-to-many: peptide-to-feature: {}".format(n3))
f1 = list(set(f1))
f2 = list(set(f2))
f3 = list(set(f3))
p1 = list(set(p1))
p2 = list(set(p2))
p3 = list(set(p3))

for k in p1:
    [pep, z] = k.split("_")
    psms.loc[(psms["Peptide"] == pep) & (psms["charge"] == z) & (psms["featureIndex"] >= 0), "category"] = "one-to-one"
for k in p2:
    [pep, z] = k.split("_")
    psms.loc[(psms["Peptide"] == pep) & (psms["charge"] == z) & (psms["featureIndex"] >= 0), "category"] = "many-to-one"
for k in p3:
    [pep, z] = k.split("_")
    psms.loc[(psms["Peptide"] == pep) & (psms["charge"] == z) & (psms["featureIndex"] >= 0), "category"] = "one-to-many"

# psms.loc[(psms["Peptide"].isin(p1)) & (psms["featureIndex"] >= 0), "category"] = "one-to-one"
# psms.loc[(psms["Peptide"].isin(p2)) & (psms["featureIndex"] >= 0), "category"] = "many-to-one"
# psms.loc[(psms["Peptide"].isin(p3)) & (psms["featureIndex"] >= 0), "category"] = "one-to-many"

psms.to_csv("psms.txt", sep="\t", index=False)