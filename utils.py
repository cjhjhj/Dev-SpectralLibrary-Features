import re, sys, os, logging, numpy as np, pandas as pd


def getParams(paramFile):
    parameters = dict()
    with open(paramFile, 'r') as file:
        for line in file:
            if re.search(r'^#', line) or re.search(r'^\s', line):
                continue
            line = re.sub(r'#.*', '', line)  # Remove comments (start from '#')
            line = re.sub(r'\s*', '', line)  # Remove all whitespaces

            # Exception for "feature_files" parameter
            if "feature_files" in parameters and line.endswith(".feature"):
                parameters["feature_files"].append(line)
            elif "library" in parameters and line.endswith(".db"):
                parameters["library"].append(line)
            else:
                key = line.split('=')[0]
                val = line.split('=')[1]
                if key == "feature_files" or key == "library":
                    parameters[key] = [val]
                else:
                    parameters[key] = val
    return parameters


def readFeatures(featureFile):
    with open(featureFile, 'r') as file:
        features = file.readlines()
    return features


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


def summarizeFeatures(full, params):
    # Input arguments
    # 1. full: numpy recarray of fully-aligned features
    # 2. params: dictionary of parameters

    ####################################
    # Summarize fully-aligned features #
    ####################################
    mzCols = [col for col in full.dtype.names if col.lower().endswith("_mz")]
    rtCols = [col for col in full.dtype.names if col.lower().endswith("_rt")]
    intensityCols = [col for col in full.dtype.names if col.lower().endswith("_intensity")]
    chargeCols = [col for col in full.dtype.names if col.lower().endswith("_z")]
    minRtCols = [col for col in full.dtype.names if col.lower().endswith("_minrt")]
    maxRtCols = [col for col in full.dtype.names if col.lower().endswith("_maxrt")]
    snCols = [col for col in full.dtype.names if col.lower().endswith("_snratio")]

    df = pd.DataFrame.from_records(full)
    res = pd.DataFrame()
    res["feature_m/z"] = df[mzCols].mean(axis=1)
    res["feature_RT"] = df[rtCols].mean(axis=1) / 60
    res["feature_intensity"] = df[intensityCols].mean(axis=1)
    res["feature_SNratio"] = df[snCols].mean(axis=1)
    res["feature_width"] = pd.DataFrame((df[maxRtCols].values - df[minRtCols].values) / 60).mean(axis=1)
    # Handling charges
    res["feature_z"] = df[chargeCols].mode(axis=1)[0]
    res["feature_z"][res["feature_z"] == 0] = 1
    res["feature_z"] = res["feature_z"].astype(int)
    if params["mode"] == "-1":
        res["feature_ion"] = "[M-" + res["feature_z"].astype(str) + "H]" + res["feature_z"].astype(str) + "-"
        res["feature_ion"] = res["feature_ion"].replace("[M-1H]1-", "[M-H]-")
    elif params["mode"] == "1":
        res["feature_ion"] = "[M+" + res["feature_z"].astype(str) + "H]" + res["feature_z"].astype(str) + "+"
        res["feature_ion"] = res["feature_ion"].replace("[M+1H]1+", "[M+H]+")
    for intensityCol in intensityCols:
        res[intensityCol] = df[intensityCol]
    colNames = ["feature_ion", "feature_z", "feature_m/z", "feature_RT", "feature_width",
                "feature_SNratio"] + intensityCols
    res = res[colNames]
    return res


'''
def generateSummarizedFeatureFile(nFeatures, full, ms2, params):
    filePath = os.path.join(os.getcwd(), "align_" + params["output_name"])

    #############################################################
    # Summarize fully-aligned features and write them to a file #
    #############################################################
    # This file contains "summarized" information of fully-aligned features
    # e.g. mean m/z, mean intensity, mean RT of fully-aligned features and so on
    #      width and SNratio are from the reference run
    mzCols = [col for col in full.dtype.names if col.lower().endswith("_mz")]
    rtCols = [col for col in full.dtype.names if col.lower().endswith("_rt")]
    intensityCols = [col for col in full.dtype.names if col.lower().endswith("_intensity")]
    chargeCols = [col for col in full.dtype.names if col.lower().endswith("_z")]
    minRtCols = [col for col in full.dtype.names if col.lower().endswith("_minrt")]
    maxRtCols = [col for col in full.dtype.names if col.lower().endswith("_maxrt")]
    snCols = [col for col in full.dtype.names if col.lower().endswith("_snratio")]

    df = pd.DataFrame.from_records(full)
    res = pd.DataFrame()
    res["feature_m/z"] = df[mzCols].mean(axis=1)
    res["feature_RT"] = df[rtCols].mean(axis=1)
    res["feature_intensity"] = df[intensityCols].mean(axis=1)
    res["feature_SNratio"] = df[snCols].mean(axis=1)
    res["feature_width"] = pd.DataFrame((df[maxRtCols].values - df[minRtCols].values) / 60).mean(axis=1)
    # Handling charges
    res["feature_z"] = df[chargeCols].mode(axis=1)[0]
    res["feature_z"][res["feature_z"] == 0] = 1
    res["feature_z"] = res["feature_z"].astype(int)
    if params["mode"] == "-1":
        res["feature_ion"] = "[M-" + res["feature_z"].astype(str) + "H]" + res["feature_z"].astype(str) + "-"
        res["feature_ion"] = res["feature_ion"].replace("[M-1H]1-", "[M-H]-")
    elif params["mode"] == "1":
        res["feature_ion"] = "[M+" + res["feature_z"].astype(str) + "H]" + res["feature_z"].astype(str) + "+"
        res["feature_ion"] = res["feature_ion"].replace("[M+1H]1+", "[M+H]+")
    for intensityCol in intensityCols:
        res[intensityCol] = df[intensityCol]

    # Add the mean m/z of feature and its charge state to the beginning of MS2 spectrum (similar to .dta file)
    for i in range(nFeatures):
        if ms2[i] is not None:
            ms2[i]["mz"] = np.insert(ms2[i]["mz"], 0, res["feature_m/z"].iloc[i])
            ms2[i]["intensity"] = np.insert(ms2[i]["intensity"], 0, res["feature_z"].iloc[i])
    res["MS2"] = ms2

    # Write the summarized fully-aligned features to a file
    res = res.sort_values(by="feature_m/z", ignore_index=True)
    res[
        "feature_num"] = res.index + 1  # Update "feature_num" according to the ascending order of "feature_m/z" (as sorted)
    outColumns = ["feature_num", "feature_ion", "feature_z", "feature_m/z", "feature_RT",
                  "feature_width", "feature_SNratio"] + intensityCols
    resOut = res[outColumns].copy()
    resOut["feature_RT"] = resOut["feature_RT"] / 60  # Change the unit to minute
    fullName = os.path.join(filePath, params["output_name"] + "_summarized_fully_aligned.feature")
    resOut.to_csv(fullName, index=False, sep="\t")

    # Write MS2 spectra to files
    ms2Path = os.path.join(filePath, "MS2")
    if not os.path.exists(ms2Path):
        os.mkdir(ms2Path)
    for i in range(res.shape[0]):
        if res["MS2"].iloc[i] is not None:
            fileName = os.path.join(ms2Path, "f" + str(i + 1) + ".MS2")
            dfMS2 = pd.DataFrame.from_dict(res["MS2"].iloc[i])
            dfMS2.to_csv(fileName, index=False, header=False, sep="\t")

    # Save fully-aligned features with their MS2 spectra (i.e. res) for debugging purpose
    # When the pipeline gets mature, this part needs to be removed
    pickle.dump(res, open(os.path.join(filePath, ".fully_aligned_feature.pickle"), "wb"))  # Make the file be hidden

    return res
'''


def processQuantityData(df, params):
    print("  Loading-bias summary")
    print("  ====================")
    logging.info("  Loading-bias summary")
    logging.info("  ====================")
    intensityCols = [col for col in df.columns if col.lower().endswith("_intensity")]
    expr = df[intensityCols]

    # Calculation and print-out loading bias information
    sampleNames = expr.columns
    rowMeans = expr.mean(axis=1)
    lexpr = np.log2(expr.div(rowMeans, axis=0))
    nFeatures, nSamples = lexpr.shape
    idx = pd.Series([True] * nFeatures)
    for i in intensityCols:
        idx = idx & (lexpr[i] > lexpr[i].quantile(q=0.1)) & (lexpr[i] < lexpr[i].quantile(q=0.9))
    meanIntensity = 2 ** (lexpr.loc[idx, :].mean(axis=0)) * 100
    sdVal = lexpr.loc[idx, :].std(axis=0)
    sdIntensity = ((2 ** sdVal - 1) + (1 - 2 ** (-sdVal))) / 2 * 100
    semIntensity = sdIntensity / np.sqrt(len(idx))
    print("  Sample_name\tMean[%]\tSD[%]\tSEM[%]\t#features")
    logging.info("  Sample_name\tMean[%]\tSD[%]\tSEM[%]\t#features")
    for i in range(expr.shape[1]):
        print("  {}\t{:.2f}\t{:.2f}\t{:.2f}\t{}".format(intensityCols[i].replace("_intensity", ""), meanIntensity[i], sdIntensity[i], semIntensity[i], len(idx)))
        logging.info("  {}\t{:.2f}\t{:.2f}\t{:.2f}\t{}".format(intensityCols[i].replace("_intensity", ""), meanIntensity[i], sdIntensity[i], semIntensity[i], len(idx)))

    # Normalization using trimmed-mean values (loading-bias correction)
    lexpr = np.log2(expr)
    if params["skip_loading_bias_correction"] == 0:
        # Parameters for normalization
        cutoff = np.nanquantile(lexpr.to_numpy(), q=0.1)  # 10% quantile of overall intensities
        # This is the original implementation (as the same as Perl-pipeline), but it may be too stringent
        idx = pd.Series([True] * nFeatures)
        for i in intensityCols:
            idx = idx & (lexpr[lexpr > cutoff][i] > lexpr[lexpr > cutoff][i].quantile(q=0.1)) & (
                        lexpr[lexpr > cutoff][i] > lexpr[lexpr > cutoff][i].quantile(q=0.9))
        meanIntensity = lexpr.loc[idx, :].mean(axis=0)
        meanIntensity = lexpr[lexpr > cutoff].mean(axis=0)
        normFactor = meanIntensity - np.mean(meanIntensity)
        lexpr = lexpr - normFactor

    # Replace missing values (i.e. nan) with the 1/2 of the global minimum intensity
    minIntensity = np.nanmin(lexpr)
    lexpr = lexpr.fillna(minIntensity - 1)  # Half of the global minimum intensity (at log2-scale)
    lexpr += lexpr.isna() * pd.DataFrame(1e-2 * np.random.uniform(0, 1, size=(lexpr.shape)),
                                         columns=lexpr.columns)  # Add small random number for numerical stability
    df[intensityCols] = 2 ** lexpr
    return df

def generateFeatureFile(full, partial, unaligned, params):
    # Input arguments
    # full: fully-aligned features (numpy recarray)
    # partial: partially-aligned features (numpy recarray, or None)
    # unaligned: un-aligned features (list of numpy recarray, or None)
    # params: dictionary of parameters
    filePath = os.path.join(os.getcwd(), "align_" + params["output_name"])
    if not os.path.exists(filePath):
        os.mkdir(filePath)

    # Organize fully-aligned features
    fullName = os.path.join(filePath, params["output_name"] + "_fully_aligned.feature")
    dfFull = pd.DataFrame(full)
    dfFull["meanMz"] = dfFull.filter(regex=(".*mz$")).mean(axis=1)
    dfFull = dfFull.sort_values(by="meanMz", ignore_index=True)  # Features are sorted by mean m/z
    colNames = dfFull.columns.tolist()
    colNames = colNames[-1:] + colNames[:-1]
    dfFull = dfFull[colNames]

    # Organize "summarized" fully-aligned features
    fullName2 = os.path.join(filePath, params["output_name"] + "_summarized_fully_aligned.feature")
    dfFull2 = summarizeFeatures(full, params)
    dfFull2 = dfFull2.sort_values(by="feature_m/z", ignore_index=True)  # Features are sorted by "feature_m/z"
    dfFull2.insert(loc=0, column="feature_num", value=dfFull2.index + 1)

    # Processing of quantity data (missing value imputation, normalization, etc.)
    dfFull2 = processQuantityData(dfFull2, params)
    intensityCols = [col for col in dfFull2.columns if col.lower().endswith("_intensity")]
    for col in intensityCols:
        dfFull[col] = dfFull2[col]  # Replace intensity values of dfFull with those of dfFull2 (preprocessed)

    # This file contains fully-aligned features with run-specific information
    # Since the run-specific information is required for MS2 processing, it should be kept
    dfFull.to_csv(fullName, index=False, sep="\t")
    # This file contains "summarized" fully-aligned features
    # Except intensity, all feature information is summarized over runs
    dfFull2.to_csv(fullName2, index=False, sep="\t")

    ##############################################
    # Write partially-aligned features, if exist #
    ##############################################
    if partial is not None:
        dfPartial = pd.DataFrame(partial)
        if dfPartial.shape[0] > 0:
            partialName = os.path.join(filePath, params["output_name"] + "_partially_aligned.feature")
            dfPartial.to_csv(partialName, index=False, sep="\t")
    else:
        dfPartial = None

    #######################################
    # Write un-aligned features, if exist #
    #######################################
    dfArrayUnaligned = None
    if unaligned is not None:
        dfArrayUnaligned = []
        for un in unaligned:
            unName = [col for col in un.dtype.names if col.endswith('_mz')][0]
            unName = os.path.join(filePath, re.sub("_mz", "", unName) + "_unaligned.feature")
            dfUn = pd.DataFrame(un)
            dfUn.to_csv(unName, index=False, sep="\t")
            dfArrayUnaligned.append(dfUn)
    else:
        dfArrayUnaligned = None

    return dfFull, dfPartial, dfArrayUnaligned
