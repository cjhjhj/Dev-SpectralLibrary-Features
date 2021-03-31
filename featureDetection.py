#!/usr/bin/python

import os, sys, utils, logging, numpy as np, pandas as pd
from pyteomics import mzxml


def detectPeaks(spec, params):
    # Parameters
    if params["data_acquisition_mode"] == "1":
        isCentroided = 1
    elif params["data_acquisition_mode"] == "2":
        isCentroided = 0
    else:
        print("Please set the proper 'data_acquisition_mode' parameter")
        sys.exit("")
    # intensityThreshold = 0  # May come from a parameter file
    intensityThreshold = float(params["min_peak_intensity"])

    # m/z and intensity arrays from a spectrum object
    mzArray = spec["m/z array"]
    intensityArray = spec["intensity array"]
    nPeaks = len(mzArray)
    newMzArray = np.array([])
    newIntensityArray = np.array([])

    # Detect peaks (i.e. centroidization of MS1 spectrum)
    if isCentroided == 0:  # i.e. Profile mode MS1
        for i in range(2, nPeaks - 2):
            if intensityArray[i] > 0:
                # Consider 2 points before and after the point of interest x, i.e. 5 point window
                b2, b1, x, a1, a2 = intensityArray[(i - 2):(i + 3)]
                if x >= intensityThreshold:
                    if isMax(b2, b1, x, a1, a2):
                        # If x is the local maximum in a 5-point window, lower and upper bounds for a peak will be explored
                        # Refer Figure 1a and b in the paper, Cox and Mann, Nature Biotech. 2008; 26: 1367-22
                        minInd = findMinPeakIndex(i, intensityArray)
                        maxInd = findMaxPeakIndex(i, intensityArray)
                        if (maxInd - minInd) > 2:
                            newMz, newIntensity = findPeakCenter(minInd, i, maxInd, mzArray, intensityArray)
                            newMzArray = np.append(newMzArray, newMz)
                            newIntensityArray = np.append(newIntensityArray, newIntensity)
    else:  # i.e. Centroid mode MS1
        idx = np.where(intensityArray >= intensityThreshold)[0]
        newMzArray = mzArray[idx]
        newIntensityArray = intensityArray[idx]

    # Update "spec" object
    spec["m/z array"] = newMzArray
    spec["intensity array"] = newIntensityArray
    return spec


def isMax(b2, b1, x, a1, a2):
    if x > b1 and x > a1:
        return True
    if x > b2 and x == b1 and x > a1:
        return True
    if x > b1 and x == a1 and x > a2:
        return True
    return False


def findMinPeakIndex(ind, array):
    while ind > 0 and array[ind] != 0 and array[ind - 1] <= array[ind]:
        ind -= 1
    return ind + 1


def findMaxPeakIndex(ind, array):
    count = len(array)
    while ind < count and array[ind] != 0 and array[ind + 1] <= array[ind]:
        ind += 1
    return ind - 1


def findPeakCenter(minInd, centerInd, maxInd, mz, intensity):
    # Find the center of a peak composed of five data points
    centerMz = 0
    centerIntensity = 0
    for i in range(minInd, maxInd + 1):
        if intensity[i] >= centerIntensity:
            centerIntensity = intensity[i]  # Take the maximum intensity as a peak intensity

    # There"s a plateau, bu others are zeros
    if minInd == maxInd:
        centerMz = mz[maxInd]
        return centerMz, centerIntensity

    # Right-angled triangle-shaped peak
    if minInd == centerInd:
        centerMz = estimate2(mz[centerInd], mz[centerInd + 1], intensity[centerInd], intensity[centerInd + 1])
        return centerMz, centerIntensity

    # Left-angled triangle-shaped peak
    if maxInd == centerInd:
        centerMz = estimate2(mz[centerInd - 1], mz[centerInd], intensity[centerInd - 1], intensity[centerInd])
        return centerMz, centerIntensity

    # Typical bell(triangle)-shaped peak
    centerMz = estimate3(mz[centerInd - 1], mz[centerInd], mz[centerInd + 1], intensity[centerInd - 1],
                         intensity[centerInd], intensity[centerInd + 1])
    return centerMz, centerIntensity


def estimate2(m1, m2, i1, i2):
    centerVal = (m1 * i1 + m2 * i2) / (i1 + i2)  # Intensity-weighted average of m/z
    return centerVal


def estimate3(m1, m2, m3, i1, i2, i3):
    l1 = np.log(i1)
    l2 = np.log(i2)
    l3 = np.log(i3)
    centerVal = 0.5 * ((l1 - l2) * (m3 ** 2 - m1 ** 2) - (l1 - l3) * (m2 ** 2 - m1 ** 2)) / (
            (l1 - l2) * (m3 - m1) - (l1 - l3) * (m2 - m1))
    return centerVal


def findPeakMatch(array, value, tolerance):
    mzArray = np.asarray(array["m/z array"])
    ind = np.abs(mzArray - value).argmin()  # Index of the closest element to "value" in the array
    difference = np.abs(mzArray[ind] - value) / value * 1e6
    if difference <= tolerance:
        return True, ind
    else:
        return False, 0


def reduceMS1(spec, noise, array):
    # Input
    # spec: spectrum object read by pyteomics
    # array: index of "m/z array" (and "intensity array") to be retained
    array = array.astype(int)

    # Noise level estimation
    ind = np.setdiff1d(range(len(spec["m/z array"])), array)
    if ind.size == 0:
        noiseLevel = 500    # When noiseLevel cannot be estimated, a default noiseLevel is set to 500
    else:
        noiseLevel = np.percentile(spec["intensity array"][ind], 25)    # 25 percentile = 1st quartile = median of bottom 50%
    # noiseLevel = np.percentile([spec["intensity array"][i] for i in ind], 25)   # 25 percentile = 1st quartile = median of bottom 50%
    noise[spec["num"]] = noiseLevel

    # Reduce m/z array and intensity array of spec and replace
    spec["m/z array"] = spec["m/z array"][array]
    spec["intensity array"] = spec["intensity array"][array]
    # rMz = [spec["m/z array"][i] for i in array]
    # rIntensity = [spec["intensity array"][i] for i in array]
    # spec["m/z array"] = rMz
    # spec["intensity array"] = rIntensity

    return spec, noise


def getClosest(spec, mz, tol):
    ind = np.argmin(abs(spec["m/z array"] - mz))
    diff = abs(mz - spec["m/z array"][ind]) / mz * 1e6
    if diff < tol:
        return 1, ind
    else:
        return 0, None


def dechargeFeatures(features):
    # np.seterr(all='raise')  # Necessary when detecting RuntimeError caused by numpy overflow, underflow, etc
    # features = features of np.recarray format
    features = np.sort(features, order = "intensity")[::-1]  # Sort features in descending order of intensity
    delC = 1.00335  # Mass difference between 13C and 12C
    tolPpm = 10  # Tolerance for decharging
    maxCharge = 4
    for i in range(features.shape[0]):
        mz = features["mz"][i]
        intensity = features["intensity"][i]
        lL = mz - mz * tolPpm / 1e6
        uL = delC + mz + mz * tolPpm / 1e6
        scan = features["MS1"][i]
        ind = np.where((features["MS1"] > (scan - 50)) & (features["MS1"] < (scan + 50)))[0]
        charge = 0
        for j in ind:
            if j == i:
                continue
            else:
                mz_j = features["mz"][j]
                intensity_j = features["intensity"][j]
                # A presumable isotopic peak intensity should be greater than 20% of feature intensity (to prevent the inclusion of small/noisy peaks)
                # and smaller than 500% (i.e. inverse of 20%) of feature intensity (to prevent that a monoisotopic peak is merged for decharging of a small/noisy peak)
                # if lL <= mz_j < uL and (0.2 * intensity) <= intensity_j < (5 * intensity):
                if lL <= mz_j < uL:
                    # Check the overlap (in RT-dimension) between two features
                    min1 = features["minRT"][i]
                    max1 = features["maxRT"][i]
                    min2 = features["minRT"][j]
                    max2 = features["maxRT"][j]
                    minRange = min(max1 - min1, max2 - min2)
                    rtOverlap = findRtOverlap(min1, max1, min2, max2)
                    if minRange == 0 or (rtOverlap / minRange) < 0.5:
                        continue

                    # Look for potential isotopic peak and decharge
                    diff = np.around(1 / abs(mz - mz_j)).astype(int)
                    if diff == 0 or diff > maxCharge:
                        continue
                    dev = abs(abs(mz - mz_j) - (delC / diff))
                    if dev > (mz * tolPpm / 1e6):
                        continue
                    else:
                        charge = diff
                        if intensity > intensity_j:
                            features["isotope"][j] = 1
                        break
                else:
                    continue
        features["z"][i] = charge

    # Remove the features from isotopic peaks
    ind = np.where(features["isotope"] == 1)[0]
    features = np.delete(features, ind)

    return features


def findRtOverlap(min1, max1, min2, max2):
    if min1 < min2:
        # Examples
        # |----------------|    range of input1
        #       |----|          range of input2
        #          |----------| range of input2
        if max1 < max2:
            overlap = max1 - min2
        else:
            overlap = max2 - min2
    else:
        # Examples
        #      |--------|       range of input1
        # |----------|          range of input2
        # |-----------------|   range of input2
        if max1 < max2:
            overlap = max1 - min1
        else:
            overlap = max2 - min1

    return overlap


def detectFeatures(inputFile, params):
    ##############
    # Parameters #
    ##############
    firstScan = int(params["first_scan_extraction"])
    lastScan = int(params["last_scan_extraction"])
    gap = int(params["skipping_scans"])
    scanWindow = gap + 1
    matchPpm = float(params["mass_tolerance_peak_matching"])

    ##################
    # Initialization #
    ##################
    reader = mzxml.read(inputFile)
    f = []  # Feature array
    nFeatures = -1
    cache = []
    noise = {}  # Empty dictionary for noise level information
    oldMinInd = -1
    oldMaxInd = -1

    ############################
    # Get MS1 scan information #
    ############################
    ms = []
    with reader:
        msCount = 0
        # filename = os.path.basename(inputFile)
        # print("  Extraction of MS1 spectra from %s" % filename)
        for spec in reader:
            msLevel = int(spec["msLevel"])
            scanNum = int(spec["num"])
            if msLevel == 1 and firstScan <= scanNum <= lastScan:
                ms.append(spec)
                msCount += 1
            elif scanNum > lastScan:
                break
        # print("  Done")

    ################################
    # Feature (3D-peak) generation #
    ################################
    filename = os.path.basename(inputFile)
    print("  Feature detection from %s" % filename)
    logging.info("  Feature detection from " + filename)
    progress = utils.progressBar(msCount)
    for i in range(msCount):
        progress.increment()
        minInd = max(0, i - gap - 1)
        maxInd = min(msCount - 1, i + gap + 1)
        if i == 0:
            for j in range(maxInd + 1):
                spec = detectPeaks(ms[j], params)
                spec["index"] = j
                cache.append(spec)
        else:
            for j in range(oldMinInd, minInd):
                cache.pop(0)  # Remove the first element in cache
            for j in range(oldMaxInd + 1, maxInd + 1):
                spec = detectPeaks(ms[j], params)
                spec["index"] = j
                cache.append(spec)

        ##################
        # Reduction step #
        ##################
        p = cache[i - minInd]
        pCount = len(p["m/z array"])
        valids = np.array([])
        count = 0
        for j in range(pCount):
            cm = p["m/z array"][j]
            match = 0
            nTry = 0
            # Backward search
            for k in range(i - 1, minInd - 1, -1):
                q = cache[k - minInd]
                if q["m/z array"].size == 0:
                    continue
                else:
                    match, ind = getClosest(q, cm, matchPpm)
                if match == 1:
                    break
                nTry += 1
                if nTry > scanWindow:
                    break
            if match == 0:  # Forward search
                nTry = 0
                for k in range(i + 1, maxInd + 1):
                    q = cache[k - minInd]
                    if q["m/z array"].size == 0:
                        continue
                    else:
                        match, ind = getClosest(q, cm, matchPpm)
                    if match == 1:
                        break
                    nTry += 1
                    if nTry > scanWindow:
                        break
            if match == 1:
                valids = np.append(valids, j)

        # Peak reduction and noise-level estimation
        p, noise = reduceMS1(p, noise, valids)

        #####################
        # Peak merging step #
        #####################
        cache[i - minInd] = p
        pCount = len(p["m/z array"])
        for j in range(pCount):
            cm = p["m/z array"][j]
            match = 0
            nTry = 0
            matchedPeakInd = []
            # Backward search
            for k in range(i - 1, minInd - 1, -1):
                q = cache[k - minInd]
                if q["m/z array"].size == 0:
                    continue
                else:
                    matchIndicator, ind = getClosest(q, cm, matchPpm)
                    # $matchIndicator = 1 means that the j-th (reduced) peak in the i-th scan
                    # can form a 3D-peak with $ind-th (reduced) peak in the previous scan (%q)
                    if matchIndicator == 1:
                        matchedPeakInd.append(q["featureIndex"][ind])
                        match = 1
            if match == 1:
                matchedPeakInd = list(set(matchedPeakInd))  # Make the list unique
                fInd = None
                if len(matchedPeakInd) > 1:  # There are multiple matches to the peaks in previous scans
                    fInd = min(matchedPeakInd)
                    for m in matchedPeakInd:
                        # Merge to the lowest indexed feature and remove the "merged" features
                        if m != fInd:
                            f[fInd]["mz"].extend(f[m]["mz"])
                            f[fInd]["intensity"].extend(f[m]["intensity"])
                            f[fInd]["num"].extend(f[m]["num"])
                            f[fInd]["rt"].extend(f[m]["rt"])
                            f[fInd]["index"].extend(f[m]["index"])

                            # Revise cache array
                            for s in f[m]["index"]:
                                for t in range(len(cache)):
                                    if cache[t]["index"] == s:
                                        for u in range(len(cache[t]["featureIndex"])):
                                            if cache[t]["featureIndex"][u] == m:
                                                cache[t]["featureIndex"][u] = fInd
                            f[m] = None  # Keep the size of feature array
                else:
                    fInd = matchedPeakInd[0]
                if "featureIndex" in cache[i - minInd]:
                    cache[i - minInd]["featureIndex"].append(fInd)
                else:
                    cache[i - minInd]["featureIndex"] = [fInd]
                f[fInd]["mz"].append(p["m/z array"][j])
                f[fInd]["intensity"].append(p["intensity array"][j])
                f[fInd]["num"].append(p["num"])
                f[fInd]["rt"].append(p["retentionTime"])
                f[fInd]["index"].append(p["index"])

            if match != 1:
                if i < msCount:
                    nFeatures += 1
                    if "featureIndex" in cache[i - minInd]:
                        cache[i - minInd]["featureIndex"].append(nFeatures)
                    else:
                        cache[i - minInd]["featureIndex"] = [nFeatures]
                    f.append({"mz": [p["m/z array"][j]],
                              "intensity": [p["intensity array"][j]],
                              "num": [p["num"]],
                              "rt": [p["retentionTime"]],
                              "index": [i]})

        oldMinInd = minInd
        oldMaxInd = maxInd

    # Remove empty features
    f = [i for i in f if i is not None]

    #################################
    # Filtering features (3D-peaks) #
    #################################
    # A feature may contain multiple peaks from one scan
    # In this case, one with the largest intensity is chosen
    gMinRt, gMaxRt = 0, 0  # Global minimum and maximum RT over all features
    for i in range(len(f)):
        if len(f[i]["num"]) != len(list(set(f[i]["num"]))):
            temp = {}
            for j in range(len(f[i]["num"])):
                if f[i]["num"][j] in temp:
                    currIntensity = f[i]["intensity"][j]
                    if currIntensity > temp[f[i]["num"][j]]["intensity"]:
                        temp[f[i]["num"][j]]["intensity"] = currIntensity
                        temp[f[i]["num"][j]]["index"] = j
                else:
                    temp[f[i]["num"][j]] = {}
                    temp[f[i]["num"][j]]["intensity"] = f[i]["intensity"][j]
                    temp[f[i]["num"][j]]["index"] = j
            uInd = []
            for key in sorted(temp.keys()):
                uInd.append(temp[key]["index"])
            f[i]["mz"] = [f[i]["mz"][u] for u in uInd]
            f[i]["intensity"] = [f[i]["intensity"][u] for u in uInd]
            f[i]["num"] = [f[i]["num"][u] for u in uInd]
            f[i]["rt"] = [f[i]["rt"][u] for u in uInd]
            f[i]["index"] = [f[i]["index"][u] for u in uInd]

        if i == 0:
            gMinRt = min(f[i]["rt"])
            gMaxRt = max(f[i]["rt"])
        else:
            if min(f[i]["rt"]) < gMinRt:
                gMinRt = min(f[i]["rt"])
            if max(f[i]["rt"]) > gMaxRt:
                gMaxRt = max(f[i]["rt"])

    if gMaxRt.unit_info == "minute":
        gMaxRt = gMaxRt * 60
        gMinRt = gMinRt * 60

    ###################################
    # Organization of output features #
    ###################################
    n = 0
    ms1ToFeatures = {}
    for i in range(len(f)):
        # 1. mz: mean m/z of a feauture = weighted average of m/z and intensity
        mz = np.sum(np.multiply(f[i]["mz"], f[i]["intensity"])) / np.sum(f[i]["intensity"])

        # 2. intensity: intensity of a feature (maximum intensity among the peaks consist of the feature)
        intensity = max(f[i]["intensity"])

        # 3. z: charge of the feature, set to 1 now, but modified later
        z = 1
        isotope = 0  # Will be used later

        # 4. RT: RT of the representative peak (i.e. strongest peak) of a feature
        ind = np.argmax(f[i]["intensity"])
        rt = f[i]["rt"][ind]

        # 5. minRT and maxRT
        minRt = min(f[i]["rt"])
        maxRt = max(f[i]["rt"])

        # Conversion of RT to the unit of second
        if rt.unit_info == "minute":
            rt = rt * 60  # Convert to the unit of second
            minRt = minRt * 60
            maxRt = maxRt * 60

        # 6. MS1 scan number of the representative peak of a feature
        ms1 = f[i]["num"][ind]

        # 7. minMS1 and maxMS1
        minMs1 = min(list(map(int, f[i]["num"])))
        maxMs1 = max(list(map(int, f[i]["num"])))

        # 8. SNratio (signal-to-noise ratio of the feature)
        if ms1 in noise:
            noiseLevel = noise[ms1]
        else:
            noiseLevel = 500
        snRatio = intensity / noiseLevel
        featureIntensityThreshold = noiseLevel * float(params["signal_noise_ratio"])

        if intensity >= featureIntensityThreshold:
            # 9. Percentage of true feature
            pctTF = (maxRt - minRt) / (gMaxRt - gMinRt) * 100
            # Organize features in a structured numpy array form
            if n == 0:
                features = np.array([(mz, intensity, z, rt, minRt, maxRt, ms1, minMs1, maxMs1, snRatio, pctTF, isotope)],
                                    dtype="f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8")
                n += 1
            else:
                features = np.append(features,
                                     np.array([(mz, intensity, z, rt, minRt, maxRt, ms1, minMs1, maxMs1, snRatio, pctTF, isotope)],
                                              dtype=features.dtype))
            for j in range(len(f[i]["num"])):
                num = f[i]["num"][j]
                if num not in ms1ToFeatures:
                    ms1ToFeatures[num] = {"mz": [f[i]["mz"][j]],
                                          "intensity": [f[i]["intensity"][j]]}
                else:
                    ms1ToFeatures[num]["mz"].append(f[i]["mz"][j])
                    ms1ToFeatures[num]["intensity"].append(f[i]["intensity"][j])
        else:
            continue

    features.dtype.names = ("mz", "intensity", "z", "RT", "minRT", "maxRT", "MS1", "minMS1", "maxMS1", "SNratio", "PercentageTF", "isotope")

    ##########################
    # Decharging of features #
    ##########################
    features = dechargeFeatures(features)

    ############################################
    # Convert the features to pandas dataframe #
    # Write features to a file                 #
    ############################################
    df = pd.DataFrame(features)
    df = df.drop(columns = ["isotope"])    # "isotope" column was internally used, and no need to be transferred

    # Create a subdirectory and save features to a file
    baseFilename = os.path.splitext(os.path.basename(filename))[0]  # i.e. filename without extension
    featureDirectory = os.path.join(os.getcwd(), baseFilename)
    if not os.path.exists(featureDirectory):
        os.mkdir(featureDirectory)
    featureFilename = baseFilename + "_" + str(params["skipping_scans"]) + ".feature"
    featureFilename = os.path.join(featureDirectory, featureFilename)
    df.to_csv(featureFilename, index = False, sep = "\t")

    return df  # Pandas DataFrame
