#!/usr/bin/python

import sys, os, utils, logging, numpy as np
import rpy2.robjects as ro
from rpy2.robjects.vectors import IntVector, FloatVector
from numpy.lib.recfunctions import merge_arrays, stack_arrays
ro.r['options'](warn=-1)


def calibrateFeatures(ref, comp, params):
    # ref and comp are "ndarray"s with the following column name
    # 'index' = nominal feature index
    # 'mz' = m/z value of a feature
    # 'z' = charge of the feature (0 for undetermined)
    # 'MS1ScanNumber' = representative MS1 scan number of the feature
    # 'minMS1ScanNumber' = minimum MS1 scan number of the feature (spanned)
    # 'maxMS1ScanNumber' = maximum MS1 scan number of the feature (spanned)
    # 'RT' = representative RT (should correspond to MS1ScanNumber)
    # 'minRT' = minimum RT
    # 'maxRT' = maximum RT
    # "intensity" = intensity of the feature
    # 'SN' = signal-to-noise ratio
    # 'PercentageTF' = percentage of true feature (= feature width over RT/scan)

    initMzTol = int(params["tol_initial"])
    sdWidth = int(params["sd_width"])
    ref = np.sort(ref, order="intensity")[::-1]  # Sort features in descending order of intensity
    comp = np.sort(comp, order="intensity")[::-1]

    # Calibration of RT and m/z globally
    print("  Global calibration of features is being performed")
    logging.info("  Global calibration of features is being performed")
    rtShifts, mzShifts = globalCalibration(ref, comp, initMzTol)
    print("    Based on the matched features within %d ppm" % initMzTol)
    print("    The global RT-shift is %.4f second" % np.median(rtShifts))
    print("    The global m/z-shift is %.4f ppm" % np.median(mzShifts))
    print("    RT and m/z of the compared features are calibrated according to the above information")
    logging.info("    Based on the matched features within %d ppm" % initMzTol)
    logging.info("    The global RT-shift is %.4f second" % np.median(rtShifts))
    logging.info("    The global m/z-shift is %.4f ppm" % np.median(mzShifts))
    logging.info("    RT and m/z of the compared features are calibrated according to the above information")
    comp["RT"] = comp["RT"] - np.median(rtShifts)
    comp["mz"] = comp["mz"] / (1 + np.median(mzShifts) / 1e6)

    # Calibration of RT and m/z locally using LOESS (stepwise)
    rLoess = loess()
    rPredict = ro.r("predict")
    print("  Local calibration of features is being performed (through LOESS modeling)")
    print("    RT- and m/z-tolerances will be dynamically estimated over RT- and m/z-range as follows,")
    print("    RT- and m/z-tolerance = %d x dynamically estimated SD of RT- and m/z-shifts" % sdWidth)
    print("    LOESS modeling may take some time. Please be patient ...")
    logging.info("  Local calibration of features is being performed (through LOESS modeling)")
    logging.info("    RT- and m/z-tolerances will be dynamically estimated over RT- and m/z-range as follows,")
    logging.info("    RT- and m/z-tolerance = %d x dynamically estimated SD of RT- and m/z-shifts" % sdWidth)
    logging.info("    LOESS modeling may take some time. Please be patient ...")
    rtSd = np.maximum(1e-3, np.std(rtShifts, ddof=1))
    mzSd = np.maximum(1e-3, np.std(mzShifts, ddof=1))
    ref, comp, rtSd, mzSd = localCalibration(ref, comp, rtSd, mzSd, rLoess, rPredict, params, "RT")
    print("    The 1st round of RT-calibration is done")
    print("      min SD of RT-shifts = %.4f second" % np.amin(rtSd))
    print("      max SD of RT-shifts = %.4f second" % np.amax(rtSd))
    logging.info("    The 1st round of RT-calibration is done")
    logging.info("      min SD of RT-shifts = %.4f second" % np.amin(rtSd))
    logging.info("      max SD of RT-shifts = %.4f second" % np.amax(rtSd))
    ref, comp, rtSd, mzSd = localCalibration(ref, comp, rtSd, mzSd, rLoess, rPredict, params, "RT")
    print("    The 2nd round of RT-calibration is done")
    print("      min SD of RT-shifts = %.4f second" % np.amin(rtSd))
    print("      max SD of RT-shifts = %.4f second" % np.amax(rtSd))
    logging.info("    The 2nd round of RT-calibration is done")
    logging.info("      min SD of RT-shifts = %.4f second" % np.amin(rtSd))
    logging.info("      max SD of RT-shifts = %.4f second" % np.amax(rtSd))
    ref, comp, rtSd, mzSd = localCalibration(ref, comp, rtSd, mzSd, rLoess, rPredict, params, "mz")
    print("    The 1st round of m/z-calibration is done")
    print("      min SD of m/z-shifts = %.4f ppm" % np.amin(mzSd))
    print("      max SD of m/z-shifts = %.4f ppm" % np.amax(mzSd))
    logging.info("    The 1st round of m/z-calibration is done")
    logging.info("      min SD of m/z-shifts = %.4f ppm" % np.amin(mzSd))
    logging.info("      max SD of m/z-shifts = %.4f ppm" % np.amax(mzSd))
    ref, comp, rtSd, mzSd = localCalibration(ref, comp, rtSd, mzSd, rLoess, rPredict, params, "mz")
    print("    The 2nd round of m/z-calibration is done")
    print("      min SD of m/z-shifts = %.4f ppm" % np.amin(mzSd))
    print("      max SD of m/z-shifts = %.4f ppm" % np.amax(mzSd))
    print()
    logging.info("    The 2nd round of m/z-calibration is done")
    logging.info("      min SD of m/z-shifts = %.4f ppm" % np.amin(mzSd))
    logging.info("      max SD of m/z-shifts = %.4f ppm" % np.amax(mzSd))
    logging.info("")
    return comp, rtSd, mzSd  # "comp" is the set of calibrated features


def globalCalibration(ref, comp, mzTol=20):
    nPeaks = round(0.05 * ref.shape[0])  # Number of peaks to be considered for global calibration
    rtShifts = []  # Array for RT-shifts between reference and compared runs
    mzShifts = []  # Array for mz-shifts (ppm)
    i, j = 0, 1
    while j <= nPeaks:
        z = ref["z"][i]  # From the 1st feature of reference run (i.e. strongest feature)
        mz = ref["mz"][i]
        lL = mz - mz * mzTol / 1e6
        uL = mz + mz * mzTol / 1e6
        rt = ref["RT"][i]
        if z == 0:
            # For a reference feature with undetermined charge, consider all possible charges in compared feature
            rowInd = np.where((comp["mz"] >= lL) & (comp["mz"] < uL))[0]
        else:
            rowInd = np.where((comp["mz"] >= lL) & (comp["mz"] < uL) & (comp['z'] == z))[0]

        if len(rowInd) > 0:
            rowInd = rowInd[0]
            rtShifts.append(comp["RT"][rowInd] - rt)
            mzShifts.append((comp["mz"][rowInd] - mz) / mz * 1e6)
            comp = np.delete(comp, rowInd, 0)
            j += 1

        i += 1

    # For more robust calculation, top and bottom 10% values are trimmed
    rtShifts = np.array(rtShifts)
    rtShifts = rtShifts[(rtShifts >= np.percentile(rtShifts, 10)) &
                        (rtShifts <= np.percentile(rtShifts, 90))]
    mzShifts = np.array(mzShifts)
    mzShifts = mzShifts[(mzShifts >= np.percentile(mzShifts, 10)) &
                        (mzShifts <= np.percentile(mzShifts, 90))]
    return rtShifts, mzShifts


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


def localCalibration(ref, comp, rtSd, mzSd, rLoess, rPredict, params, type):
    refInd, compInd = findComparableFeatures(ref, comp, rtSd, mzSd, params)

    # Perform LOESS regression to calibrated RT and m/z
    refRt = ref["RT"][refInd]
    compRt = comp["RT"][compInd]
    refMz = ref["mz"][refInd]
    compMz = comp["mz"][compInd]

    if type == "RT":  # RT-calibration
        if (refRt == compRt).all():
            rtShifts = 1e-6 * np.random.normal(len(refRt))
        else:
            rtShifts = compRt - refRt
        mod = rLoess(FloatVector(compRt), FloatVector(rtShifts))
        compRt = compRt - np.array(mod.rx2("fitted"))  # Calibrated RT based on the model

        # Calculate a new (dynamic) RT-tolerance
        if (refRt == compRt).all():
            rtShifts = 1e-6 * np.random.normal(len(refRt))
        else:
            rtShifts = compRt - refRt  # Calibrated RT-shifts
        ind = np.where((rtShifts >= np.percentile(rtShifts, 10)) &
                       (rtShifts <= np.percentile(rtShifts, 90)))[0]
        modRtSd = rLoess(FloatVector(compRt[ind]), FloatVector(rtShifts[ind] ** 2))
        rtSd = np.sqrt(np.maximum(0, rPredict(modRtSd, FloatVector(ref["RT"]))))

        # Calculate a new (dynamic) m/z-tolerance
        if (refMz == compMz).all():
            mzShifts = 1e-6 * np.random.normal(len(refMz))
        else:
            mzShifts = (compMz - refMz) / refMz * 1e6
        # Sometimes, the variation of mzShifts cannot be captured when trimming is applied
        # so, the trimming is not used for mzShifts
        modMzSd = rLoess(FloatVector(compMz), FloatVector(mzShifts ** 2), 1, "aicc", "gaussian")
        mzSd = np.sqrt(np.maximum(0, rPredict(modMzSd, FloatVector(ref["mz"]))))

        # Calibration of the entire comp["RT"]
        comp["RT"] = comp["RT"] - rPredict(mod, FloatVector(comp["RT"]))
    elif type == "mz":  # m/z-calibration
        if (refMz == compMz).all():
            mzShifts = 1e-6 * np.random.normal(len(refMz))
        else:
            mzShifts = (compMz - refMz) / refMz * 1e6
        mod = rLoess(FloatVector(compMz), FloatVector(mzShifts), 1, "aicc", "gaussian")
        compMz = compMz * (1 + np.array(mod.rx2("fitted")) / 1e6)  # Calibrated m/z basd on the model

        # Calculate a new (dynamic) m/z-tolerance
        if (refMz == compMz).all():
            mzShifts = 1e-6 * np.random.normal(len(refMz))
        else:
            mzShifts = (compMz - refMz) / refMz * 1e6  # Calibrated m/z-shifts
        modMzSd = rLoess(FloatVector(compMz), FloatVector(mzShifts ** 2), 1, "aicc", "gaussian")
        mzSd = np.sqrt(np.maximum(0, rPredict(modMzSd, FloatVector(ref["mz"]))))

        # Calculate a new (dynamic) RT-tolerance
        if (refRt == compRt).all():
            rtShifts = 1e-6 * np.random.normal(len(refRt))
        else:
            rtShifts = compRt - refRt
        ind = np.where((rtShifts >= np.percentile(rtShifts, 10)) &
                       (rtShifts <= np.percentile(rtShifts, 90)))[0]
        modRtSd = rLoess(FloatVector(compRt[ind]), FloatVector(rtShifts[ind] ** 2), 1, "aicc", "gaussian")
        rtSd = np.sqrt(np.maximum(0, rPredict(modRtSd, FloatVector(ref["RT"]))))

        # Calibration of the entire comp["mz"]
        comp["mz"] = comp["mz"] / (1 + np.array(rPredict(mod, FloatVector(comp["mz"]))) / 1e6)

    return ref, comp, rtSd, mzSd


def findComparableFeatures(ref, comp, rtSd, mzSd, params):
    ref = np.sort(ref, order="intensity")[::-1]  # Sort features in descending order of intensity
    comp = np.sort(comp, order="intensity")[::-1]

    n = ref.shape[0]
    if not isinstance(rtSd, (list, np.ndarray)):
        rtSd = np.repeat(rtSd, n)
    if not isinstance(mzSd, (list, np.ndarray)):
        mzSd = np.repeat(mzSd, n)
    rtSd[rtSd == 0] = min(rtSd[rtSd > 0])
    mzSd[mzSd == 0] = min(mzSd[mzSd > 0])
    sdWidth = float(params["sd_width"])
    rtTol = rtSd * sdWidth
    mzTol = mzSd * sdWidth

    # Look for matching features between "ref" and "comp"
    refInd, compInd = [], []
    for i in range(n):  # For each feature in "ref", look for a matching one in "comp"
        z = ref["z"][i]
        mz = ref["mz"][i]
        rt = ref["RT"][i]
        rtDev = comp["RT"] - rt
        mzDev = (comp["mz"] - mz) / mz * 1e6  # Unit of PPM
        if z == 0:  # Undetermined charge
            # For the feature with undetermined charge,
            # look for a matching one without considering charge state
            rowInd = np.where((abs(rtDev) <= rtTol[i]) &
                              (abs(mzDev) <= mzTol[i]))[0]
        else:
            # For the feature with a charge state,
            # look for a matching one with considering charge state
            rowInd = np.where((abs(rtDev) <= rtTol[i]) &
                              (abs(mzDev) <= mzTol[i]) &
                              (comp['z'] == z))[0]
        if len(rowInd) > 0:
            # When multiple features in "comp" are matched to a feature in "ref",
            # choose the one with the highest intensity
            # Since "comp" is sorted by descending order of intensity,
            # the first one has the highest intensity
            rowInd = rowInd[0]
            if rowInd in compInd:
                continue
            else:
                refInd.append(i)
                compInd.append(rowInd)

    # # Optional/advanced function for rescuing some "unaligned" features by loosening tolerances
    # if step != "calibration" and params["rescue"] == "1":
    #     rtTolUnit = params["rt_tolerance_unit"].split(",")
    #     rtTol = params["rt_tolerance_value"].split(",")
    #     mzTolUnit = params["mz_tolerance_unit"].split(",")
    #     mzTol = params["mz_tolerance_value"].split(",")
    #     for i in range(len(rtTol)):
    #         refInd, compInd = rescueComparableFeatures(ref, comp, refInd, compInd, rtSd, mzSd, rtTolUnit[i], rtTol[i], mzTolUnit[i], mzTol[i])

    return refInd, compInd


def rescueComparableFeatures(ref, comp, refInd, compInd, rtSd, mzSd, rtTolUnit, rtTol, mzTolUnit, mzTol):
    nRescue = 0
    nUnaligned = comp.shape[0] - len(compInd)
    print("    There are %d unaligned features" % nUnaligned)
    logging.info("    There are %d unaligned features" % nUnaligned)

    # 1. Reduce variables to ones of unaligned features
    n = ref.shape[0]
    uRefInd = np.setdiff1d(range(n), refInd)
    rtSd = rtSd[uRefInd]
    mzSd = mzSd[uRefInd]
    uRef = ref[uRefInd]
    nc = comp.shape[0]
    uCompInd = np.setdiff1d(range(nc), compInd)
    uComp = comp[uCompInd]

    # 2. Criteria for rescuing unaligned features
    #    - Absolute intensity level: hard-coded (grand median intensity of aligned features)
    #    - Intensity-ratio between ref and comp runs: hard-coded (within 95% of the ratios of aligned features)
    #    - RT- and m/z-shifts should be within specified tolerance (e.g. 10SD or 10ppm)
    intLevel = np.ma.median(comp["intensity"][compInd])
    ratioPct = 95
    ratioAligned = np.log2(ref["intensity"][refInd]) - np.log2(comp["intensity"][compInd])
    lRatio = np.percentile(ratioAligned, (100 - ratioPct) / 2)
    uRatio = np.percentile(ratioAligned, 100 - (100 - ratioPct) / 2)
    print("    - Intensity higher than %d (median intensity of aligned features)" % intLevel)
    print("    - Ratio between reference and compared runs within %d %%" % ratioPct)
    logging.info("    - Intensity higher than %d (median intensity of aligned features)" % intLevel)
    logging.info("    - Ratio between reference and compared runs within %d %%" % ratioPct)
    if rtTolUnit == "1":
        print("    - RT-shifts within %s x SD of estimated RT-shifts from aligned features" % rtTol)
        logging.info("    - RT-shifts within %s x SD of estimated RT-shifts from aligned features" % rtTol)
        rtTol = float(rtTol) * rtSd
    elif rtTolUnit == "2":
        print("    - RT-shifts less than %s seconds" % rtTol)
        logging.info("    - RT-shifts less than %s seconds" % rtTol)
        rtTol = np.repeat(float(rtTol), len(uRefInd))
    else:
        print("    WARNING: check your parameters for RT-tolerance unit. It should be either 1 or 2")
        print("    Due to incorrect parameter settings, the rescue step is skipped")
        logging.info("    WARNING: check your parameters for RT-tolerance unit. It should be either 1 or 2")
        logging.info("    Due to incorrect parameter settings, the rescue step is skipped")
        return refInd, compInd

    if mzTolUnit == "1":
        print("    - m/z-shifts within %s x SD of estimated m/z-shifts from aligned features" % mzTol)
        logging.info("    - m/z-shifts within %s x SD of estimated m/z-shifts from aligned features" % mzTol)
        mzTol = float(mzTol) * mzSd
    elif mzTolUnit == "2":
        print("    - m/z-shifts less than %s seconds" % mzTol)
        logging.info("    - m/z-shifts less than %s seconds" % mzTol)
        mzTol = np.repeat(float(mzTol), len(uRefInd))
    else:
        print("    WARNING: check your parameters for RT-tolerance unit. It should be either 1 or 2")
        print("    Due to incorrect parameter settings, the rescue step is skipped")
        logging.info("    WARNING: check your parameters for RT-tolerance unit. It should be either 1 or 2")
        logging.info("    Due to incorrect parameter settings, the rescue step is skipped")
        return refInd, compInd

    # 3. Apply the criteria to unaligned features
    for i in range(len(uRefInd)):
        intRatios = np.log2(uComp["intensity"]) - np.log2(uRef["intensity"][i])
        rtShifts = uComp["RT"] - uRef["RT"][i]
        mzShifts = (uComp["mz"] - uRef["mz"][i]) / uRef["mz"][i] * 1e6
        ind = np.where((uComp["intensity"] > intLevel) & (intRatios >= lRatio) & (intRatios <= uRatio) &
                       (abs(mzShifts) <= mzTol[i]) & (abs(rtShifts) <= rtTol[i]))[0]
        if len(ind) > 0:
            charge = uRef["z"][i]
            if charge == 0:
                ind = ind[0]
                if uCompInd[ind] not in compInd:
                    refInd.append(uRefInd[i])
                    compInd.append(uCompInd[ind])
                    nRescue += 1
            else:
                ind2 = np.where(uComp["z"][ind] == charge)[0]
                if len(
                        ind2) > 0:  # When there is/are feature(s) matched to the reference one in terms of mz, rt and charge
                    ind = ind[ind2[0]]
                    if uCompInd[ind] not in compInd:
                        refInd.append(uRefInd[i])
                        compInd.append(uCompInd[ind])
                        nRescue += 1
                else:
                    ind2 = np.where(uComp["z"][ind] == 0)[0]
                    if len(ind2) > 0:
                        ind = ind[ind2[0]]
                        if uCompInd[ind] not in compInd:
                            refInd.append(uRefInd[i])
                            compInd.append(uCompInd[ind])
                            nRescue += 1
    print("    Through the rescue procedure %d features are additionally aligned" % nRescue)
    logging.info("    Through the rescue procedure %d features are additionally aligned" % nRescue)
    return (refInd, compInd)


def matchFeatures(ref, comp, rtSd, mzSd, params):
    ref = np.sort(ref, order="intensity")[::-1]  # Sort features in descending order of intensity
    comp = np.sort(comp, order="intensity")[::-1]

    n = ref.shape[0]
    if not isinstance(rtSd, (list, np.ndarray)):
        rtSd = np.repeat(rtSd, n)
    if not isinstance(mzSd, (list, np.ndarray)):
        mzSd = np.repeat(mzSd, n)
    rtSd[rtSd == 0] = min(rtSd[rtSd > 0])
    mzSd[mzSd == 0] = min(mzSd[mzSd > 0])
    sdWidth = float(params["sd_width"])
    rtTol = rtSd * sdWidth
    mzTol = mzSd * sdWidth

    # Look for matching features between "ref" and "comp"
    refInd, compInd = [], []
    for i in range(n):  # For each feature in "ref", look for a matching one in "comp"
        z = ref["z"][i]
        mz = ref["mz"][i]
        rt = ref["RT"][i]
        rtDev = comp["RT"] - rt
        mzDev = (comp["mz"] - mz) / mz * 1e6  # Unit of PPM
        rowInd = np.where((abs(rtDev) <= rtTol[i]) & (abs(mzDev) <= mzTol[i]))[0]

        if len(rowInd) > 0:
            if z == 0:  # When the reference feature's charge is undetermined
                rowInd = rowInd[0]  # Choose the strongest feature as matched
                if rowInd not in compInd:
                    refInd.append(i)
                    compInd.append(rowInd)
            else:
                rowInd2 = np.where(comp["z"][rowInd] == z)[0]
                if len(
                        rowInd2) > 0:  # When there is/are feature(s) matched to the reference one in terms of rt, mz and charge
                    rowInd = rowInd[rowInd2[0]]  # Choose the strongest among matched features
                    if rowInd not in compInd:
                        refInd.append(i)
                        compInd.append(rowInd)
                else:  # When there's no feature having the same charge as the reference feature
                    rowInd2 = np.where(comp["z"][rowInd] == 0)[
                        0]  # Find "comp" feature(s) with charge = 0 (and rtDev and mzDev are within tolerances)
                    if len(rowInd2) > 0:
                        rowInd = rowInd[rowInd2[0]]  # Choose the strongest among matched features with charge = 0
                        if rowInd not in compInd:
                            refInd.append(i)
                            compInd.append(rowInd)

    return refInd, compInd


def findMatchedFeatures(refNo, fArray, rtSdArray, mzSdArray, fNames, params):
    n = len(fNames)
    indArray = - np.ones((fArray[refNo].shape[0], n), dtype=int)
    indArray[:, refNo] = range(fArray[refNo].shape[0])
    for i in range(n):
        if i != refNo:
            refName = os.path.basename(fNames[refNo])
            compName = os.path.basename(fNames[i])
            print("  " + refName + ": %d features (reference run)" % fArray[refNo].shape[0])
            print("  " + compName + ": %d features (compared run)" % fArray[i].shape[0])
            logging.info("  " + refName + ": %d features (reference run)" % fArray[refNo].shape[0])
            logging.info("  " + compName + ": %d features (compared run)" % fArray[i].shape[0])
            refInd, compInd = matchFeatures(fArray[refNo], fArray[i], rtSdArray[i], mzSdArray[i], params)

            # Alignment indication
            # indArray is [# reference features x # feature files] matrix
            # If reference is the 3rd run (i.e. feature file),
            # indArray[i, j] indicates the feature index of j-th run matched to i-th reference feature
            # indArray[i, j] = -1 means that there's no feature matched to i-th reference feature in j-th run
            # rowInd = np.nonzero(np.in1d(indArray[:, refNo], refInd))[0]
            # indArray[rowInd, i] = compInd
            # print("  %d features are aligned between runs" % len(rowInd))
            indArray[refInd, i] = compInd
            print("  %d features are aligned between runs" % len(refInd))
            logging.info("  %d features are aligned between runs" % len(refInd))

            # Optional/advanced function for rescuing some "unaligned" features by loosening tolerances
            if params["rescue"] == "1":
                rtTolUnit = params["rt_tolerance_unit"].split(",")
                rtTol = params["rt_tolerance_value"].split(",")
                mzTolUnit = params["mz_tolerance_unit"].split(",")
                mzTol = params["mz_tolerance_value"].split(",")
                for j in range(len(rtTol)):
                    refInd, compInd = rescueComparableFeatures(fArray[refNo], fArray[i], refInd, compInd, rtSdArray[i],
                                                               mzSdArray[i],
                                                               rtTolUnit[j], rtTol[j], mzTolUnit[j], mzTol[j])
                    # rowInd = np.nonzero(np.in1d(indArray[:, refNo], refInd))[0]
                    # indArray[rowInd, i] = compInd
                    indArray[refInd, i] = compInd
    print()
    logging.info("")

    # Indexes of fully- and partially-aligned features
    fullInd = range(indArray.shape[0])
    partialInd = []
    colNames = list(fArray[0].dtype.names)
    for i in range(n):
        fArray[i].dtype.names = [os.path.splitext(fNames[i])[0] + "_" + c for c in colNames]
        if i != refNo:
            fullInd = np.intersect1d(fullInd, np.where(indArray[:, i] >= 0)[0])
            partialInd = np.union1d(partialInd, np.where(indArray[:, i] >= 0)[0])
    partialInd = np.setdiff1d(partialInd, fullInd)
    partialInd = partialInd.astype(int)

    # Fully-aligned features
    full = None
    for i in range(n):
        full_i = fArray[i][indArray[fullInd, i]]
        # full_i.dtype.names = [fNames[i] + "_" + c for c in full_i.dtype.names]
        if i == 0:
            full = full_i
        else:
            full = merge_arrays((full, full_i), asrecarray=True, flatten=True)

    # Partially-aligned features
    partial = None
    if partialInd.size > 0:
        for i in range(len(partialInd)):
            for j in range(n):
                colInd = indArray[partialInd[i], j]
                if j == 0:
                    if colInd == -1:
                        nullFeature = np.zeros((1,), dtype=fArray[j].dtype)
                        nullFeature[:] = np.nan
                        line = nullFeature
                    else:
                        line = fArray[j][colInd]
                else:
                    if colInd == -1:
                        nullFeature = np.zeros((1,), dtype=fArray[j].dtype)
                        nullFeature[:] = np.nan
                        line = merge_arrays((line, nullFeature), asrecarray=True, flatten=True)
                    else:
                        line = merge_arrays((line, fArray[j][colInd]), asrecarray=True, flatten=True)
            if i == 0:
                partial = line
            else:
                partial = stack_arrays((partial, line), asrecarray=True, usemask=False)

    # Un-aligned features
    unaligned = []
    for i in range(n):
        if i == refNo:
            rowInd = np.setdiff1d(range(fArray[refNo].shape[0]), np.union1d(fullInd, partialInd))
            unaligned.append(fArray[refNo][rowInd,])
        else:
            rowInd = np.setdiff1d(range(fArray[i].shape[0]), indArray[:, i])
            unaligned.append(fArray[i][rowInd,])

    # Alignment summary
    print("  Alignment/match summary")
    print("  =======================")
    print("  After alignment/feature matching, fully-, partially- and un-aligned features are as follows")
    print("  Filename\t\tfully-aligned\tpartially-aligned\tun-aligned")
    logging.info("  Alignment/match summary")
    logging.info("  =======================")
    logging.info("  After alignment/feature matching, fully-, partially- and un-aligned features are as follows")
    logging.info("  Filename\t\tfully-aligned\tpartially-aligned\tun-aligned")
    nFull = len(fullInd)
    for i in range(n):
        if i == refNo:
            nPartial = len(partialInd)
        else:
            nPartial = np.sum(indArray[:, i] >= 0) - len(fullInd)
        nUn = unaligned[i].shape[0]
        print("  %s\t\t%d\t%d\t%d" % (fNames[i], nFull, nPartial, nUn))
        logging.info("  %s\t\t%d\t%d\t%d" % (fNames[i], nFull, nPartial, nUn))

    # Depending on the parameter, "pct_full_alignment", some partially-aligned features can be included in fully-aligned ones
    pctFullAlignment = float(params["pct_full_alignment"])
    if pctFullAlignment < 100:
        if n > 2:
            colNames = [col for col in partial.dtype.names if col.endswith('mz')]
            # nRuns = np.sum(np.array(partial[colNames].tolist()) > 0, axis=1)  # For each partially-aligned feature, the number of aligned runs (i.e. feature files)
            nRuns = np.sum(~np.isnan(np.array(partial[colNames].tolist())), axis = 1)
            rowInd = np.where(nRuns >= np.ceil(pctFullAlignment / 100 * n))[0]

            # Add some partially-aligned features to fully-aligned features
            full = stack_arrays((full, partial[rowInd]), asrecarray=True, usemask=False)
            partial = np.delete(partial, rowInd, axis=0)
            print("  According to the parameter setting, %d partially-aligned features are regarded as fully-aligned" % len(rowInd))
            logging.info("  According to the parameter setting, %d partially-aligned features are regarded as fully-aligned" % len(rowInd))
        else:
            print("  Since you have only two runs, there's no partially-aligned features")
            logging.info("  Since you have only two runs, there's no partially-aligned features")
    else:
        print("  According to the parameter setting, no feature is added to the set of fully-aligned ones")
        logging.info("  According to the parameter setting, no feature is added to the set of fully-aligned ones")
    print()
    logging.info("")
    return full, partial, unaligned


def alignFeatures(fArray, xmlFiles, paramFile):
    nFiles = len(xmlFiles)

    # Pandas dataframe to numpy structured array for internal computation
    for i in range(nFiles):
        fArray[i] = fArray[i].to_records(index=False)

    ###################
    # Load parameters #
    ###################
    params = utils.getParams(paramFile)
    # Features derived from feature files are stored in fArray. For example,
    # xmlFiles = [file1, file2, file3]
    # fArray[0] = features from file1 (which has column names like 'index', 'mz', etc.)
    # fArray[1] = features from file2
    # ...
    # The array of m/z values from the first feature file can be accessed by fArray[0]['mz']

    if nFiles > 1:  # Multiple feature files -> alignment is required
        print("  Feature calibration")
        print("  ===================")
        logging.info("  Feature calibration")
        logging.info("  ===================")

        ###################################
        # Selection of a reference sample #
        ###################################
        if params["reference_feature"] == "0":
            # A run with the largest median of top 100 intensities is set to a reference run
            refNo = 0
            refIntensity = 0
            for i in range(nFiles):
                tmpIntensity = np.median(sorted(fArray[i]["intensity"], reverse=True)[0: 100])
                if tmpIntensity >= refIntensity:
                    refNo = i
                    refIntensity = tmpIntensity
        elif params["reference_feature"] == "1":
            # A run with the most number of features is set to a reference run
            refNo = 0
            refN = 0
            for i in range(nFiles):
                tmpN = len(fArray[i])
                if tmpN >= refN:
                    refNo = i
                    refN = tmpN
        else:
            try:
                refNo = xmlFiles.index(params["reference_feature"])
            except:
                sys.exit("  'reference_feature' parameter should be correctly specified")
        print("  %s is chosen as the reference run" % os.path.basename(xmlFiles[refNo]))
        logging.info("  %s is chosen as the reference run" % os.path.basename(xmlFiles[refNo]))

        ############################################################
        # Calibration of features against those in a reference run #
        ############################################################
        rtSdArray, mzSdArray = [], []
        featureNames = []
        for i in range(nFiles):
            featureName = os.path.basename(xmlFiles[i])
            featureNames.append(featureName)
            if i != refNo:
                print("  " + featureName + " is being aligned against the reference run (it may take a while)")
                logging.info("  " + featureName + " is being aligned against the reference run (it may take a while)")
                fArray[i], rtSd, mzSd = calibrateFeatures(fArray[refNo], fArray[i], params)
                rtSdArray.append(rtSd)
                mzSdArray.append(mzSd)
            else:
                rtSdArray.append("NA")
                mzSdArray.append("NA")

        print("  Calibration summary")
        print("  ===================")
        print("  After calibration, RT- and m/z-shifts of each run (against the reference run) are centered to zero")
        print("  Variations (i.e. standard deviation) of RT- and m/z-shifts are as follows,")
        print("  Filename\t\t\t#features\tSD of RT-shifts [second]\tSD of m/z-shifts [ppm]")
        logging.info("  Calibration summary")
        logging.info("  ===================")
        logging.info("  After calibration, RT- and m/z-shifts of each run (against the reference run) are centered to zero")
        logging.info("  Variations (i.e. standard deviation) of RT- and m/z-shifts are as follows,")
        logging.info("  Filename\t\t\t#features\tSD of RT-shifts [second]\tSD of m/z-shifts [ppm]")
        for i in range(nFiles):
            nFeatures = str(fArray[i].shape[0])
            if i != refNo:
                meanRtSd = "%.6f" % np.mean(rtSdArray[i])
                meanMzSd = "%.6f" % np.mean(mzSdArray[i])
            else:
                meanRtSd = "NA"
                meanMzSd = "NA"
            print("  " + featureNames[i] + "\t\t\t" + nFeatures + "\t" + meanRtSd + "\t" + meanMzSd)
            logging.info("  " + featureNames[i] + "\t\t\t" + nFeatures + "\t" + meanRtSd + "\t" + meanMzSd)
        print()
        logging.info("")

        #################################################################
        # Identification of fully-aligned features for further analysis #
        #################################################################
        print("  Feature alignment")
        print("  =================")
        logging.info("  Feature alignment")
        logging.info("  =================")
        fullFeatures, partialFeatures, unalignedFeatures = findMatchedFeatures(refNo, fArray, rtSdArray, mzSdArray,
                                                                               featureNames, params)
    else:
        print("  Since a single feature is used, the feature alignment is skipped")
        logging.info("  Since a single feature is used, the feature alignment is skipped")
        fullFeatures = np.copy(fArray[0])  # Masked array to 2D numpy array
        colNames = list(fullFeatures.dtype.names)
        featureName = os.path.splitext(os.path.basename(xmlFiles[0]))[0]
        fullFeatures.dtype.names = [featureName + "_" + c for c in colNames]
        partialFeatures, unalignedFeatures = None, None

    ################################################################
    # Write fully-, partially- and/or un-aligned features to files #
    ################################################################
    # At this step, fully-, partially- and unaligned features are written to files and saved
    # Also, those features are converted to pandas DataFrame format and returned
    dfFull, dfPartial, dfArrayUnaligned = utils.generateFeatureFile(fullFeatures, partialFeatures, unalignedFeatures,
                                                                    params)

    return dfFull, dfPartial, dfArrayUnaligned
