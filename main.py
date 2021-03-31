#!/usr/bin/python

import sys, os, re, logging, pandas as pd
import utils
from featureDetection import detectFeatures
from datetime import datetime


##################
# Initialization #
##################
# For desktop debugging,
paramFile = r"jumpm.params"
inputFiles = [r"/home/jcho/dev/spectralLibrary/FTLD_Batch2_F50.mzXML", 
              r"/home/jcho/dev/spectralLibrary/FTLD_Batch2_F51.mzXML",
              r"/home/jcho/dev/spectralLibrary/FTLD_Batch2_F52.mzXML"]
params = utils.getParams(paramFile)

skipScans = [1, 3, 5, 7, 10]
for skipScan in skipScans:
    params["skipping_scans"] = skipScan
    logFile = "jump_m.log"
    if os.path.exists(logFile):
        os.system("rm " + logFile)
    logging.basicConfig(format='%(message)s', filename=logFile, level=logging.INFO)

    print()
    print("  Jump -m started")
    logging.info("  Jump -m started")
    now = datetime.now()
    nowString = now.strftime("%Y/%m/%d %H:%M:%S")
    print("  " + nowString)
    logging.info("  " + nowString)

    #####################
    # Feature detection #
    #####################
    print()
    print("  #####################")
    print("  # Feature detection #")
    print("  #####################")
    logging.info("")
    logging.info("  #####################")
    logging.info("  # Feature detection #")
    logging.info("  #####################")
    featureArray = []
    nFiles = len(inputFiles)
    for i in range(nFiles):
        f = detectFeatures(inputFiles[i], params)

    print()
    logging.info("")

    print("  Jump -m finished")
    logging.info("  Jump -m finished")
    now = datetime.now()
    nowString = now.strftime("%Y/%m/%d %H:%M:%S")
    print("  " + nowString)
    logging.info("  " + nowString)

