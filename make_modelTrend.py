#!/usr/bin/env python

# to run the script with the correct version of uvcdat:
#  source /usr/local/uvcdat/1.4.0/bin/setup_runtime.sh

import cdms2
from cdms2 import MV2
import numpy
import glob
import sys
import os
from os import path
import re
import string
import random
import gc
import logging
import logging.handlers

# ____________________________
def usage():
    textUsage='make_modelTrend.py.\n\tComputes models linear trend (y = a . time + y0). Developped in first place for estimating variables projections divergence from control runs (ideally, a should be null).\n'
    return textUsage
# ____________________________
def exitMessage(msg, exitCode='1'):
    thisLogger.critical(msg)
    print msg
    print
    print usage()
    sys.exit(exitCode)
# __________________________
def getTrendType(argString):    
    if argString == "esm" or argString == "esmcontrol":
        return 'esmControl'
    elif argString == 'pi' or argString == 'pictrl' or argString == 'picontrol':
        return 'piControl'
    else return ''
# __________________________
def getListFromFile(infile):
    modeList=[]
    try:
        with open(infile, "r") as f:
            for textline in f:
                thisStr = textLine.replace(" ", "").replace('\n','')
                if not (thisStr==""):
                    modelList.append( thisStr )
    except IOError as e:
        exitMessage('I/O Error {1} while processing text file {0}:{2}. Exit(10).'.format(modelListFile, e.errno, e.strerror), 10)
    except:
        exitMessage('Unexpected error while processing text file {0}. Exit(11).'.format(modeListFile), 11)
# __________________________
# returns basefilenames, in alphabetical order
def selectModelFiles(indir,variable, frequency, iModel, trendType, rip):
    searchString = ''.join([ '{0}_'.format(s) for s in [a,b,c,d] if s != '']) + '*.nc'

    # get file list, non empty, sorted alphabetically
    listFile = [ os.path.basename(f) for f in glob.glob( os.path.join(indir, searchString) ) if (os.stat(f).st_size and pattern.match(os.path.basename(f) ) ) ]

    if len(listFile)=0: return None
    return sorted(listFile)
#___________________________
if __name__=="__main__":

    variable = None
    frequency = 'Omon'
    trendType = 'esmControl'
    rip='r1i1p1'
    indir = None
    tmpdir = None
    outdir = None
    modelListFile=None

    ii = 1
    while ii < len(sys.argv):
        arg = sys.argv[ii].lower()
        
        if arg == '-path':
            ii = ii + 1
            indir = sys.argv[ii]
        elif arg == '-outdir':
            ii = ii + 1
            outdir = sys.argv[ii]
        elif arg == '-tmpdir':
            ii = ii + 1
            tmpdir = sys.argv[ii]
        elif arg == '-v':
            ii = ii + 1
            variable = sys.argv[ii]
        elif arg == '-trendType':
            ii = ii + 1
            trendType = getTrendType(sys.argv[ii].lower())
        elif arg== '-rip':
            ii = ii + 1
            rip = sys.argv[ii].lower()
        elif arg=='-log':
            ii = ii + 1
            logFile = sys.argv[ii]
        ii = ii + 1

    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    thisLogger = logging.getLogger('MyLogger')
    thisLogger.setLevel(logging.DEBUG)
    handler = logging.handlers.RotatingFileHandler(logFile, maxBytes=1024*500, backupCount=5)
    thisLogger.addHandler(handler)

    # for netcdf3: set flag to 0
    cdms2.setNetcdfShuffleFlag(1)
    cdms2.setNetcdfDeflateFlag(1)
    cdms2.setNetcdfDeflateLevelFlag(3)

    # for each model
    listModels = getListFromFile(modelListFile)
    for iModel in listModels:
        # get the list of input files
        lstFiles = selectModelFiles(indir, variable, frequency, iModel, trendType, rip)
        if lstFiles is None: continue
        # sort them in chronological order
        # call the trend estimator
        (a, b)=do_trend(files)
        # save a and b for this model

# end of file
