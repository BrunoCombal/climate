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
    textUsage=''
    return textUsage

# ____________________________
def exitMessage(msg, exitCode='1'):
    thisLogger.critical(msg)
    print msg
    print
    print usage()
    sys.exit(exitCode)

#___________________________
if __name__=="__main__":

    variable = None
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



# end of file
