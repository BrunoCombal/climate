#!/usr/bin/env python

# \author Bruno Combal, IOC-UNESCO
# \date December 2013
# many thanks to Mark Carson (ZMAW) for advising me on the computation and usage of control run with CMIP5 datasets.

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
    textUsage='make_modelTrend.py.\n\tComputes models linear trend (y = a . time + y0).\nDevelopped in first place for estimating variables projections divergence from control runs (ideally, a should be null).\n'
    textUsage='The output is an array of polynomial coefficients (DEGREEPOLYFIT+1 coefficients) for each grid point. The time step used for interpolating is given in year(float): year+(month-1)/12.0, starts counting at 1'
    textUsage='SYNOPSIS:\n\tmake_modelTrend.py -path|-p INPATH -outdir|-o OUTPATH [-tmpdir WRKPATH] -v VARIABLE -trendType TRENDTYPE [-rip RIP] [-log LOGFILE] [-degree|deg DEGREEPOLYFIT] [-annualAvg ANNUALAVG]'
    textUsage=textUsage+'\tDEGREEPOLYFIT: degree of the polynom to fit; 1 for linear, 2 for quadratic, etc.'
    textUsage=textUsage+'\tANNUALAVG: False/True computes the annual average from the time series before calling polyfit'
    return textUsage
# ____________________________
def exitMessage(msg, exitCode='1'):
    thisLogger.critical(msg)
    print msg
    print
    print usage()
    sys.exit(exitCode)
# ____________________________
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for x in range(size))
# __________________________
def getTrendType(argString):    
    if argString == "esm" or argString == "esmcontrol":
        return 'esmControl'
    elif argString == 'pi' or argString == 'pictrl' or argString == 'picontrol':
        return 'piControl'
    else: return ''
# __________________________
def getAnnualAvg(aa):
    if aa == 'none' or aa == 'false' or aa == 'no' or aa=='0':
        return False
    else: 
        return True
# __________________________
def getListFromFile(infile):
    modelList=[]
    try:
        with open(infile, "r") as f:
            for textLine in f:
                thisStr = textLine.replace(" ", "").replace('\n','')
                if not (thisStr==""):
                    modelList.append( thisStr )
    except IOError as e:
        exitMessage('I/O Error {1} while processing text file {0}:{2}. Exit(10).'.format(infile, e.errno, e.strerror), 10)
    except:
        exitMessage('Unexpected error while processing text file {0}. Exit(11).'.format(infile), 11)

    return modelList
# __________________________
# returns basefilenames, in alphabetical order
def selectModelFiles(indir,variable, frequency, iModel, trendType, rip):

    searchString = ''.join([ '{0}_'.format(s) for s in [variable, frequency, iModel, trendType, rip] if s != '']) + '*.nc'

    # get file list, non empty, sorted alphabetically
    listFile = [ os.path.basename(f) for f in glob.glob( os.path.join(indir, searchString) ) if (os.stat(f).st_size) ]

    if len(listFile)==0: return None
    return sorted(listFile)
# ___________________________
# average a series by steps of 12
# only keep complete years
def yearlyAvg(timeAxis, dataIn):
    dataOut=[]
    timeOut=[]
    for iyear in range(0, len(dataIn), 12):
        dataOut.append(numpy.average( dataIn[iyear:iyear+12] ))
        timeOut.append(numpy.average( timeAxis[iyear+6] ))

    return (timeOut, dataOut)
# ___________________________
# for this version, assume the list is sorted in chronological order
# the total amount of data won't fit in memory. Process by line
# Pyhton shows memory management issues: force deleting arrays and call garbage collector in loops
def do_trend(indir, fileList, variable, outfile, degree, annualAVG):

    newTimeAxis=None

    # open all files
    lstFID = []
    print 'fileList ',fileList
    for ii in fileList:
        try:
            fid = cdms2.open(os.path.join(indir, ii), 'r')
        except:
            thisLogger.warn('Could not open file {0}. Continue.'.format(ii))
            continue
        lstFID.append(fid)
    # if no file was open, return None
    if len(lstFID)==0: return (None,None)

    dims=lstFID[0][variable].shape # assume t, z,y,x or t, y,x
    
     # create time axis
    # timeAxis: year counts (on floating point), starting at 1
    # timeAxisOrg: original time values
    timeAxis=[]
    timeAxisOrg=[]
    for ifid in lstFID:
        thisTime = [ t.year + (t.month-1.0)/12.0 for t in ifid['time'].asComponentTime() ]
        thisTimeOrg = [ t.year*10000 + t.month*100 + t.day for t in ifid['time'].asComponentTime() ]
        timeAxis = numpy.concatenate( (timeAxis, thisTime), axis=0)
        timeAxisOrg = numpy.concatenate( (timeAxisOrg, thisTimeOrg), axis=0)

    # Some datasets do not correctly encode masks: quick fix=if no change (along time dimension), then mask
    lstTime=lstFID[0][variable].getTime().asComponentTime()
    wtk = cdms2.MV2.array(lstFID[0][variable].subRegion(time=lstTime[0])).mask.squeeze()
    redoMask = False
    
    if wtk.all(): redoMask=True
    if wtk is None: redoMask=True
    if wtk.ndim==0: redoMask=True
    # if the mask has only 'True' values, let's compute the mask
    if redoMask:
        thisLogger.info('Mask not found in the dataset, computing mask from the time series. Continue.')
        # 3 consecutives values should be enough
        test = numpy.array(lstFID[0][variable].subRegion(time=(lstTime[0], lstTime[2],'cc') ))
        wtk = ( numpy.max(test, axis=0) - numpy.min(test, axis=0) ) < 0.00001
    lstIdx = numpy.ndindex(dims[1:])

   # create trend coefficient matrix
    coeff = numpy.zeros( dims[1:]+(degree + 1,) ) + 1.e20

    if annualAVG:
        thisLogger.info('Annual average computed before calling polyfit (deg={0}).'.format(degree))
    else:
        thisLogger.info('NO annual average, raw data used for calling polyfit (deg={0}).'.format(degree))

    thisLogger.info('Iterating over {0}/{1} points.'.format(lstIdx.total, len(dims[1:]) ))
    for idx in lstIdx:
        if wtk[idx] == True:
            continue # if this pixel is in the mask, jump to the next iteration
        # get data from all files for this position idx
        data=None
        for ifid in lstFID: # concatenate thisData into data
            thisData = numpy.array(ifid[variable][:, idx[0], idx[1]]).ravel()
            if data is None:
                data = thisData.copy()
            else:
                try:
                    data = numpy.concatenate( (data, thisData), axis=0)
                except Exception, e:
                    print 'thisData.shape: {0}'.format(thisData.shape)
                    print 'data.shape: {0}'.format(data.shape)
                    print 'netcdf variable {0} shape: {1}'.format(variable, ifid[variable].shape)
                    exitMessage('Unexpected error while processing file in the series {0}, with error {1}. Exit(101).'.format(fileList[0],e), 101)

        if annualAVG:
            thisLogger.info('Computing yearly averages, from data.shape: {0} and len(timeAxis): {1}'.format(data.shape, len(timeAxis)))
            (newTimeAxis, yearlyData) = yearlyAvg(timeAxis, data)
            thisLogger.debug('yearly average produced len(newTimeAxis): {0} and len(yearlyData): {1}.'.format( len(newTimeAxis), len(yearlyData)))
            coeff[idx[0], idx[1],:] = numpy.polyfit(newTimeAxis, yearlyData, degree)
        else:
            coeff[idx[0], idx[1],:] = numpy.polyfit(timeAxis, data, degree)
        #print idx, coeff[idx[0], idx[1],:]

        del thisData
        del data
        gc.collect()

    # save result
    thisLogger.info('Saving results')
    outfid=cdms2.open(outfile, 'w')

    outvar=cdms2.createVariable(coeff, id='coeff')
    outfid.write(outvar)
    if newTimeAxis is None:
        outtime=cdms2.createVariable(timeAxis, id='timeAxis')
    else:
        outtime=cdms2.createVariable(newTimeAxis, id='timeAxis')
    outfid.write(outtime)

    outtimeorg=cdms2.createVariable(timeAxisOrg, id='timeAxisOrg')
    outfid.write(outtimeorg)
    outfid.close()

    # close fid
    for ifid in lstFID: ifid.close()

    # some cleaning: python is not very good with collections
    del lstIdx
    del wtk
    gc.collect()
# ___________________________
if __name__=="__main__":

    variable = None
    frequency = 'Omon'
    trendType = 'esmControl'
    rip='r1i1p1'
    indir = None
    tmpdir = None
    outdir = None
    modelListFile=None
    logFile='{0}.log'.format(__file__)
    degree=1
    annualAVG=False

    ii = 1
    while ii < len(sys.argv):
        arg = sys.argv[ii].lower()
        
        if arg == '-path' or arg=='-p':
            ii = ii + 1
            indir = sys.argv[ii]
        elif arg == '-outdir' or arg=='-o':
            ii = ii + 1
            outdir = sys.argv[ii]
        elif arg == '-tmpdir':
            ii = ii + 1
            tmpdir = sys.argv[ii]
        elif arg == '-v':
            ii = ii + 1
            variable = sys.argv[ii]
        elif arg=='-deg' or arg=='-degree':
            ii = ii + 1
            degree = int(sys.argv[ii])
        elif arg=='-annualavg':
            ii=ii+1
            annualAVG=getAnnualAvg(sys.argv[ii].lower())
        elif arg == '-modellist':
            ii = ii + 1
            modelListFile = sys.argv[ii]
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

    # logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    thisLogger = logging.getLogger('MyLogger')
    thisLogger.setLevel(logging.DEBUG)
    handler = logging.handlers.RotatingFileHandler(logFile, maxBytes=1024*500, backupCount=5)
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(formatter)
    thisLogger.addHandler(handler)

    thisLogger.info('============ Starting logger ===========')

    # process input parameter
    if variable is None:
        exitMessage('Missing a variable name to process. Exit(1).', 1)
    if indir is None:
        exitMessage('Missing an input data path. Exit(5)',5)
    if outdir is None:
        exitMessage('Missing an output directory. Exit(2).',2)
    if modelListFile is None:
        exitMessage('Missing a model list file, use option -modellist. Exit(3).',3)
    if not os.path.exists(modelListFile):
        exitMessage('Model list file {0} not found. Exit(4).',4)
    if tmpdir is None:
        tmpdir = '{0}/tmp_{1}'.format(outdir, id_generator() )

    if not os.path.exists(outdir): os.makedirs(outdir)
    if not os.path.exists(tmpdir): os.makedirs(tmpdir) 

    # for netcdf3: set flag to 0
    cdms2.setNetcdfShuffleFlag(1)
    cdms2.setNetcdfDeflateFlag(1)
    cdms2.setNetcdfDeflateLevelFlag(3)

    # for each model
    listModels = getListFromFile(modelListFile)
    for iModel in listModels:
        # get the list of input files
        thisLogger.info( 'processing model {0}'.format(iModel))
        lstFiles = selectModelFiles(indir, variable, frequency, iModel, trendType, rip)
        if lstFiles is None: continue
        # sort them in chronological order
        # call the trend estimator
        outfile='{0}/trend_{1}_{2}_{3}_{4}_{5}.nc'.format(outdir, variable, frequency, iModel, trendType, rip)
        do_trend(indir, lstFiles, variable, outfile, degree, annualAVG)

    thisLogger.info('========== end of processing ==========')

# end of file
