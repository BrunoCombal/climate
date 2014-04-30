#!/usr/bin/env python

# \author Bruno Combal, IOC-UNESCO
# \date November 2013

# to run the script with the correct version of uvcdat:
#  source /usr/local/uvcdat/1.2.0/bin/setup_cdat.sh

# computation of DHM
# 1./ compute a model climatology (ready made) and get an actual climato
# 1.a/ compute current model (SST - climato) -> delta
# 1.b/ apply delta to actual climato -> corrected SST
# 2./ From Corrected SST, compute DHM
# 2.a/ DHM: compare sum(past 4 months to now) to max(actual climato): 2 degrees above-> level 2

import cdms2
#from cdms2 import MV
import numpy
#import glob
import sys
import os
import string
import random
#from os import path
#import re
from scipy import interpolate
# import shutil
import logging
import logging.handlers

# ____________________________
def usage():
    textUsage = 'SYNOPSIS:\n\t{0} -o|-outdir  PATHOUT [-outpref OUTPREFIX] -input|-in|-i PATHIN PREFIXIN [-tmpdir TMPPATH]\n'.format(__file__)
    textUsage = textUsage + '\t[-var VARIABLE]\n'
    textUsage = textUsage + '\t-clim CLIMATO CLIMATOMAX CLIMRMSATMAX\n'
    textUsage = textUsage + '\t-modelClim MODELCLIMPATH MODELCLIMPREF\n'
    textUsage = textUsage + '\t-decad DECAD [-log LOGFILE]\n'
    textUsage = textUsage + '\n\tPATHOUT: output directory, created if does not exist;\n'
    textUsage = textUsage + '\tOUTPREF: prefix for the output name, default: dhm_;\n'
    textUsage = textUsage + '\tPATHIN: input data directory (does not support sub-directories);\n'
    textUsage = textUsage + '\tPREFIXIN: prefix of the input files;\n'
    textUsage = textUsage + '\tTMPPATH: temporary path. Default: a random pathname is defined at runtime, as a leaf of PATHOUT;\n'
    textUsage = textUsage + '\tVARIABLE: netcdf variable name to use for processing. Default is tos;\n'
    textUsage = textUsage + '\tCLIMATO: a climatology file, same grid as the SST, with 12 months;\n'
    textUsage = textUsage + '\tCLIMATOMAX: max value for the climatology;\n'
    textUsage = textUsage + '\tCLIMTSMATMAX: RMS observed at the climatology max;\n'
    textUsage = textUsage + '\tMODELCLIMPATH: path to models historical values (before projections);\n'
    textUsage = textUsage + '\tMODELCLIMPREF: root name (prefix) for the historical files;\n'
    textUsage = textUsage + '\tDECAD: the Year at which start counting a decad for the final synthesis;\n'
    textUsage = textUsage + '\tLOGFILE: a logfile name. Default: {0}.log'.format(__file__)

    print textUsage
# ____________________________
def exitMessage(msg, exitCode='1'):
    thisLogger.critical(msg)
    print msg
    print
    print usage()
    sys.exit(exitCode)
# _________________________________
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for x in range(size))
# _________________________________
# return the month count from year 0
def yyyymm2count(year, month):
    return year*12+month
# ________________
def count2yyyymm(count):
    year=int((count-1)/12)
    month = (count-12*year)
    return (year, month)
def makeGrid():
    xstart=0
    xend=360
    xstep=0.5
    ystart=-85
    yend=85
    ystep=0.5

    lon_bnds=[]
    lon=[]
    for ii in numpy.arange(xstart, xend, xstep):
        lon_bnds.append( [ii, ii+xstep] )
        lon.append(ii+0.5*xstep)
    lon_bnds=numpy.array(lon_bnds)
    lon=numpy.array(lon)

    lat_bnds=[]
    lat=[]
    for ii in numpy.arange(ystart, yend, ystep):
        lat_bnds.append([ii, ii+ystep])
        lat.append(ii+0.5*ystep)
    lat_bnds=numpy.array(lat_bnds)
    lat=numpy.array(lat)

    latAxis = cdms2.createAxis(lat, lat_bnds)
    latAxis.designateLatitude(True)
    latAxis.units='degrees_north'
    latAxis.long_name='Latitude'
    latAxis.id='latitude'

    lonAxis = cdms2.createAxis(lon, lon_bnds)
    lonAxis.designateLongitude(True, 360.0)
    lonAxis.units='degrees_east'
    lonAxis.id='longitude'
    lonAxis.long_name='Longitude'

    return((cdms2.createGenericGrid(latAxis, lonAxis, lat_bnds, lon_bnds), latAxis, lonAxis, lat_bnds, lon_bnds))
# _______________
def saveData(outfilename, data, typecode, id, fill_value, grid, copyaxes, attribute1, attribute2, latAxis, lonAxis):
    
    # for netcdf3: set flags to 0
    cdms2.setNetcdfShuffleFlag(0) #1
    cdms2.setNetcdfDeflateFlag(0) #1
    cdms2.setNetcdfDeflateLevelFlag(0) #3
    
    thisLogger.info('Saving data {0} to file {1}'.format(id, outfilename))
    if os.path.exists(outfilename): os.remove(outfilename)
    outfile = cdms2.open( outfilename, 'w')
    var = cdms2.createVariable(data, typecode=typecode, id=id, fill_value=fill_value, grid=grid, copyaxes=copyaxes, attributes=dict(long_name=attribute1, units=attribute2) )
    var.setAxisList((latAxis, lonAxis))
    outfile.write(var)
    outfile.close()
# _______________
# data have y first, then x
def do_resize(var, infile):
    yvar=infile[var][:] # y, x
    dim=yvar.shape

    points=numpy.zeros( (dim[0], dim[1], 2))
    for jj in range(dim[0]):
        for ii in range(dim[1]):
            points[jj, ii, 0] = jj
            points[jj, ii, 1] = ii

    outGY, outGX = numpy.mgrid[ -85 + 85 :85 + 85 :0.5 , 0:360:0.5 ]

    outvar = interpolate.griddata(numpy.reshape(points, (dim[0]*dim[1],2)), numpy.ravel(yvar), (outGY, outGX), method='linear')

    return outvar
# _______________
def readVar(var, infile):
    fh = cdms2.open(infile, 'r')
    data = fh[var][:]
    fh.close()
    return data
# _______________
# arrays were interpolated: creation of false 'nodata' values along the coast.
# no data are all reset
def resetNoData(data, limit, nodata):
    wtnd = data>limit
    if wtnd.any():
        data[wtnd] = 1.e20
    return data
# _______________
# note: real Climato has a coarser resolution
def do_dhm(var, inhist, modelClimatoRootName, indir, sstRootName, realClimato, maxRealClimato, realClimRMSAtMaxSST, outdir, dhmRootName, yearList):

    nodata = 1.e20
    (referenceGrid, latAxis, lonAxis, latBounds, lonBounds) = makeGrid()

    # open climatoMas
    climMax = readVar('sst', maxRealClimato)
    climMax = resetNoData(climMax, 50, nodata)
    RMSatMaxSST=readVar('rms_at_max',realClimRMSAtMaxSST)
    RMSatMaxSST=resetNoData(RMSatMaxSST, 50, nodata)
    climMax = climMax + RMSatMaxSST

    # open realClimato
    realClim = readVar('sst', realClimato)
    realClim = resetNoData(realClim, 50, 1.e20) # realClim a le meilleur masque de nodata

    # open the 12 model climatos
    modelClim=[]
    # initial value for frequency is set to nodata=1.e20
    frequencyLvl2=numpy.zeros( realClim[0].shape[0] * realClim[0].shape[1] ) # + nodata
    for imonth in range(1,13):
        thisName=os.path.join(inhist, '{0}{1:02}.nc'.format(modelClimatoRootName, imonth))
        if not os.path.exists(thisName): 
            print 'missing model climato {0}. Exit.'.format(thisName)
            sys.exit(1)
        modelClim.append(cdms2.open(thisName,'r'))

    # processing loop
    for iyear in yearList:

        dhmYearly = numpy.zeros( realClim[0].shape[0] * realClim[0].shape[1] ) # defined as the max dhm for this year
        for imonth in range(1,13):

            dhm = numpy.zeros( realClim[0].shape[0] * realClim[0].shape[1] )
            for ishift in 0, -1, -2, -3:
                shiftYear, shiftMonth = count2yyyymm( yyyymm2count(iyear, imonth) + ishift)
                #print 'rolling window: ', iyear, imonth, ishift, ' >> ', shiftYear, shiftMonth
                sstFName = os.path.join(indir, '{0}{1}{2:02}.nc'.format(sstRootName, shiftYear, shiftMonth))
                if not os.path.exists(sstFName):
                    thisLogger.critical('in do_dhm: file {0} not found. Exit(100).'.format(sstFName))
                    exitMessage('in do_dhm: file {0} not found. Exit(100).'.format(sstFName),100)

                thisModelSST = cdms2.open( sstFName, 'r' )

                # thisModelSST: K ; modelClim: K ; realClimato: C
                # thisModelSST: continents=1.e20, modelClim=continents=1.e20
                tosCorrected = realClim[ shiftMonth-1, :, : ] + ( thisModelSST[var][:] - modelClim[ shiftMonth-1 ]['Band1'][:])

                if ishift == 0:
                    delta = thisModelSST[var][:] - modelClim[ shiftMonth-1 ]['Band1'][:]

                thisModelSST.close()

                ##
                # wtk = tosCorrected.ravel() > ( numpy.ravel(climMax + RMSatMaxSST) )
                ## use an already corrected version of climMax
                wtk = numpy.ravel(tosCorrected) > numpy.ravel(climMax)
                if wtk.any():
                    ##
                    # dhm[wtk] = dhm[wtk] + numpy.ravel(tosCorrected - climMax - RMSatMaxSST)[wtk]
                    ##
                    dhm[wtk] = dhm[wtk] + numpy.ravel(tosCorrected - climMax )[wtk]

            # reset nodata: from realClim and modelClim
            wnodata = numpy.ravel(modelClim[ 0 ]['Band1'][:]) >= nodata
            if wnodata.any():
                dhm[wnodata] = nodata
                dhmYearly[wnodata] = nodata
            wnodata = realClim[0,:,:].ravel() >= nodata
            if wnodata.any():
                dhm[wnodata] = nodata
                dhmYearly[wnodata] = nodata

            wtk = dhm > dhmYearly
            if wtk.any():
                dhmYearly[wtk] = dhm[wtk]

            # write output
            outfilename=os.path.join( outdir , '{0}{1}{2:02}.nc'.format(dhmRootName, iyear, imonth) )
            saveData(outfilename, dhm.reshape( (realClim[0].shape[0] , realClim[0].shape[1])  ), 'f', 'dhm', 1.e20, referenceGrid, 1, 'dhm', 'None', latAxis, lonAxis)
        # all months done for this year    
        saveData( os.path.join( outdir , '{0}{1}.nc'.format(dhmRootName, iyear)), dhmYearly.reshape( (realClim[0].shape[0] , realClim[0].shape[1])  ), 'f', 'dhm', 1.e20, referenceGrid, 1, 'dhm', 'None', latAxis, lonAxis)

        # update frequency
        # let's update frequencies having real values
        wtk = dhmYearly > 2 # * (frequencyLvl2 < nodata) 
        if wtk.any():
            frequencyLvl2[wtk] = frequencyLvl2[wtk] + 1
        ## let's consider frequencies never updated so far (and identified as nodata)
        #wtk = (dhmYearly > 2 ) * (frequencyLvl2 >= nodata)
        #if wtk.any():
        #    frequencyLvl2[wtk] = 1


    wnodata = numpy.ravel(modelClim[ 0 ]['Band1'][:]) >= nodata
    if wnodata.any():
        frequencyLvl2[wnodata] = -1
        wnodata = realClim[0,:,:].ravel() >= nodata
        if wnodata.any():
            frequencyLvl2[wnodata] = -1

    # reset nodata to -1
    wrecode = frequencyLvl2 >= nodata
    if wrecode.any():
        frequencyLvl2[wrecode] = -1

    # save frequency
    saveData(os.path.join( outdir , '{0}_{1}.nc'.format('frequency_lvl2_', yearList[0])), frequencyLvl2.reshape( (realClim[0].shape[0] , realClim[0].shape[1])  ), 'i', 'lvl2_freq', -1, referenceGrid, 1, 'dhm', 'None', latAxis, lonAxis)
    for ii in modelClim:
        ii.close()
# _______________
if __name__=="__main__":

    indir=None #'/data/cmip5/rcp/rcp8.5/tos_ensemble/'#    indir='/data/cmip5/rcp/rcp8.5/tos4.5_ensemble/'
    outdir=None #'/data/cmip5/rcp/rcp8.5/tos_ensemble/'#    outdir='/data/cmip5/rcp/rcp8.5/tos4.5_ensemble/'
    tmpdir=None #'/data/tmp/'

    inhist=None #'/data/cmip5/rcp/rcp8.5/toshist_ensemble/'
    modelClimatoRootName=None #'climato_tos_1971_2000_' # climato 1980-2000 based on model output ensemble mean, make_tos_climato.sh

    realClimato=None #'/data/sst/reynolds_climatology/noaa_oist_v2/resized_fitted/sst.ltm.1971-2000_resized.nc' # has continent=1.e20
    maxRealClimato=None #'/data/sst/reynolds_climatology/noaa_oist_v2/resized_fitted/max_sst.ltm.1971-2000_resized.nc' # has continent=1.e20
    realClimRMSAtMaxSST=None #'/data/sst/reynolds_climatology/noaa_oist_v2/resized_fitted/rms_at_maxsst_resized.nc' # has continent outline

    sstRootName=None #'modelmean_tos_' # ensemble mean of projections
    dhmRootName='dhm_'
    dekad=None #2050
    var='tos'

    logFile='{0}.log'.format(__file__)


    ii = 1
    while ii < len(sys.argv):
        arg=sys.argv[ii].lower()

        if arg == '-outdir' or arg == '-o':
            ii = ii + 1
            outdir=sys.argv[ii]

        elif arg=='-outpref':
            ii = ii + 1
            dhmRootName=sys.argv[ii]

        elif arg == '-input' or arg=='-i' or arg=='-in':
            ii = ii + 1
            indir=sys.argv[ii]
            ii = ii + 1
            sstRootName = sys.argv[ii]

        elif arg == '-tmpdir':
            ii = ii + 1
            tmpdir=sys.argv[ii]

        elif arg == '-var':
            ii = ii + 1
            var=sys.argv[ii]

        elif arg=='-decad':
            ii = ii + 1
            dekad = int(sys.argv[ii])

        elif arg=='-modelclim':
            ii = ii + 1
            inhist = sys.argv[ii]
            ii = ii + 1
            modelClimatoRootName = sys.argv[ii]

        elif arg=='-clim':
            ii = ii + 1
            realClimato = sys.argv[ii]
            ii = ii + 1
            maxRealClimato = sys.argv[ii]
            ii = ii  + 1
            realClimRMSAtMaxSST = sys.argv[ii]

        elif arg=='-log':
            ii = ii + 1
            logFile = sys.argv[ii]
        ii = ii + 1

    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    thisLogger = logging.getLogger('MyLogger')
    thisLogger.setLevel(logging.DEBUG)
    handler = logging.handlers.RotatingFileHandler(logFile, maxBytes=1024*500, backupCount=5)
    thisLogger.addHandler(handler)

    # check input parameters
    if indir is None:
        exitMessage('Missing input directory, use option -input. Exit(2).',2)
    if sstRootName is None:
        exitMessage('Missing input files prefix, use option -input. Exit(2).',2)
    if outdir is None:
        exitMessage('Missing output directory, use option -outdir. Exit(3).', 3)
    if inhist is None:
        exitMessage('Missing directory for historical data, use option -modelClim. Exit(4).',4)
    if modelClimatoRootName is None:
        exitMessage('Missing prefix for historical dataset, use option -modelClim. Exit(4).',4)
    if realClimato is None:
        exitMessage('Missing file for real climatology, use option -clim. Exit(5).',5)
    if maxRealClimato is None:
        exitMessage('Missing file for real climatology max, use option -clim. Exit(5).',5)
    if realClimRMSAtMaxSST is None:
        exitMessage('Missing file for climato RMS at max SST, use option -clim. Exit(5).',5)
    if dhmRootName is None:
        exitMessage('Missing prefix for dhm output files, use option -outprefix. Exit(6).',6)
    if dekad is None:
        exitMessage('Undefined decad, use option -decad. Exit(10).',10)

    # check dirs existing
    if not os.path.isdir(indir):
        exitMessage('{0} does not exist or is not a directory. Exit(20).'.format(indir),20)
    if not os.path.isdir(outdir):
        if not os.path.exists(outdir):
            thisLogger.warning('Creating directory {0}. Continue.'.format(outdir))
            os.makedirs(outdir)
        else:
            exitMessage('{0} exists and is not a directory. Exit(21).'.format(outdir),21)
    if not os.path.isdir(inhist):
        exitMessage('{0} does not exist or is not a directory. Exit(20).'.format(), 20)

    # check files existing
    if not os.path.exists(realClimato):
        exitMessage('Climatology file {0} not found. Exit(22).'.format(realClimato),22)
    if not os.path.exists(maxRealClimato):
        exitMessage('Max Real Climatology file {0} not found. Exit(22).'.format(maxRealClimato),22)
    if not os.path.exists(realClimRMSAtMaxSST):
        exitMessage('climato RMS at max SST file {0} not found. Exit(22).'.format(realClimRMSAtMaxSST),22)

    # manage tmpdir
    if tmpdir is None:
        tmpdir = '{0}/tmp_{1}'.format(outdir, id_generator() )
    if not os.path.exists(tmpdir): os.makedirs(tmpdir)


    yearList=range(dekad, dekad+10)

    thisLogger.info('Processing information:')
    thisLogger.info('input directory: {0}'.format(indir))
    thisLogger.info('Input variable: {0}'.format(var))
    thisLogger.info('Model Climatology, path: {0}'.format(inhist))
    thisLogger.info('Model Climatology, file prefix: {0}'.format(modelClimatoRootName))
    thisLogger.info('Observation Climatology, file: {0}'.format(realClimato))
    thisLogger.info('Observation Climatology, max: {0}'.format(maxRealClimato))
    thisLogger.info('Observation Climatology, rms at max: {0}'.format(realClimRMSAtMaxSST))
    thisLogger.info('Output directory: {0}'.format(outdir))

    do_dhm(var, inhist, modelClimatoRootName, indir, sstRootName, realClimato, maxRealClimato, realClimRMSAtMaxSST, outdir, dhmRootName, yearList)


# end of file
