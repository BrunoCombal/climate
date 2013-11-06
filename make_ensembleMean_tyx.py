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

# ____________________________
def usage():
    textUsage='SYNOPSIS:\n\tmake_ensemble_Mean_tyx.py -v VARIABLE -path PATHIN -outdir PATHOUT\n\t-minVar MINVAL -maxVar MAXVAL\tn-model MODELLIST -startYear STARTYEAR -endYear ENDYEAR\n'
    textUsage=textUsage+'\tVARIABLE: a netcdf CMIP5 variable name, such as tos, zos, so, thetao;\n'
    textUsage=textUsage+'\tPATHIN: input data directory (does not support sub-directories);\n'
    textUsage=textUsage+'\tPATHOUT: output directory, created if does not exist;\n'
    textUsage=textUsage+'\tMINVAL: any value below minVar is considered as nodata;\n'
    textUsage=textUsage+'\tMAXVAL: any value above maxVar is considered as nodata;\n'
    textUsage=textUsage+'\tMODELLIST: a text file with a model name per name, the model name is used to select the files to process;\n'
    textUsage=textUsage+'\tSTARTYEAR: first year in the series of dates to process;\n'
    textUsage=textUsage+'\tENDYEAR: last year in the series of date to process'
    textUsage=textUsage+'In first place, the programme will average model output per model (if a model output has several rXiYpZ ensemble, they are averaged. Then, the averages are averaged to produce the ensemble mean.\n'
    textUsage=textUsage+'Averages are computed for each month of the year.\n'
    return textUsage
# ____________________________
def exitMessage(msg, exitCode='1'):
    print msg
    print
    print usage()
    sys.exit(exitCode)
# ____________________________
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for x in range(size))

def dateTime2Year(datetime):
    result=[]
    for ii in datetime:
        result.append(ii.year)
    return(numpy.array(result), sorted(set(result)))

def makeOutfileName(infile, outdir, prefix, year):
    return('{0}/{1}_{2}_{3}.nc'.format(outdir,prefix,os.path.basename(infile)[:-17], year))

def writeToFileOne(outfilename, data):
    if os.path.exists(outfilename): os.remove(outfilename)
    print 'writing {0}'.format(os.path.basename(outfilename))
    outfile = cdms2.open(outfilename, 'w')
    outfile.write(data)
    outfile.history='Created with '+__file__.encode('utf8')
    outfile.close()

def writeToFile(outfilename, data1, data2):
    if os.path.exists(outfilename): os.remove(outfilename)
    print 'writing {0}'.format(os.path.basename(outfilename))
    outfile = cdms2.open(outfilename, 'w')
    outfile.write(data1)
    outfile.write(data2)
    outfile.history='Created with '+__file__.encode('utf8')
    outfile.close()

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
# ___________________________
def updateCountersNew(accum, N, mini, maxi, data, minVar, maxVar, nodata=1.e20):

    print dir(accum)
    myOnes = MV2.ones(accum.array().shape)
    if accum is None:
        accum = data
        mini = data
        maxi = data
        print 'xxxx'
        print dir(accum)
        N = MV2.zeros( accum.shape , typecode='i') 
    else:

        accum = accum + MV2.masked_outside(data, minVar, maxVar)
        miniTmp = MV2.minimum(mini, accum)
        maxiTmp = MV2.maximum(maxi, accum)
        N = N 

    return [accum, N, mini, maxi]  
# ____________________________
def updateCounters(accum, N, mini, maxi, data, minVar, maxVar, nodata=1.e20):

    dim = numpy.squeeze(data[:]).shape

    if accum is None:
        accum = numpy.zeros(dim) + nodata
        N = numpy.zeros(dim) + nodata
        mini = data.copy()
        maxi = data.copy()

    wtadd = (data >= minVar ) * (data <= maxVar) * (accum < nodata) # add where not nodata
    wtreplace = (data >= minVar) * (data <= maxVar) * (accum >= nodata) # replace if no data
    wmax = (data >= maxi) * (data<nodata) * (data >= minVar) * (data <= maxVar)
    wmaxReplace = (mini >= nodata) * (data < nodata) * (data >= minVar)
    wmin = (data <= mini) * (data >= minVar) * ( data <= maxVar) * ( maxi < nodata )
    wminReplace = (mini >= nodata) * (data<nodata) * (data >= minVar)
    if wtadd.any():
        accum[wtadd] = accum[wtadd] + data[wtadd]
        N[wtadd] = N[wtadd] + 1 #numpy.ones(dim)
    if wtreplace.any():
        accum[wtreplace] = data[wtreplace]
        N[wtreplace] = 1 #numpy.ones(dim)
    if wmax.any():
        maxi[wmax] = data[wmax]
    if wmin.any():
        mini[wmin] = data[wmin]
    if wmaxReplace.any():
        maxi[wmaxReplace] = data[wmaxReplace]
    if wminReplace.any():
        mini[wminReplace] = data[wminReplace]

    return [accum, N, mini, maxi]
# ____________________________
def shiftGrid(infile, outfile, variable, latShift=0, lonShift=-280):

    if os.path.exists(outfile): os.remove(outfile)
    
    thisfile=cdms2.open(infile)

    thisLon = thisfile[variable].getLongitude()
    thisLon[:].data[:] = (thisLon[:].data[:] - lonShift)%360
    newvar = cdms2.createVariable(MV.array(thisfile[variable][:]), id=variable, fill_value=1.e20)
    newvar.setAxisList((thisfile[variable].getTime(), thisfile[variable].getLatitude(), thisLon))

    if os.path.exists(outfile): os.remove(outfile)
    outFile = cdms2.open(outfile,'w')
    outFile.write(newvar)
    thisfile.close()
    outFile.close()
    sys.exit(1)
# ___________________________
def do_regrid(variable, indir, lstInFile, outdir, stringBefore):

    createdFiles={}    
    nodata=1.e20

    # for netcdf3: set flag to 0
    cdms2.setNetcdfShuffleFlag(1)
    cdms2.setNetcdfDeflateFlag(1)
    cdms2.setNetcdfDeflateLevelFlag(3)

    if len(lstInFile)==0:
        exitMessage('Found no file to process, consider revising search pattern. Exit 6.',6)

    (newGrid, latAxis, lonAxis, lat_bnds, lon_bnds) = makeGrid()
    print lstInFile
    for fileName in lstInFile:
        print 'processing file ', fileName
        thisFile = cdms2.open(fileName)
        data = cdms2.createVariable(thisFile[variable])
        regrided = data.regrid(newGrid)
#        regrided = cdms2.createVariable(thisFile[variable]).regrid(newGrid)
        tmp = cdms2.createVariable(regrided, copyaxes=1, grid=newGrid)
        outfilename = '{0}/{1}{2}.nc'.format(outdir, stringBefore, os.path.basename(fileName))
        if os.path.exists(outfilename): os.remove(outfilename)
        outfile = cdms2.open(outfilename, 'w')
        outfile.write(regrided)
        outfile.close()

# ___________________________
# for a list of files: open all files, go from date 1 to date 2, compute avg for thisdate, save thisdate
# if a new grid is passed: regrid
def do_stats(variable, indir, lstInFile, outdir, stringBefore, outnameBase, minVar=-1.e20, maxVar=1.e20):
    
    if validYearList is None:
        exitMessage('List of years to process is undefined, edit code. Exit 5.',5)

    createdFiles={}   
    nodata=1.e20

    # for netcdf3: set flag to 0
    cdms2.setNetcdfShuffleFlag(1)
    cdms2.setNetcdfDeflateFlag(1)
    cdms2.setNetcdfDeflateLevelFlag(3)

    if lstInFile is None:
        exitMessage('Found no file to process, consider revising search pattern. Exit 6.',6)

    print lstInFile

    # open all files
    listFID=[]
    for ifile in lstInFile: listFID.append(cdms2.open(ifile, 'r'))
    
    # go through the list of dates, compute ensemble average
    for iyear in validYearList:
        print 'Processing year {0}'.format(iyear)
        for imonth in range(1,3):
            print '.',
            accumVar=None
            accumN=None
            mini=None
            maxi=None
            refGrid=None
            dims=None
            units=None
            for ifile in listFID:
                # get the actual file-date for this year/month (because day numbering vary, sometimes 1, 15 or 16...)
                thisTime = [ii for ii in ifile[variable].getTime().asComponentTime() if (ii.year==iyear and ii.month==imonth)] 
                if len(thisTime)==1:
                    if refGrid is None:
                        refGrid = ifile[variable].getGrid().toGenericGrid()
                        print refGrid.getLongitude(), refGrid.getLongitude().shape
                        dims = numpy.squeeze(ifile[variable].subRegion(time=thisTime[0])).shape
                        units= ifile[variable].units

                    [accumVar, accumN, mini, maxi]= updateCounters(accumVar, accumN, mini, maxi,
                                                                   numpy.array( ifile[variable].subRegion(time=thisTime[0])).ravel(),
                                                                   minVar, maxVar, nodata )
                        
            # compute average
            wtdivide = (accumN < nodata) * (accumN > 0)

            if wtdivide.any():
                accumVar[wtdivide] = accumVar[wtdivide] / accumN[wtdivide]
            accumVar = accumVar/accumN

            # compute std

            # save variables
            xx = MV2.array(accumVar.reshape(dims))
            newGrid = cdms2.createGenericGrid(refGrid.getLatitude(), refGrid.getLongitude(), latBounds=refGrid.getLatitude().getBounds(), lonBounds=refGrid.getLongitude().getBounds())
            xx.setGrid(newGrid)
            meanVar = cdms2.createVariable( xx, typecode='f', id='mean',  fill_value=nodata, attributes=dict(long_name='mean', units=units) )

            meanVar.setAxis(0, refGrid.getLatitude())
            meanVar.setAxis(1, refGrid.getLongitude())

            counter = cdms2.createVariable(accumN.reshape(dims), typecode='i', id='count', fill_value=nodata, attributes=dict(long_name='count', units='None') )
            miniVar = cdms2.createVariable(mini.reshape(dims), typecode='f', id='minimum', fill_value=nodata, attributes=dict(long_name='minimum', units=units) )
            maxiVar = cdms2.createVariable(maxi.reshape(dims), typecode='f', id='maximum', fill_value=nodata, attributes=dict(long_name='maximum', units=units) )
            outfilename = '{0}/{1}_{2}_{3}{4}.nc'.format(outdir, stringBefore, outnameBase, iyear, imonth )
            if os.path.exists(outfilename): os.remove(outfilename)
            outfile = cdms2.open(outfilename, 'w')
            outfile.write(meanVar)
            outfile.write(counter)
            outfile.write(miniVar)
            outfile.write(maxiVar)
            outfile.close()

            createdFiles['{0}{1}'.format(iyear,imonth)]= outfilename

    # close input files
    for ii in listFID: ii.close()

    return(createdFiles)
# ____________________________
# output: '{0}/{1}_{2}{3:02}.nc'.format(outdir, os.path.basename(ifile)[:-17], itime.year, itime.month)
def monthlyRegrid(variable, indir, outdir, validYearList=None, selectDATATYPE='Omon',selectMODEL='.*', selectRCP='.*', selectRIP='.*', selectTIMEFRAME='[0-9]*-[0-9]*'):

    #pattern=re.compile('.*_BNU-ESM_.*') # problem on this grid: longitude shift of about 280degrees, to be corrected with shiftGrid function

    # build file selection pattern
    # print '{0}_{1}_{2}_{3}_{4}_{5}.nc'.format(variable, selectDATATYPE, selectMODEL,selectRCP,selectRIP,selectTIMEFRAME)
    pattern = re.compile('{0}_{1}_{2}_{3}_{4}_{5}.nc'.format(variable, selectDATATYPE, selectMODEL,selectRCP,selectRIP,selectTIMEFRAME))

    nodata=1.e20
    minVar=273.15-40
    maxVar=273.15+100
    if validYearList is None:
        validYearList = (2006, 2030,2031,2032,2033,2034,2035, 2036, 2037, 2038, 2039, 2050,2051,2052,2053,2054,2055, 2056, 2057, 2058, 2059)
    maskpattern = re.compile('.*EC-EARTH.*') # nodate was set to 273.15

    # for netcdf3: set flag to 0
    cdms2.setNetcdfShuffleFlag(1)
    cdms2.setNetcdfDeflateFlag(1)
    cdms2.setNetcdfDeflateLevelFlag(3)
    (referenceGrid, latAxis, lonAxis, latBounds, lonBounds) = makeGrid()

    lstInFile=[f for f in glob.glob('{0}/*.nc'.format(indir)) if (os.stat(f).st_size and pattern.match(os.path.basename(f) ) ) ]
    #lstInFile = [f for f in os.listdir() if re.search(pattern, f)]
    if lstInFile is None:
        print 'Error 1. Exit'
        sys.exit(1)

    # print lstInFile

    sys.exit(1)

    # loop over files, extract month, regrid it
    for ifile in lstInFile:
        print ifile
        thisfile=cdms2.open(ifile)
        var=cdms2.createVariable(thisfile[variable]) # needed if one want to update the mask (else read only variable)

        if maskpattern.match(os.path.basename(ifile)):
            print 'Correcting mask for ',os.path.basename(ifile)
            # mask pixels where there is no change (abs(change)<1.e4)
            refshape=var.shape
            tmp = numpy.reshape(var, (refshape[0], refshape[1] * refshape[2]) )
            wtnodata = (tmp.max(axis=0) - tmp.min(axis=0)) < 0.001
            print wtnodata.shape
            if wtnodata.any():
                for ii in range(refshape[0]):
                    tmp[ii, wtnodata] = nodata
                var[:] = numpy.reshape(tmp, refshape)

        # note: usually time reference is stored in the file name
        startDate = int(ifile[-16:-10])
        startYear = int(ifile[-16:-12])
        if startYear <= numpy.max(validYearList):
            itimePos=0

            regriddedAll = var.regrid(referenceGrid)
            for itime in var.getTime().asComponentTime():
                if itime.year in validYearList:

                    outfilename='{0}/{1}_{2}{3:02}.nc'.format(outdir, os.path.basename(ifile)[:-17], itime.year, itime.month)
                    #print outfilename

                    if not os.path.exists(outfilename): # if it already exists, don't do it again: save time
                        regridded = regriddedAll.subRegion(time=itime)
                        print 'regridded !'
                        # thisVar = var.subRegion(time=itime)
                        # regridded = thisVar.regrid(referenceGrid)
                        # recodedVar=MV.masked_where( thisVar <= minVar, thisVar, copy=1) # recoding of no-data was the first approach to 
                        # #regridded = recodedVar.regrid(referenceGrid)   # build EC-EARTH mask; now done in regridding function
            
                        temp1 = cdms2.createVariable(regridded, typecode='f',\
                                                         id=variable, fill_value=nodata,\
                                                         grid=referenceGrid, copyaxes=0,\
                                                         attributes=dict(long_name=var.long_name, units=var.units) )
                    
                        writeToFileOne(outfilename, temp1)
                itimePos = itimePos+1
 
        thisfile.close()

#_______________________________________
# monthly average: average models for the same YYYYMM
# compute model ensemble mean, max, min, std.
def monthlyAvg(variable, indir, outdir, minYear=2006, maxYear=2050, select='*'):
    nodata=1.e20
    minVar=273.15 - 40
    maxVar=273.15 + 100
    unitsAvg=None
    # assume data are aligned
    pattern=re.compile('.*_BNU-ESM_.*') # problem on this grid (use shiftGrid to create a new version, discard this one).
    # maskpattern = re.compile('.*EC-EARTH.*') # nodate was set to 273.15
    
    # for netcdf3: set flags to 0
    cdms2.setNetcdfShuffleFlag(1)
    cdms2.setNetcdfDeflateFlag(1)
    cdms2.setNetcdfDeflateLevelFlag(3)

    print minYear, maxYear

    dateList=[]
    for iyear in range(minYear, maxYear+1):
        for imonth in range(1,13):
            dateList.append('{0}{1:02}'.format(iyear,imonth))

    for idate in dateList:
        print 'processing date ',idate
        # get list of files for this iyear, excluding one file:
        lstFiles = [f for f in glob.glob(indir+'/{0}_*{2}*_{1}.nc'.format(variable, idate,select)) if not pattern.match(f) ]
        print indir+'/{0}_*{2}*_{1}.nc'.format(variable, idate,select)

        if len(lstFiles) > 1:
            print 'Model ensemble mean for date {0} with {1} files'.format(idate, len(lstFiles))
            # accumulate data
            accumVar=None
            accumN=None

            for iFile in lstFiles:
                print 'processing file ', iFile
                thisFile = cdms2.open(iFile)
                dimVar = numpy.squeeze(thisFile[variable][:]).shape # remove time-single dimension if exists
                thisVar = numpy.ravel(thisFile[variable][:])
                
                # recode thisVar if needed
                #if maskpattern.match(os.path.basename(iFile)):
                #    wtrecode = numpy.abs( thisVar - 273.15 ) < 0.001
                #    if wtrecode.any():
                #        thisVar[wtrecode] = nodata
                #        print 'date recoded for ',os.path.basename(iFile)

                if accumVar is None:
                    accumVar  = numpy.zeros( dimVar[0]*dimVar[1] ) + nodata
                    accumN    = numpy.zeros( dimVar[0]*dimVar[1] )
                    unitsAvg  = thisFile[variable].units
                    oneMatrix = numpy.ones(dimVar[0]*dimVar[1])
                    maxEnsemble = thisVar.copy()
                    minEnsemble = thisVar.copy()

                # add to accumVar if accumVar is not no-data, and incoming data are in the range
                wtadd = (thisVar >= minVar ) * (thisVar <= maxVar) * (accumVar < nodata)
                # if the value in accumVar is no-data, replace it.
                wtreplace = (thisVar >= minVar ) * (thisVar <= maxVar) * ( accumVar >= nodata)
                # min, max
                wmax = (thisVar >= maxEnsemble) * (thisVar < nodata) * (thisVar >= minVar) * (thisVar <= maxVar)
                wmaxReplace = (maxEnsemble >= nodata) * (thisVar < nodata) * (thisVar >= minVar)
                wmin = (thisVar <= minEnsemble) * (thisVar >= minVar) * (thisVar <= maxVar) * (maxEnsemble < nodata)
                wminReplace = (minEnsemble >= nodata) * (thisVar < nodata) * (thisVar >= minVar)
                if wtadd.any():
                    accumVar[wtadd] = accumVar[wtadd] + thisVar[wtadd]
                    accumN[wtadd] = accumN[wtadd] + oneMatrix[wtadd]
                if wtreplace.any():
                    accumVar[wtreplace] = thisVar[wtreplace]
                    accumN[wtreplace] = oneMatrix[wtreplace]
                if wmax.any():
                    maxEnsemble[wmax] = thisVar[wmax]
                if wmin.any():
                    minEnsemble[wmin] = thisVar[wmin]
                if wmaxReplace.any():
                    maxEnsemble[wmaxReplace] = thisVar[wmaxReplace]
                if wminReplace.any():
                    minEnsemble[wminReplace] = thisVar[wminReplace]

                thisFile.close()

            # now compute the average, where accumN is not 0
            wnz = accumN > 0
            average = numpy.zeros(dimVar[0] * dimVar[1]) + nodata
            if wnz.any():
                average[wnz] = accumVar[wnz] / accumN[wnz]

            # and let's compute the std
            std = numpy.zeros(dimVar[0] * dimVar[1]) + nodata
            stdN = numpy.zeros(dimVar[0] * dimVar[1])
            for iFile in lstFiles:
                thisFile = cdms2.open(iFile)
                thisVar = numpy.ravel(thisFile[variable][:])
                wtadd = (thisVar < nodata ) * (average < nodata ) * (thisVar >= minVar) * (thisVar <= maxVar) # average should be clean, no need to implement a 'replace'
                if wtadd.any():
                    std[wtadd] = (average[wtadd] - thisVar[wtadd]) * (average[wtadd] - thisVar[wtadd])
                    stdN[wtadd] = stdN[wtadd] + 1.0
                thisFile.close()

            wtstd = (stdN > 0) * (std < nodata)
            std[wtstd] = numpy.sqrt( std[wtstd]/stdN[wtstd] )

            # save to disk
            outfilename='{0}/modelmean_{1}_{2}.nc'.format(outdir, variable, idate)
            (referenceGrid, latAxis, lonAxis, latBounds, lonBounds) = makeGrid()
            avgOut = cdms2.createVariable(numpy.reshape(average,dimVar), typecode='f', id=variable, fill_value=1.e20, grid=referenceGrid, copyaxes=1, attributes=dict(long_name='model average for {0} at date {1}'.format(variable, idate), units=unitsAvg))
            avgOut.setAxisList((latAxis, lonAxis))

            accumOut = cdms2.createVariable(numpy.reshape(accumN,dimVar), typecode='i', id='count_{0}'.format(variable), fill_value=1.e20, grid=referenceGrid, copyaxes=1, attributes=dict(long_name='count of valid for {0} at date {1}'.format(variable, idate), units=None))
            accumOut.setAxisList((latAxis, lonAxis))

            maxEns = cdms2.createVariable(numpy.reshape(maxEnsemble,dimVar), typecode='f', id='max {0}'.format(variable), fill_value=1.e20, grid=referenceGrid, copyaxes=1, attributes=dict(long_name='max ensemble for {0} at date {1}'.format(variable, idate), units=unitsAvg))
            maxEns.setAxisList((latAxis, lonAxis))

            minEns = cdms2.createVariable(numpy.reshape(minEnsemble,dimVar), typecode='f', id='min {0}'.format(variable), fill_value=1.e20, grid=referenceGrid, copyaxes=1, attributes=dict(long_name='min ensemble for {0} at date {1}'.format(variable, idate), units=unitsAvg))
            minEns.setAxisList((latAxis, lonAxis))

            stdVar = cdms2.createVariable(numpy.reshape(std,dimVar), typecode='f', id='std_{0}'.format(variable), fill_value=1.e20, grid=referenceGrid, copyaxes=1, attributes=dict(long_name='model std for {0} at date {1}'.format(variable, idate), units=unitsAvg))
            stdVar.setAxisList((latAxis, lonAxis))

            if os.path.exists(outfilename): os.remove(outfilename)
            print 'saving to file ', outfilename
            outfile = cdms2.open(outfilename, 'w')
            outfile.write(avgOut)
            outfile.write(accumOut)
            outfile.write(minEns)
            outfile.write(maxEns)
            outfile.write(stdVar)
            outfile.history='Created with '+__file__.encode('utf8')
            outfile.close()
#___________________________
if __name__=="__main__":

    print 'To make this script properly work, ensure to source cdat setup file first (source /usr/local/uvcdat/VERSION/bin/setup_runtime.sh)'
    print

    variable = None
    indir = None
    tmpdir = None
    outdir = None
    # model list: model name (for output), regexp is based on .*_MODEL_.*
    #['ACCESS1-0', 'ACCESS1-3', 'bcc-csm1-1', 'bcc-csm1-1-m', 'BNU-ESM','CanESM2', 'CCSM4', 'CESM1-BGC', 'CESM1-CAM5', 'CESM1-WACCM','CMCC-CESM', 'CMCC-CM', 'CMCC-CMS', 'CNRM-CM5', 'CSIRO-Mk3-6-0','EC-EARTH', 'FIO-ESM', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M','GISS-E2-H', 'GISS-E2-H-CC', 'GISS-E2-R', 'GISS-E2-R-CC', 'HadGEM2-AO', 'HadGEM2-CC', 'HadGEM2-ES', 'inmcm4', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR', 'MIROC5', 'MIROC-ESM', 'MPI-ESM-LR','MPI-ESM-MR', 'MRI-CGCM3', 'NorESM1-M', 'NorESM1-ME']
    modelListFile=None
    startYear=None
    endYear=None

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
        elif arg=='-minVar':
            ii = ii + 1
            minVar = float(sys.argv[ii])
        elif arg == '-maxVar':
            ii = ii + 1
            maxVar = float(sys.argv[ii])
        elif arg =='-modellist':
            ii = ii + 1
            modelListFile = sys.argv[ii]
        elif arg=='-startyear':
            ii = ii + 1
            startYear = int(sys.argv[ii])
        elif arg=='-endyear':
            ii = ii + 1
            endYear = int(sys.argv[ii]) + 1
        ii = ii + 1

    if variable is None:
        exitMessage('Missing variable name, use option -v. Exit(1).', 1)
    if indir is None:
        exitMessage('Missing input directory, use option -path. Exit(2).',2)
    if outdir is None:
        exitMessage('Missing output directory, use option -outdir. Exit(3).', 3)
    if modelListFile is None:
        exitMessage('Missing a model list file, use option -modellist. Exit(12).',12)
    if startYear is None:
        exitMessage('Please define a starting year, use option -startyear. Exit(13).',13)
    if endYear is None:
        exitMessage('Please define an ending year, use option -endyear. Exit(14).',14)

    if tmpdir is None:
        tmpdir = '{0}/tmp_{1}'.format(outdir, id_generator() )

    if not os.path.exists(tmpdir): os.makedirs(tmpdir)
    
    # models list
    modelList=[]
    try:
        with open(modelListFile,"r") as f:
            for textLine in f:
                thisStr = textLine.replace(" ","").replace('\n','')
                if not (thisStr==""):
                    modelList.append( thisStr )
    except IOError as e:
        exitMessage('I/O Error {1} while processing text file {0}:{2}. Exit(10).'.format(modelListFile, e.errno, e.strerror), 10)
    except:
        exitMessage('Unexpected error while processing text file {0}. Exit(11).'.format(modeListFile), 11)

    validYearList=range(startYear, endYear)
    if len(validYearList)==0:
        exitMessage('No date to process, startYear={0}, endYear{1}. Exit(20).'.format(startYear, endYear),20)

    processedFiles=[]
    for thisModel in modelList:

        pattern=re.compile('{0}_{1}_{2}_{3}_{4}_{5}.nc'.format(variable, 'Omon', thisModel, 'rcp85', 'r.*i.*p.*', '.*') )
        lstInFile=[f for f in glob.glob('{0}/*.nc'.format(indir)) if (os.stat(f).st_size and pattern.match(os.path.basename(f) ) ) ]

        regridedFiles = do_regrid(variable, indir, lstInFile, tmpdir, 'regrid_')

#        pattern=re.compile('{0}_{1}_{2}_{3}_{4}_{5}.nc'.format(variable, 'Omon', thisModel, 'rcp85', 'r.*i.*p.*', '.*') )
#        lstInFile=[f for f in glob.glob('{0}/*.nc'.format(indir)) if (os.stat(f).st_size and pattern.match(os.path.basename(f) ) ) ]

        #thisModelFiles = do_stats(variable, indir, lstInFile, tmpdir, 'stats_', '{0}_{1}'.format(thisModel, 'rcp85', 1, 300) )


#        processedFiles.append(thisModelFiles)
#        print processedFiles

    print '____________'
    print processedFiles

#    do_regrid(new_grid, thismodelfiles)

    # monthlyRegrid(variable, indir, tmpdirNew, validYearList, 'Omon', 'MPI.*', 'rcp85', '.*', '.*')
        
    # note: to correctly compute the ensemble means, consider deleting unwanted files!!!

    # for historical data set, selecting only r1i1p1:
    # monthlyAvg(variable, tmpdirHist, outdirHist, 1959, 2005, 'historical_r1i1p1')
    # for projections data set, selecting only r1i1p1:

    # file naming convention is:
    # VARIABLE_Omon_MODEL_RIP_date1_date2.nc
    # monthlyAvg(variable, tmpdirNew, outdir, 2006, 2059, 'Omon','.*','rcp85','r[0-9]*i[0-9]*p[0-9]*','.*')


