#!/usr/bin/env python

# to run the script with the correct version of uvcdat:
#  source /usr/local/uvcdat/1.2.0/bin/setup_cdat.sh

import cdms2
from cdms2 import MV
import numpy
import glob
import sys
import os
from os import path
import re
import string

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
#_____________________________________
def avgYearRegrid(variable, indir, tmpdir, minYear=2006, maxYear=2050):
    nodata=1.e20
    minVar=274-40
    maxVar=274+100

    # for netcdf3:
    cdms2.setNetcdfShuffleFlag(0)
    cdms2.setNetcdfDeflateFlag(0)
    cdms2.setNetcdfDeflateLevelFlag(0)
    (referenceGrid, latAxis, lonAxis, latBounds, lonBounds)=makeGrid()

    lstInFile=[f for f in glob.glob('{0}/*.nc'.format(indir)) if os.stat(f).st_size ]
    if lstInFile is None:
        print 'Error 1. Exit'
        sys.exit(1)

    # first: compute annual means: there must be 12 values/means
    for ifile in lstInFile:
        thisfile=cdms2.open(ifile)
        var=thisfile[variable]
        # note: usually time reference is stored in the file name
        startDate=int(ifile[-16:-10])
        startYear=int(ifile[-16:-12])
        if startYear <= maxYear:
            listYears, uniqYears=dateTime2Year(var.getTime().asComponentTime())
            listIndex=numpy.array(range(len(listYears)))

            for iuniq in uniqYears:
                wtp = (listYears == iuniq) * (iuniq >= minYear) * (iuniq <=maxYear)
                if len(listIndex[wtp]) == 12:
                    thisList=listIndex[wtp]
                    print 'year {0}, file {1}, averaging for {2}-{3}'.format(iuniq, os.path.basename(ifile),thisList[0],thisList[-1])
                    recodedVar=MV.masked_where( var<=minVar, var, copy=1)
                    thisaverage = MV.average(recodedVar[ thisList[0]:thisList[-1]+1, :, :], axis=0)
                    # let's regrid
                    regridded = thisaverage.regrid(referenceGrid)
                    
                    temp1 = cdms2.createVariable(regridded, typecode='f', id=variable, fill_value=1.e20, grid=referenceGrid, copyaxes=0, attributes=dict(long_name=var.long_name, units=var.units) )
                    outfilename=makeOutfileName(ifile, tmpdir, 'avg',iuniq)

                    writeToFileOne(outfilename, temp1)
 
        thisfile.close()

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

# ____________________________
# output: '{0}/{1}_{2}{3:02}.nc'.format(outdir, os.path.basename(ifile)[:-17], itime.year, itime.month)
def monthlyRegrid(variable, indir, outdir, validYearList=None, select='*'):

    pattern=re.compile('.*_BNU-ESM_.*') # problem on this grid: longitude shift of about 280degrees, to be corrected with shiftGrid function

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
    (referenceGrid, latAxis, lonAxis, latBounds, lonBounds)=makeGrid()

    lstInFile=[f for f in glob.glob('{0}/tos*_{1}_*.nc'.format(indir, select)) if os.stat(f).st_size ]
    if lstInFile is None:
        print 'Error 1. Exit'
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

#____________________________
def yearlyAvg(variable, indir, outdir, minYear=2006, maxYear=2050):
    nodata=1.e20
    minVar=274-40
    maxVar=274+100
    unitsAvg=None
    # assume data are aligned
    pattern=re.compile('.*_BNU-ESM_.*')
    
    for iyear in range(minYear, maxYear+1):
        # get list of files for this iyear
        # let's exclude one file:
        lstFiles = [f for f in glob.glob(indir+'/avg_{0}_*_{1}.nc'.format(variable, iyear)) if not pattern.match(f) ]
    
        if len(lstFiles) > 1:
            print 'Averaging year {0} with {1} files'.format(iyear, len(lstFiles))
            # accumulate data
            accumVar=None
            accumN=None

            for iFile in lstFiles:
                #print 'processing file ', iFile
                thisFile = cdms2.open(iFile)
                dimVar=thisFile[variable][:].shape
                thisVar = numpy.ravel(thisFile[variable][:])
                if accumVar is None:
                    print 'INITIALISATION !'
                    accumVar  = numpy.zeros( dimVar[0]*dimVar[1] ) + nodata
                    accumN    = numpy.zeros( dimVar[0]*dimVar[1] ) 
                    unitsAvg  = thisFile[variable].units
                    oneMatrix = numpy.ones(dimVar[0]*dimVar[1])

                # add to accumVar if accumVar is not no-data, and incoming data are in the range
                wtadd = (thisVar >= minVar ) * (thisVar <= maxVar) * (accumVar < nodata)
                # if the value in accumVar is no-data, replace it.
                wtreplace = (thisVar >= minVar ) * (thisVar <= maxVar) * ( accumVar >= nodata)
                if wtadd.any():
                    #print 'adder ',thisVar[wtadd].max(), accumVar[wtadd].max(), accumN[wtadd].max()
                    accumVar[wtadd] = accumVar[wtadd] + thisVar[wtadd]
                    accumN[wtadd] = accumN[wtadd] + oneMatrix[wtadd]
                    #print 'adder ',thisVar[wtadd].max(), accumVar[wtadd].max(), accumN[wtadd].max()
                if wtreplace.any():
                    #print 'replace'
                    accumVar[wtreplace] = thisVar[wtreplace]
                    accumN[wtreplace] = oneMatrix[wtreplace]
                
                thisFile.close()

            # now compute the average, where accumN is not 0
            wnz = accumN > 0
            average = numpy.zeros(dimVar[0] * dimVar[1]) + nodata
            if wnz.any():
                average[wnz] = accumVar[wnz] / accumN[wnz]

            # save to disk
            # for netcdf3: 
            cdms2.setNetcdfShuffleFlag(0)
            cdms2.setNetcdfDeflateFlag(0)
            cdms2.setNetcdfDeflateLevelFlag(0)

            outfilename='{0}/avgyearly_{1}_{2}.nc'.format(outdir, variable, iyear)
            referenceGrid=makeGrid()
            avgOut = cdms2.createVariable(numpy.reshape(average,dimVar), typecode='f', id=variable, fill_value=1.e20, grid=referenceGrid, copyaxes=0, attributes=dict(long_name='average {0} for {1}'.format(variable, iyear), units=unitsAvg))
            accumOut = cdms2.createVariable(numpy.reshape(accumN,dimVar), typecode='i', id='count_{0}'.format(variable), fill_value=1.e20, grid=referenceGrid, copyaxes=0, attributes=dict(long_name='count of valid {0} for {1}'.format(variable, iyear), units=None))
            writeToFile(outfilename, avgOut, accumOut)

#_______________________________________
# monthly average: average models for the same YYYYMM
# compute model ensemble mean, max, min, std.
def monthlyAvgOLD(variable, indir, outdir, minYear=2006, maxYear=2050):
    print 'XXX'
    sys.exit(0)
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

    print 'here'
    dateList=[]
    for iyear in range(minYear, maxYear+1):
        for imonth in range(1,13):
            dateList.append('{0}{1:02}'.format(iyear,imonth))

    print dateList

    for idate in dateList:
        # get list of files for this iyear, excluding one file:
        lstFiles = [f for f in glob.glob(indir+'/{0}_*_{1}.nc'.format(variable, idate)) if not pattern.match(f) ]
        print 'for date ',idate, 'looking for files in ','{2}/{0}_*_{1}.nc'.format(variable, idate,indir)

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
            outfile = cdms2.open(outfilename, 'w')
            outfile.write(avgOut)
            outfile.write(accumOut)
            outfile.write(minEns)
            outfile.write(maxEns)
            outfile.write(stdVar)
            outfile.history='Created with '+__file__.encode('utf8')
            outfile.close()

#___________________________
def buildTimeSeries(indir, outdir, variableList, rootname, startYear, endYear):

    # for netcdf3: set flags to 0
    cdms2.setNetcdfShuffleFlag(1)
    cdms2.setNetcdfDeflateFlag(1)
    cdms2.setNetcdfDeflateLevelFlag(3)

    dateList=[]
    for iyear in range(startYear, endYear+1):
        for imonth in range(1, 12+1):
            dateList.append('{0}{1:02}'.format(iyear, imonth))


    outfilename='{0}/ts_{1}.nc'.format(outdir,variableList[0])
    if os.path.exists(outfilename): os.remove(outfilename)
    outfile = cdms2.open(outfilename, 'w')

    for idate in dateList:
        fname='{0}/{1}_{2}_{3}.nc'.format(indir,rootname, variable, idate)
        print 'appending ', fname
        thisfile = cdms2.open(fname)
        for ivar in variableList:
            print 'XX ',ivar
            outfile.write(thisfile[ivar])
        thisfile.close()

    outfile.close()

# __________________________
def do_precision(indir, fnameStart, var, yearList, outfile):

   # for netcdf3: set flag to 0
    cdms2.setNetcdfShuffleFlag(1)
    cdms2.setNetcdfDeflateFlag(1)
    cdms2.setNetcdfDeflateLevelFlag(3)
    (referenceGrid, latAxis, lonAxis, latBounds, lonBounds)=makeGrid()

    lstInFile=glob.glob('{0}/*.nc'.format(indir))
    if lstInFile is None:
        print 'Error 1. Exit'
        sys.exit(1)

    # loop over files, extract month, regrid it
    for ifile in lstInFile:
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
                    regridded = regriddedAll.subRegion(time=itime)
                    #thisVar = var.subRegion(time=itime)
                    #regridded = thisVar.regrid(referenceGrid)
                    #recodedVar=MV.masked_where( thisVar <= minVar, thisVar, copy=1) # recoding of no-data was the first approach to 
                    ##regridded = recodedVar.regrid(referenceGrid)   # build EC-EARTH mask; now done in regridding function
            
                    temp1 = cdms2.createVariable(regridded, typecode='f',\
                                                     id=variable, fill_value=nodata,\
                                                     grid=referenceGrid, copyaxes=0,\
                                                     attributes=dict(long_name=var.long_name, units=var.units) )
                    outfilename='{0}/{1}_{2}{3:02}.nc'.format(outdir, os.path.basename(ifile)[:-17], itime.year, itime.month)
                    print outfilename
                    writeToFileOne(outfilename, temp1)
                itimePos = itimePos+1
 
        thisfile.close()

#____________________________
def yearlyAvg(variable, indir, outdir, minYear=2006, maxYear=2050):
    nodata=1.e20
    minVar=274-40
    maxVar=274+100
    unitsAvg=None
    # assume data are aligned
    pattern=re.compile('.*_BNU-ESM_.*')
    
    for iyear in range(minYear, maxYear+1):
        # get list of files for this iyear
        # let's exclude one file:
        lstFiles = [f for f in glob.glob(indir+'/avg_{0}_*_{1}.nc'.format(variable, iyear)) if not pattern.match(f) ]
    
        if len(lstFiles) > 1:
            print 'Averaging year {0} with {1} files'.format(iyear, len(lstFiles))
            # accumulate data
            accumVar=None
            accumN=None

            for iFile in lstFiles:
                #print 'processing file ', iFile
                thisFile = cdms2.open(iFile)
                dimVar=thisFile[variable][:].shape
                thisVar = numpy.ravel(thisFile[variable][:])
                if accumVar is None:
                    print 'INITIALISATION !'
                    accumVar  = numpy.zeros( dimVar[0]*dimVar[1] ) + nodata
                    accumN    = numpy.zeros( dimVar[0]*dimVar[1] ) 
                    unitsAvg  = thisFile[variable].units
                    oneMatrix = numpy.ones(dimVar[0]*dimVar[1])

                # add to accumVar if accumVar is not no-data, and incoming data are in the range
                wtadd = (thisVar >= minVar ) * (thisVar <= maxVar) * (accumVar < nodata)
                # if the value in accumVar is no-data, replace it.
                wtreplace = (thisVar >= minVar ) * (thisVar <= maxVar) * ( accumVar >= nodata)
                if wtadd.any():
                    #print 'adder ',thisVar[wtadd].max(), accumVar[wtadd].max(), accumN[wtadd].max()
                    accumVar[wtadd] = accumVar[wtadd] + thisVar[wtadd]
                    accumN[wtadd] = accumN[wtadd] + oneMatrix[wtadd]
                    #print 'adder ',thisVar[wtadd].max(), accumVar[wtadd].max(), accumN[wtadd].max()
                if wtreplace.any():
                    #print 'replace'
                    accumVar[wtreplace] = thisVar[wtreplace]
                    accumN[wtreplace] = oneMatrix[wtreplace]
                
                thisFile.close()

            # now compute the average, where accumN is not 0
            wnz = accumN > 0
            average = numpy.zeros(dimVar[0] * dimVar[1]) + nodata
            if wnz.any():
                average[wnz] = accumVar[wnz] / accumN[wnz]

            # save to disk
            # for netcdf3: 
            cdms2.setNetcdfShuffleFlag(0)
            cdms2.setNetcdfDeflateFlag(0)
            cdms2.setNetcdfDeflateLevelFlag(0)

            outfilename='{0}/avgyearly_{1}_{2}.nc'.format(outdir, variable, iyear)
            referenceGrid=makeGrid()
            avgOut = cdms2.createVariable(numpy.reshape(average,dimVar), typecode='f', id=variable, fill_value=1.e20, grid=referenceGrid, copyaxes=0, attributes=dict(long_name='average {0} for {1}'.format(variable, iyear), units=unitsAvg))
            accumOut = cdms2.createVariable(numpy.reshape(accumN,dimVar), typecode='i', id='count_{0}'.format(variable), fill_value=1.e20, grid=referenceGrid, copyaxes=0, attributes=dict(long_name='count of valid {0} for {1}'.format(variable, iyear), units=None))
            writeToFile(outfilename, avgOut, accumOut)

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
def buildTimeSeries(indir, outdir, variableList, rootname, startYear, endYear):

    # for netcdf3: set flags to 0
    cdms2.setNetcdfShuffleFlag(1)
    cdms2.setNetcdfDeflateFlag(1)
    cdms2.setNetcdfDeflateLevelFlag(3)

    dateList=[]
    for iyear in range(startYear, endYear+1):
        for imonth in range(1, 12+1):
            dateList.append('{0}{1:02}'.format(iyear, imonth))


    outfilename='{0}/ts_{1}.nc'.format(outdir,variableList[0])
    if os.path.exists(outfilename): os.remove(outfilename)
    outfile = cdms2.open(outfilename, 'w')

    for idate in dateList:
        fname='{0}/{1}_{2}_{3}.nc'.format(indir,rootname, variable, idate)
        print 'appending ', fname
        thisfile = cdms2.open(fname)
        for ivar in variableList:
            print 'XX ',ivar
            outfile.write(thisfile[ivar])
        thisfile.close()

    outfile.close()

# __________________________
def do_precisionDHM(indir, fnameStart, var, yearList, outfile):

    # for netcdf3: set flag to 0
    cdms2.setNetcdfShuffleFlag(0)
    cdms2.setNetcdfDeflateFlag(0)
    cdms2.setNetcdfDeflateLevelFlag(0)

    refYear=2006
    nodata = 1.e20

    lstFile = glob.glob('{0}/{1}_*.nc'.format(indir, fnameStart))
    
    # list of dates
    dateList=[]
    dateListRef=[]
    sumStd = None
    thisGrid=None
    thisUnits=None
    for iy in yearList:
        for im in range(1, 12+1):
            dateList.append('{0}{1:02}'.format(iy, im))
            dateListRef.append(['{0}{1:02}'.format(refYear, im),\
                                    '{0}{1:02}'.format(refYear, (im-1-1)%12+1),\
                                    '{0}{1:02}'.format(refYear, (im-2-1)%12+1)]) 

    for (idate, idateRef) in zip(dateList, dateListRef):
        print 'Summing precision for dates ', idate, idateRef
        print '{0}/{1}_{2}.nc'.format( indir, fnameStart, idate)
        thisFile     = cdms2.open('{0}/{1}_{2}.nc'.format( indir, fnameStart, idate) )
        if var not in thisFile.variables.keys():
            print 'Variable {0} not found in {1}_{2}.nc. Exit 1.'.format(var, fnameStart, idate)
            sys.exit(1)
        if thisGrid is None: thisGrid = thisFile[var].getGrid()
        if thisUnits is None: thisFile[var].units
        thisFileRef0 = cdms2.open('{0}/{1}_{2}.nc'.format( indir, fnameStart, idateRef[0]) )
        thisFileRef1 = cdms2.open('{0}/{1}_{2}.nc'.format( indir, fnameStart, idateRef[1]) )
        thisFileRef2 = cdms2.open('{0}/{1}_{2}.nc'.format( indir, fnameStart, idateRef[2]) )
        
        # add everything
        if sumStd is None:
            sumStd = MV.array(thisFile[var])
        sumStd[:] = sumStd[:] + MV.array(thisFileRef0[var][:]) + MV.array(thisFileRef1[var][:]) + MV.array(thisFileRef2[var][:])

        thisFile.close()
        thisFileRef0.close()
        thisFileRef1.close()
        thisFileRef2.close()

    # save
    if os.path.exists(outfile): os.remove(outfile)
    # average the stds
    sumStd[:] = sumStd[:]/len(dateList)    
    # then correct the mask
    wmask = MV.array(thisFile[var]) >= nodata
    if wmask.any():
        print 're- masking'
        sumStd[wmask] = nodata
    dataOut = cdms2.createVariable(sumStd, typecode='f', id='average_std', \
                                       fill_value=nodata, grid=thisGrid, \
                                       copaxes=1, \
                                       attributes=dict(long_name='average of std is an estimate of precision', units=thisUnits ))
    outFile = cdms2.open(outfile,'w')
    outFile.write(sumStd)
    outFile.close()

#___________________________
if __name__=="__main__":

    variable='tos'
    #indir='/data/cmip5/rcp/rcp8.5/tos/'
    indir='/databis/cmip5_bis/rcp/rcp4.5/tos/'
    dirHistorical='/data/cmip5/rcp/tos_historical/'
    tmpdir='/home/bruno/Documents/tmp/tos/'
    #tmpdirNew='/home/bruno/Documents/tmp/tos_monthly/'
    tmpdirNew='/home/bruno/Documents/tmp/tos4.5_monthly/'
    tmpdirHist='/home/bruno/Documents/tmp/tos_hist/'
    #outdir='/data/cmip5/rcp/rcp8.5/tos_ensemble/'
    outdir='/data/cmip5/rcp/rcp8.5/tos4.5_ensemble/'
    outdirHist='/data/cmip5/rcp/rcp8.5/toshist_ensemble/'

    if not os.path.exists(tmpdir): os.makedirs(tmpdir)
    if not os.path.exists(tmpdirNew): os.makedirs(tmpdirNew)
    if not os.path.exists(tmpdirHist): os.makedirs(tmpdirHist)

    processing='climato'

    if processing == 'tos':
        # avgYearRegrid(variable, indir, tmpdir, 2006, 2010)
        # yearlyAvg(variable, tmpdir, outdir, 2006, 2050)

        shiftGrid('/data/cmip5/rcp/rcp8.5/tos/tos_Omon_BNU-ESM_rcp85_r1i1p1_200601-210012.nc','/data/cmip5/rcp/rcp8.5/tos/tos_Omon_NEWbnu-esm_rcp85_r1i1p1_200601-210012.nc','tos',0,-280)

        #    monthlyRegrid(variable, indir, tmpdirNew)
        #    listYear=[]
        #    for ii in range(2010, 2029 + 1): listYear.append(ii)
        #    for ii in range(2040, 2049 + 1): listYear.append(ii)
        #    monthlyRegrid(variable, indir, tmpdirNew, listYear)
        #    monthlyAvg(variable, tmpdirNew, outdir, 2006, 2059)
        
        #    buildTimeSeries(outdir, outdir, [variable,'min tos','max tos','std_tos'], 'modelmean', 2010,2059)
        
        # computing final precision
        #    do_precisionDHM(outdir, 'modelmean_tos', 'std_tos', range(2030, 2040), '{0}/{1}'.format(outdir, 'dhm_precision_2030.nc'))
    elif processing == 'climato':
        #validYearList = (1850, 1851, 1852, 1853, 1854, 1855, 1856, 1857, 1858, 1859, 1860, 1861, 1862, 1863, 1864, 1865, 1866, 1867, 1868, 1869, 1870, 1871, 1872, 1873, 1874, 1875, 1876, 1877, 1878, 1879, 1880, 1881, 1882, 1883, 1884, 1885, 1886, 1887, 1888, 1889, 1890, 1891, 1892, 1893, 1894, 1895, 1896, 1897, 1898, 1899, 1900)
        validYearList=range(1959,2005)
        # regrid for history data, only r1i1p1, if already done, not reprocess
        #monthlyRegrid(variable, dirHistorical, tmpdirHist, validYearList, 'historical_r1i1p1')

        # regrid for projections, only r1i1p1
        validYearList=range(2006, 2060)
        #monthlyRegrid(variable, indir, tmpdirNew, validYearList)
        
        # note: to correctly compute the ensemble means, consider deleting unwanted files!!!

        # for historical data set, selecting only r1i1p1:
        # monthlyAvg(variable, tmpdirHist, outdirHist, 1959, 2005, 'historical_r1i1p1')
        # for projections data set, selecting only r1i1p1:
        monthlyAvg(variable, tmpdirNew, outdir, 2006, 2059, 'r1i1p1')


