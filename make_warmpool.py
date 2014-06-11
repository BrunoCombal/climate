#!/usr/bin/env python
# to call from bash, use:
# source /usr/local/uvcdat/1.2.0/bin/setup_cdat.sh

## \author Bruno Combal
## \date June 2013
import cdms2
from cdms2 import MV
import numpy
import glob
import sys
import os
from os import path
import re
import string

# ___________________________
def usage():
    text='SYNOPSIS:\n\t{0} -indir sstdir -fileBasename basename -var varId -start startYear -end endYear [-bounds xmin xmax] -outdir outdir'.format(os.path.basename(__file__))
    return text
# ___________________________
def exitMessage(msg, exitCode='1'):
    print msg
    print
    print usage()
    sys.exit(exitCode)
#____________________________
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
    latAxis.units='degree_north'
    latAxis.long_name='Latitude'

    lonAxis = cdms2.createAxis(lon, lon_bnds)
    lonAxis.designateLongitude(True, 360.0)

    return(cdms2.createGenericGrid(latAxis, lonAxis, lat_bnds, lon_bnds))
# ___________________________
def do_convolve(data, thisFilter, lower, upper, xs, xe, ys, ye):

    filterDims=thisFilter.shape

    if len(filterDims) > 2:
        print 'wrong filter dimensions. Exit(200).'
        sys.exit(200)

    dataOut = data

    for il in range(ys, ye+1):
        for ic in range(xs, xe+1):
            total=0
            count=0
            for kl in range(filterDims[0]):
                for kc in range(filterDims[1]):
                    if (data[il][ic] > lower) and (data[il][ic]<upper):
                        total = total + thisFilter[kl][kc] * data[il][ic]
                        count = count + thisFilter[kl][kc]

            if count > 0:
                dataOut[il][ic] = total / float(count)

    return dataOut
# ___________________________
def do_interp(data, lower, upper, xs,xe,ys,ye):

    dataOut = data
    
    for il in range(ys, ye+1):
        if (data[il][xs-1] > lower) and (data[il][xs-1] < upper):
            v0 = data[il][xs-1]

            xfinal=xe+1 # where does the interpolation finish?
            if (data[il][xfinal] > lower) and (data[il][xfinal] < upper):
                v1 = data[il][xfinal]
            else:
                for xfinal in range(xe+1, xs+1, -1):
                    if (data[il][xfinal] > lower) and (data[il][xfinal] <upper):
                        break # we keep xfinal for the final position for interpolating

            v1 = data[il][xfinal]

            for ic in range(xs, xfinal):
                if dataOut[il][ic] > lower and dataOut[il][ic] < upper: # do not change no data
                    t = (ic-xs + 1)/(xe-xs+1)
                    dataOut[il][ic] = (1-t) * v0 + t * v1

    return dataOut
# ___________________________
# WP definition: all values must be >= threshold
def do_yearlyWPall(sstdir, sstrootname, variable, outdir, yearStart, yearEnd, latWindow=None):

    threshold = 28 + 273.15
    nodata = 1.e20
    varUnits=None
    referenceGrid=None
    
    referenceGrid = makeGrid()
    latws, lonwts = referenceGrid.getWeights()
    weights = MV.outerproduct(latws, lonwts)
    EarthSurface = 510072000

    areaWP=[]
    
    latWindowMatrix = None

    for iyear in range(yearStart, yearEnd+1):
        warmpool = None
        wnodata = None
        maxwarm = None
        minWarm = None
        print 'Processing year {0}'.format(iyear)
        for imonth in range(1, 12+1):
            idate = '{0}{1:02}'.format(iyear, imonth)
            fname = '{0}/{1}_{2}.nc'.format(sstdir, sstrootname, idate)
            thisFile = cdms2.open(fname)
            thisVar = numpy.ravel(thisFile[variable][:])

            if latWindow is not None:
                if latWindowMatrix is None:
                    latWindowMatrix = numpy.zeros(thisFile[variable].shape)
                    for ii in xrange(latWindow[0], latWindow[1]+1):
                        latWindowMatrix[:,ii] = 1
                thisVarTmp = thisVar
                thisVar = numpy.multiply( thisVarTmp, numpy.ravel(latWindowMatrix) )


            if warmpool is None:
                dimVar = numpy.squeeze(thisFile[variable][:]).shape
                varUnits = thisFile[variable].units
                warmpool = numpy.ravel(numpy.zeros(dimVar))
                wnodata  = (thisVar >= nodata)
                maxWarm = numpy.ravel(numpy.zeros(dimVar))
                minWarm = numpy.ravel(numpy.zeros(dimVar))
                monthMin =  numpy.ravel(numpy.zeros(dimVar))
                monthMax =  numpy.ravel(numpy.zeros(dimVar))
                # on first iteration, no comparison to the previous state (warmpool>=threshold)
                wwp = (thisVar >= threshold)  * (thisVar < nodata)
                wmax = (thisVar < nodata) * (thisVar > nodata) # set to false, everywhere
                wmin = (thisVar < nodata) * (thisVar > nodata) # set to False, everywhere

                maxWarm[wwp] = thisVar[wwp]
                minWarm[wwp] = thisVar[wwp]
                monthMax[wwp] = imonth
                monthMin[wwp] = imonth
            else:
                # warmpool: for all months, temperature > threshold
                # means that current value AND memo value are > threshold
                wwp = (thisVar >= threshold) * (warmpool>=threshold) * (thisVar < nodata)
                wmax = (thisVar >= threshold) * (warmpool>=threshold) * (thisVar < nodata) * (thisVar >= warmpool)
                wmin = (thisVar >= threshold) * (warmpool>=threshold) * (thisVar < nodata) * (thisVar < warmpool)

            # reset warmpool to 0, to keep only the minimal extension
            # we will encode the max and min observed 
            if wwp.any():
                maxWarm[:]=0
                minWarm[:]=0

                maxWarm[wwp] = warmpool[wwp]
                if wmax.any():
                    maxWarm[wmax] = thisVar[wmax]
                    monthMax[wmax]=imonth

                minWarm[wwp] = warmpool[wwp]
                if wmin.any():
                    minWarm[wmin] = thisVar[wmin]
                    monthMin[wmin]=imonth

                warmpool[:]=0 # reset warmpool to keep intersection between iterations
                warmpool[wwp]=thisVar[wwp]

            else:
                warmpool[:]=0
            thisFile.close()

        # ensure mask is set
        if wnodata.any():
            warmpool[wnodata] = nodata
            maxWarm[wnodata] = nodata
            minWarm[wnodata] = nodata
            monthMin[wnodata]=nodata
            monthMax[wnodata]=nodata
        
        wtonull = (warmpool ==0)
        if wtonull.any():
            monthMin[wtonull]=0
            monthMax[wtonull]=0

        warea = (warmpool >= threshold ) * (warmpool <nodata)
        surface = EarthSurface * numpy.ravel(weights)[warea].sum()
        areaWP.append([iyear, surface])

        wpOut = cdms2.createVariable(warmpool.reshape(dimVar), typecode='f', id='warmpool', \
                                         fill_value=nodata, grid=referenceGrid, copyaxes=0, \
                                         attributes=dict(long_name='warmpool, all temperatures method, year {0}'.format(iyear), units=varUnits))
        wpMax = cdms2.createVariable(maxWarm.reshape(dimVar), typecode='f', id='warmpool_max', \
                                         fill_value=nodata, grid=referenceGrid, copyaxes=0, \
                                         attributes=dict(long_name='warmpool max temperature, year {0}'.format(iyear), units=varUnits))
        wpMin = cdms2.createVariable(minWarm.reshape(dimVar), typecode='f', id='warmpool_min', \
                                         fill_value=nodata, grid=referenceGrid, copyaxes=0, \
                                         attributes=dict(long_name='warmpool min temperature, year {0}'.format(iyear), units=varUnits))
        monthMin = cdms2.createVariable(monthMin.reshape(dimVar), typecode='i', id='min_date', \
                                            fill_value=0, grid=referenceGrid, copyaxes=0, \
                                            attributes=dict(long_name='warmpool month of min(1-12), year {0}'.format(iyear), units=varUnits))
        monthMax = cdms2.createVariable(monthMax.reshape(dimVar), typecode='i', id='max_date', \
                                            fill_value=0, grid=referenceGrid, copyaxes=0, \
                                            attributes=dict(long_name='warmpool month of max(1-12), year {0}'.format(iyear), units=varUnits))
        
        outfilename='{0}/warmpool_{1}.nc'.format(outdir, iyear)
        if os.path.exists(outfilename): os.remove(outfilename)
        outfile=cdms2.open(outfilename, 'w')
        outfile.write(wpOut)
        outfile.write(wpMax)
        outfile.write(wpMin)
        outfile.write(monthMin)
        outfile.write(monthMax)
        outfile.close()

    return areaWP
# ____________________________
# WP definition: average annual temperature >= threshold
def do_yearlyWPAvg(sstdir, sstrootname, variable, outdir, yearStart, yearEnd, threshold=28.5+273.15, latWindow=None):
    nodata=1.e20
    varUnits=None
    thisFilter=numpy.ones((3,3))
    #thisFilter[1][1]=1
    
    EarthSurface = 510072000
    factor = 1 #(2*85*360) /( 360.0*180.0)

    areaWP=[]
    latWindowMatrix = None
    for iyear in range(yearStart, yearEnd+1):
        tempAvg=None
        wdata=None
        weights = None
        dimVar=None
        counter = None
        print 'Processing year {0}'.format(iyear) 
        for imonth in range(1, 12+1):
            idate = '{0}{1:02}'.format(iyear, imonth)
            fname = '{0}/{1}_{2}.nc'.format(sstdir, sstrootname, idate)
            thisFile = cdms2.open(fname)
            thisVar = numpy.ravel(thisFile[variable][:])
    
            if latWindow is not None:
                if latWindowMatrix is None:
                    latWindowMatrix = numpy.zeros(thisFile[variable].shape)
                    for ii in xrange(latWindow[0], latWindow[1]+1):
                        latWindowMatrix[:,ii] = 1
                thisVarTmp = thisVar
                thisVar = numpy.multiply( thisVarTmp, numpy.ravel(latWindowMatrix) )

            if tempAvg is None: # settings
                thisGrid = thisFile[variable].getGrid()
                if thisGrid is None:
                    thisGrid=makeGrid()
                (latws, lonws) = thisGrid.getWeights()
                weights = MV.outerproduct(latws, lonws)
                wdata = thisVar < nodata 
                tempAvg = numpy.zeros(thisVar.shape)
                tempAvg[wdata] = thisVar[wdata]

                dimVar = numpy.squeeze(thisFile[variable][:]).shape
                counter = numpy.zeros(tempAvg.shape, dtype='float')
                counter[wdata] = 1
            else:
                wdata = thisVar < nodata
                counter[wdata]=counter[wdata]+1
                tempAvg[wdata] = tempAvg[wdata] + thisVar[wdata]
        # compute average
        wdivide = counter > 0
        avg = numpy.zeros(tempAvg.shape)
        if wdivide.any:
            avg[wdivide] = tempAvg[wdivide] / counter[wdivide]
        # set to  areas < threshold
        wtzero = avg < threshold
        avg[wtzero] = 0
            
        # compute current area
        warea = (avg >= threshold) * (avg < nodata)
        area = factor * EarthSurface * numpy.ravel(weights)[warea].sum()
        areaWP.append([iyear, area])
        
        # create variables
        outAreaTmp = numpy.reshape(avg, dimVar)

        # filter before saving: grid stiching area: less data for ensemble mean here
        outAreaBis = do_interp(outAreaTmp, threshold, nodata, 156, 156+3, 0, 4*80) # rebuild
        outArea = do_convolve(outAreaBis, thisFilter,  threshold, nodata, 156+2 , 156+5, 0, 4*80) # smoothen

        wpOut = cdms2.createVariable(outArea, typecode='f', id='warmpool', \
                                         grid=thisGrid, copyaxes=1, \
                                         attributes=dict(long_name='warmpool, average temperature method, year {0}'.format(iyear), units=varUnits))
        # write to file
        outfilename='{0}/warmpool_{1}.nc'.format(outdir, iyear)
        if os.path.exists(outfilename): os.remove(outfilename)
        outfile=cdms2.open(outfilename, 'w')
        outfile.write(wpOut)
        outfile.close()
        # close files
        thisFile.close()

    return areaWP
# ____________________________
if __name__=="__main__":

    sstdir=None
    outdir=None
    startYear=None #2010
    endYear=None # 2059
    bounds=None #[70 580]
    var=None #'mean_mean_tos'
    fileBasename=None #'ensemble_tos_rcp85'

    ii = 1
    while ii < len(sys.argv):
        arg = sys.argv[ii].lower()
        
        if arg == '-indir':
            ii = ii + 1
            sstdir=sys.argv[ii]
        elif arg == '-filebasename':
            ii = ii + 1
            fileBasename=sys.argv[ii]
        elif arg == '-start':
            ii = ii + 1
            startYear = int(sys.argv[ii])
        elif arg == '-end':
            ii = ii + 1
            endYear = int(sys.argv[ii])
        elif arg == '-bounds':
            ii = ii + 1
            bounds=[]
            bounds.append(int(sys.argv[ii]))
            ii = ii + 1
            bounds.append(int(sys.argv[ii]))
        elif arg=='-var':
            ii = ii + 1
            var = sys.argv[ii]
        elif arg == '-outdir':
            ii = ii + 1
            outdir=sys.argv[ii]

        ii = ii + 1

    # check parameters exist
    if sstdir is None:
        exitMessage('Missing an input directory, use option -indir. Exit(1).',1)
    if startYear is None:
        exitMessage('Missing a starting year, use option -start. Exit(2).',2)
    if endYear is None:
        exitMessage('Missing an ending year, use option -end. Exit(3).', 3)
    if outdir is None:
        exitMessage('Missing an output directory, use option -outdir. Exit(4).',4)
    if var is None:
        exitMessage('Missing a variable identifier for the netcdf files. Exit(5).', 5)
    if fileBasename is None:
        exitMessage('Missing a filebasename, use option -fileBasename. Exit(6).', 6)

#    rcp='8'
##    sstdir='/data/tmp/new_algo/tos_rcp85'
#    sstdir='/data/tmp/new_algo/tos_rcp{0}5'.format(rcp)
##    sstdirHist='/data/cmip5/rcp/rcp8.5/toshist_ensemble'
#    outdir='/data/cmip5/rcp/rcp{0}.5/tos_warmpools'.format(rcp)
    
    # all temp, projections
    areaWP = do_yearlyWPAvg(sstdir, fileBasename, var, outdir, startYear, endYear, threshold=28.5+273.15, latWindow=bounds)
    # avg temp, projections
    #areaWP = do_yearlyWPAvg(sstdir, 'modelmean_tos', 'tos', outdir, 2006, 2059, 28+273.15)
    # avg temp, hist
    #areaWP = do_yearlyWPAvg(sstdirHist, 'modelmean_tos', 'tos', outdir, 1850, 2005, 28+273.15)
    for ii, area in areaWP:
        print '{0},{1},{2}'.format(ii, area, area / areaWP[0][1])

# end of script
