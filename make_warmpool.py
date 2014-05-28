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
# WP definition: all values must be >= threshold
def do_yearlyWPall(sstdir, sstrootname, variable, outdir, yearStart, yearEnd):
    
    threshold = 28.75 + 273.15
    nodata = 1.e20
    varUnits=None
    referenceGrid=None
    
    referenceGrid = makeGrid()
    latws, lonwts = referenceGrid.getWeights()
    weights = MV.outerproduct(latws, lonwts)
    EarthSurface = 510072000

    areaWP=[]
    
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
def do_yearlyWPAvg(sstdir, sstrootname, variable, outdir, yearStart, yearEnd, threshold=28.5+273.15):
    nodata=1.e20
    varUnits=None
    
    EarthSurface = 510072000
    factor = (2*85*360) /( 360.0*180.0)
    
    areaWP=[]
    print 'variable ',variable
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
            print 'opening file {0}'.format(fname)
            thisFile = cdms2.open(fname)
            thisVar = numpy.ravel(thisFile[variable][:])

            if tempAvg is None: # settings
                thisGrid = thisFile[variable].getGrid()
                if thisGrid is None:
                    thisGrid=makeGrid()

                (latws, lonws) = thisGrid.getWeights()
                weights = MV.outerproduct(latws, lonws)
                tempAvg = thisVar
                wdata = thisVar < nodata # the data creation has ensured that nodata are aligned accross the volume
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
        outArea = numpy.reshape(avg, dimVar)

        wpOut = cdms2.createVariable(outArea, typecode='f', id='warmpool', \
                                         fill_value=nodata, grid=thisGrid, copyaxes=1, \
                                         attributes=dict(long_name='warmpool, average temperature method, year {0}'.format(iyear), units=varUnits))
        # write to file
        outfilename='{0}/warmpool_{1}.nc'.format(outdir, iyear)
        if os.path.exists(outfilename): os.remove(outfilename)
        outfile=cdms2.open(outfilename, 'w')
#        outfile.write(wpOut)
        outfile.close()
        # close files
        thisFile.close()
    return areaWP
# ____________________________
if __name__=="__main__":

    sstdir='/data/tmp/new_algo/tos_rcp85'
    sstdirHist='/data/cmip5/rcp/rcp8.5/toshist_ensemble'
    outdir='/data/cmip5/rcp/rcp8.5/tos_warmpools'
    
    # all temp, projections
    areaWP = do_yearlyWPAvg(sstdir, 'ensemble_tos_rcp85', 'mean_mean_tos', outdir, 2010, 2059)
    # avg temp, projections
    #areaWP = do_yearlyWPAvg(sstdir, 'modelmean_tos', 'tos', outdir, 2006, 2059, 28+273.15)
    # avg temp, hist
    #areaWP = do_yearlyWPAvg(sstdirHist, 'modelmean_tos', 'tos', outdir, 1850, 2005, 28+273.15)
    for ii, area in areaWP:
        print '{0}\t{1}'.format(ii, area)


# end of script
