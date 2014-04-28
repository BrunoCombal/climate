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
from scipy import interpolate
import shutil

# this code to resample and resize real climato data
# input: 180*360, output from 85S :85N, 340*720



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
# ________________
def saveData(outfilename, data, typecode, id, fill_value, grid, copyaxes, attribute1, attribute2, latAxis, lonAxis):
    
    # for netcdf3: set flags to 0
    cdms2.setNetcdfShuffleFlag(1)
    cdms2.setNetcdfDeflateFlag(1)
    cdms2.setNetcdfDeflateLevelFlag(3)

    if os.path.exists(outfilename): os.remove(outfilename)

    outfile = cdms2.open( outfilename, 'w')
    var = cdms2.createVariable(data, typecode=typecode, id=id, fill_value=fill_value, grid=grid, copyaxes=copyaxes, attributes=dict(long_name=attribute1, units=attribute2) )
    var.setAxisList((latAxis, lonAxis))
    outfile.write(var)
    outfile.close()
# ______________________
def do_resize(var, fname, limit=100, interpol='nearest'):
    infile = cdms2.open(fname,'r')
    yvar=numpy.squeeze(infile[var][:]) # y, x
    dim=yvar.shape

    wnodata = yvar >= limit
    print 'wnodata',wnodata.shape
    if wnodata.any():
        print 'correcting'
        yvar[wnodata] = 1.e20

    print dim
    # les points sont numerotes de 0 a 179 (-90 a 90 degrees)
    # on ne veut que les points de 4 a 174 (i.e. de -85 a 85), par pas de 2

    points=numpy.zeros( (dim[0], dim[1], 2))
    for jj in range(dim[0]):
        for ii in range(dim[1]):
            points[jj, ii, 0] = jj
            points[jj, ii, 1] = ii

    outGY, outGX = numpy.mgrid[ 4 : 174 :0.5 , 0:360:0.5 ]

    outvar = interpolate.griddata(numpy.reshape(points, (dim[0]*dim[1],2)), numpy.ravel(yvar), (outGY, outGX), method=interpol)

    infile.close()

    return numpy.flipud(outvar)
# ______________________
# a general purpose function
def do_resize_multi(var, fname, limit=100, interpol='nearest'):
    infile = cdms2.open(fname,'r')
    yvar=numpy.squeeze(infile[var][:]) # y, x
    dim=yvar.shape

    wnodata = numpy.ravel(yvar) >= limit
    print 'wnodata multi',wnodata.shape,limit
    if wnodata.any():
        print 'correcting multi'
        tmp=numpy.ravel(yvar)
        tmp[wnodata]=1.e20
        print 'dim tmp',tmp.shape
        yvar=numpy.reshape(tmp, dim)

    print dim
    # les points sont numerotes de 0 a 179 (-90 a 90 degrees)
    # on ne veut que les points de 4 a 174 (i.e. de -85 a 85), par pas de 2

    points=numpy.zeros( (dim[1], dim[2], 2))
    for jj in range(dim[1]):
        for ii in range(dim[2]):
            points[jj, ii, 0] = jj
            points[jj, ii, 1] = ii

    outGY, outGX = numpy.mgrid[ 4 : 174 :0.5 , 0:360:0.5 ]
    outdim=outGY.shape

    outvar=numpy.zeros( (dim[0], outdim[0], outdim[1] )) + 1.e20
    for itime in range(dim[0]):
        tmp = interpolate.griddata(numpy.reshape(points, (dim[1]*dim[2],2)), numpy.ravel(yvar[itime]), (outGY, outGX), method=interpol)
        outvar[itime] = numpy.flipud(tmp)

    infile.close()

    return outvar
# _____________________
# replicate data by an entiere number of times
def do_resize_int(var, fname, nodata, limit=100):
    
    infile = cdms2.open(fname, 'r')
    yvar = numpy.squeeze(infile[var][:])
    dims = yvar.shape

    # change no data
    yvar = numpy.ravel(yvar)
    wnodata = yvar > 100
    if wnodata.any():
        yvar[wnodata] = nodata
    yvar = yvar.reshape(dims)

    if len(dims)==3:
        print 'with time'
        outmatrix = numpy.zeros( (dims[0], dims[1]*2, dims[2]*2)) + nodata
        for jjSrc in xrange(dims[1]):
            for iiSrc in xrange(dims[2]):
                thisData = yvar[ : , jjSrc, iiSrc]
                outmatrix[ : , jjSrc*2  , iiSrc*2 ] = thisData 
                outmatrix[ : , jjSrc*2+1, iiSrc*2 ] = thisData
                outmatrix[ : , jjSrc*2  , iiSrc*2+1 ] = thisData
                outmatrix[ : , jjSrc*2+1, iiSrc*2+1 ] = thisData

        tmpShape = outmatrix.shape
        outData = numpy.zeros( (tmpShape[0], tmpShape[1]-20, tmpShape[2]) )
        for itime in xrange(outData.shape[0]):
            for jjSrc in xrange(outData.shape[1]):
                thisData = outmatrix[itime, jjSrc + 10, :]
                outData[itime, jjSrc, :] = thisData

        # flipud
        returnData = numpy.zeros(outData.shape)
        for itime in xrange(dims[0]):
            returnData[itime]=numpy.flipud(outData[itime,:,:])
    else:
        print 'no time dimension'
        outmatrix = numpy.zeros( (dims[0]*2, dims[1]*2) ) + nodata
        for jjSrc in xrange(dims[0]):
            for iiSrc in xrange(dims[1]):
                thisData = yvar[jjSrc, iiSrc]
                outmatrix[ jjSrc*2  , iiSrc*2  ] = thisData 
                outmatrix[ jjSrc*2+1, iiSrc*2  ] = thisData
                outmatrix[ jjSrc*2  , iiSrc*2+1] = thisData
                outmatrix[ jjSrc*2+1, iiSrc*2+1] = thisData

        tmpShape = outmatrix.shape
        outData = numpy.zeros( (tmpShape[1]-20, tmpShape[2]) )
        for jjSrc in xrange(outdata.shape[0]):
            outData[jjSrc, :] = outmatrix[jjSrc+10, :]

        returnData = numpy.flipud(outData)

    infile.close()

    return returnData

# ______________________
if __name__=='__main__':
    nodata = 1.e20
    (referenceGrid, latAxis, lonAxis, latBounds, lonBounds) = makeGrid()
    # for netcdf3: set flags to 0
    cdms2.setNetcdfShuffleFlag(1)
    cdms2.setNetcdfDeflateFlag(1)
    cdms2.setNetcdfDeflateLevelFlag(3)

    indir='/data/sst/reynolds_climatology/noaa_oist_v2/'
    outdir='/data/tmp/new_algo/resizing/'#'/data/sst/reynolds_climatology/noaa_oist_v2/resized_framed/'
    
    infile=indir+'/max_sst.ltm.1971-2000.nc'
    infile=indir+'sst.ltm.1971-2000.nc'
    outvar = do_resize_int('sst', infile, nodata, 100)
    print outvar.shape
    print lonAxis.shape
    print latAxis.shape
    

    outfile=outdir+'new_resized.nc'
    if os.path.exists(outfile): os.remove(outfile)
    fh=cdms2.open(outfile, 'w')
    var = cdms2.createVariable(outvar, typecode='f', id='sst', fill_value=nodata, grid=referenceGrid, copyaxes=1 )
    fh.write(var)
    fh.close()
    sys.exit(1)
    


    #saveData(outdir+'/max_sst.ltm.1971-2000_resized.nc', outvar, typecode='f', id='sst', fill_value=nodata, grid=referenceGrid, copyaxes=1, attribute1='real Climato max',attribute2='Degrees Celsius',latAxis=latAxis,lonAxis=lonAxis)
    
    #infile=indir+'sst.ltm.1971-2000.nc'
    #outvar = do_resize_multi('sst',infile)
    #outfile=outdir+'sst.ltm.1971-2000_resized.nc'
    #if os.path.exists(outfile): os.remove(outfile)
    #fh=cdms2.open(outfile, 'w')
    #var = cdms2.createVariable(outvar, typecode='f', id='sst', fill_value=nodata, grid=referenceGrid, copyaxes=1 )
    #fh.write(var)
    #fh.close()

    # resize rms at max sst
    #infile='/data/sst/oimonth_v2/rms_at_maxsst.nc'
    #print 'processing ',infile
    #outvar = do_resize('rms_at_max',infile)
    #saveData(outdir+'/rms_at_maxsst_resized.nc', outvar, typecode='f', id='rms_at_max', fill_value=nodata,grid=referenceGrid, copyaxes=1, attribute1='rms at maximum sst',attribute2='Degrees Celsius',latAxis=latAxis, lonAxis=lonAxis)
