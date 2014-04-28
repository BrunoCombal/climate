#!/usr/bin/env python

# to run the script with the correct version of uvcdat:
#  source /usr/local/uvcdat/1.2.0/bin/setup_cdat.sh

# series 1 (yyyy, mm): original data
# series 2 (mm): data to select
# select data in series 2, for the month when max was observed in series 1

import cdms2
from cdms2 import MV
import numpy
import glob
import sys
import os
from os import path
import re
from scipy import interpolate


## \brief message displayed on exit
def messageOnExit(message=None, exitCode=1):
    if message is not None:
        print message
        print

    print "Compute the rms of a set of files"
    print 'Usage: nc_rms.py [-p datapath] -v var -o outfile file*'

    sys.exit(exitCode)


# _______________
if __name__=="__main__":

    indir='/data/tmp/oimonth_v2/'
    
    fileList=glob.glob(indir+'oiv2mon.*.nc')

    maxVal = None
    maxMonth = None
    dataShape=None
    for ifile in fileList:
        fname=os.path.basename(ifile)
        # get month from the filename
        monthStr=fname[12:14]
        
        fh=cdms2.open(ifile,'r')
        data = numpy.ravel(fh['sst'][:])
        if maxVal is None:
            maxVal = numpy.ravel(data)
            maxMonth = numpy.zeros(maxVal.shape) + int(monthStr)
            dataShape = fh['sst'][:].shape
        else:
            wtm = data >= maxVal
            if wtm.any():
                maxVal[wtm]=data[wtm]
                maxMonth[wtm]=int(monthStr)

        fh.close()

    # resize the max date array
    # dateResized = maxMonth

    # now we now when occurs the max, let's select the corresponding rms
    print 'distributing dates'
    rms_at_maxsst=numpy.zeros(maxMonth.ravel().shape)
    for imonth in range(1,13):
        print '/data/sst/oimonth_v2/rms_sst_1981_2000_{0:02}.nc'.format(imonth)
        fh=cdms2.open('/data/sst/oimonth_v2/sst_rms_{0:02}.nc'.format(imonth),'r')
        rms=numpy.ravel(numpy.flipud(fh['rms_sst'])[:])
        wtm = maxMonth.ravel() == imonth
        if wtm.any():
            print rms_at_maxsst.shape, rms.shape
            rms_at_maxsst[wtm] = rms[wtm]
        fh.close()

    print 'saving data'
    outname='/data/sst/oimonth_v2/rms_at_maxsst.nc'
    if os.path.exists(outname): os.remove(outname)
    outfile=cdms2.open(outname,'w')
    outData=cdms2.createVariable(rms_at_maxsst.reshape(dataShape), typecode='f', id='rms_at_max')
    outfile.write(outData)
    outfile.close()

    outname='/data/tmp/max.nc'
    if os.path.exists(outname): os.remove(outname)
    outfile=cdms2.open(outname,'w')
    outMax=cdms2.createVariable(maxVal.reshape(dataShape), typecode='f', id='max')
    outMaxDate=cdms2.createVariable(maxMonth.reshape(dataShape), typecode='i',id='month')
    outfile.write(outMax)
    outfile.write(outMaxDate)
    outfile.close()
