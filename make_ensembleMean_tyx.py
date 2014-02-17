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
    textUsage='SYNOPSIS:\n\tmake_ensemble_Mean_tyx.py -v VARIABLE -path PATHIN -outdir PATHOUT [-tmpdir TMPPATH] \n\t-minVar MINVAL -maxVar MAXVAL\tn-model MODELLIST -startYear STARTYEAR -endYear ENDYEAR [-monthList MONTHLIST]\n\t[-regridFirst REGRIDBOOL] [-deleteGrid DELETEBOOL] -rcp RCP\n'
    textUsage=textUsage+'\tVARIABLE: a netcdf CMIP5 variable name, such as tos, zos, so, thetao;\n'
    textUsage=textUsage+'\tPATHIN: input data directory (does not support sub-directories);\n'
    textUsage=textUsage+'\tPATHOUT: output directory, created if does not exist;\n'
    textUsage=textUsage+'\tTMPPATH: temporary path. Default: a random pathname is defined at runtime, as a leaf of PATHOUT;\n'
    textUsage=textUsage+'\tMINVAL: any value below minVar is considered as nodata;\n'
    textUsage=textUsage+'\tMAXVAL: any value above maxVar is considered as nodata;\n'
    textUsage=textUsage+'\tMODELLIST: a text file with a model name per name, the model name is used to select the files to process;\n'
    textUsage=textUsage+'\tSTARTYEAR: first year in the series of dates to process;\n'
    textUsage=textUsage+'\tENDYEAR: last year in the series of date to process;\n'
    textUsage=textUsage+'\tMONTHLIST: a comma separated list of month, such as "1,2,3" or "1,6,12". Values range is [1, 12].\n'
    textUsage=textUsage+'In first place, the programme will average model output per model (if a model output has several rXiYpZ ensemble, they are averaged. Then, the averages are averaged to produce the ensemble mean;\n'
    textUsage=textUsage+'\tREGRIDBOOL\n'
    textUsage=textUsage+'\tDELETEBOOL\n'
    textUsage=textUsage+'\tRCP a string corresponding to the RCP string to match in filenames.\n'
    textUsage=textUsage+'Averages are computed for each month of the year.\n'
    return textUsage
# ____________________________
def exitMessage(msg, exitCode='1'):
    thisLogger.critical(msg)
    print msg
    print
    print usage()
    sys.exit(exitCode)
# ___________________________
def boolConvert(code):
    if code=='0':
        return False
    if code.lower()=='false':
        return False
    if code.lower()=='no':
        return False
    if code=='1':
        return True
    if code.lower()=='true':
        return True
    if code.lower()=='yes':
        return True
# ____________________________
def decodeMonthList(parameter):

    listMonth = [int(x) for x in parameter.strip().split(',')]
    for ii in listMonth:
        if ii<1 or ii>12:
            exitMessage('month defined in the month list must be in [1, 12]. Exit(100).',100)
    return listMonth
# ____________________________
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for x in range(size))
#_____________________________
def flatten(foo):
    for x in foo:
        if hasattr(x, '__iter__'):
            for y in flatten(x):
                yield y
        else:
            yield x
# ____________________________
# dict{date:[filename]}
def agregateDict(refDict, newDict):

    if refDict is None and newDict is None:
        return None    
    # get list of all keys
    if refDict is None:
        return newDict
    if len(refDict)==0:
        return newDict
    if newDict is None:
        return refDict
    if len(newDict)==0:
        return refDict
    
    keyList = sorted(set(refDict.keys() + newDict.keys()))

    result={}
    for ikey in keyList:
        val = []
        if ikey in refDict.keys(): val.append( refDict[ikey] )
        if ikey in newDict.keys(): val.append( newDict[ikey] )
        result[ikey] = [ x for x in flatten(val) ]

    del val
    gc.collect()
    return result
# ____________________________
def makeGrid(thisStep=0.5):
    xstart=0
    xend=360
    xstep=thisStep
    ystart=-85
    yend=85
    ystep=thisStep

    lon_bnds=[]
    lon=[]
    for ii in numpy.arange(xstart, xend, xstep):
        lon_bnds.append( [ii, ii + xstep] )
        lon.append(ii+0.5*xstep)
    lon_bnds=numpy.array(lon_bnds)
    lon=numpy.array(lon)

    lat_bnds=[]
    lat=[]
    for ii in numpy.arange(ystart, yend, ystep):
        lat_bnds.append([ii, ii + ystep])
        lat.append(ii+0.5*ystep)
    lat_bnds=numpy.array(lat_bnds)
    lat=numpy.array(lat)

    latAxis = cdms2.createAxis(lat, lat_bnds)
    latAxis.designateLatitude(True)
    latAxis.units='degrees_north'
    latAxis.id='latitude'
    latAxis.long_name='Latitude'

    lonAxis = cdms2.createAxis(lon, lon_bnds)
    lonAxis.designateLongitude(True, xend)
    lonAxis.designateCircular(xend)
    lonAxis.units='degrees_east'
    lonAxis.id='longitude'
    lonAxis.long_name='Longitude'

    return((cdms2.createGenericGrid(latAxis, lonAxis, lat_bnds, lon_bnds), latAxis, lonAxis, lat_bnds, lon_bnds))
# ____________________________
def do_cleanNodataLines(var, nodata):

    oneSlice = numpy.squeeze(var[:,:,0])

    refShape=oneSlice.shape
    # where are the nodata vertical lines?
    # 1./ transform the slice: 0=data, 1=nodata
    test = numpy.zeros(oneSlice.shape)
    wto1 = oneSlice >= nodata
    if wto1.any():
        test[wto1] = 1
    else:
        thisLogger.info('do_cleanNodataLines: no-data is missing from this dataset. Return.')
        return var

    # 2./ multiplications: if there are only nodata, results is 1
    line = numpy.array(oneSlice[0, :]) # copy first line
    for il in range(oneSlice.shape[1]):
        line = line * oneSlice[il, :]

    # 3./ do we have a 1 somewhere? It means that there was only nodata along the line
    wone = line == 1
    if wone.any():
        thisLogger.info('do_cleanNodataLines: found {0} lines to correct.'.format(len(wone)))
    else:
        thisLogger.info('do_cleanNodataLines: found no line to correct.')
        return var
# ____________________________
# auto mask based on the principle that the mask does not change in-between dates
def autoMask(var, nodata):

    refshape = var.shape
    if len(refshape)==3:
        tmp = numpy.reshape(var, (refshape[0], refshape[1] * refshape[2]) )
    elif len(refshape)==4:
        tmp = numpy.reshape(var, (refshape[0], refshape[1] * refshape[2] * refshape[3]) )
    wtnodata = (tmp.max(axis=0) - tmp.min(axis=0)) < 0.001

    if wtnodata.any():
        for ii in range(refshape[0]):
            tmp[ii, wtnodata] = nodata
        var[:] = numpy.reshape(tmp, refshape)
    
    del tmp, wtnodata
    gc.collect()
    return var
# ____________________________
def updateCounters(accum, N, mini, maxi, data, minVar, maxVar, nodata=1.e20):

    if data is None:
        return [accum, N, mini, maxi]

    dim = numpy.squeeze(data[:]).shape

    if accum is None:
        accum = numpy.zeros(dim) + nodata
        N = numpy.zeros(dim) + nodata
        mini = data.copy()
        maxi = data.copy()

    wtadd = (data >= minVar ) * (data < maxVar) * (accum < nodata) # add where not nodata
    wtreplace = (data >= minVar) * (data < maxVar) * (accum >= nodata) # replace if no data
    wmax = (data >= maxi) * (data < nodata) * (data >= minVar) * (data < maxVar)
    wmaxReplace = (mini >= nodata) * (data < nodata) * (data >= minVar)
    wmin = (data <= mini) * (data >= minVar) * ( data < maxVar) * ( maxi < nodata )
    wminReplace = (mini >= nodata) * (data < nodata) * (data >= minVar)

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
        
    del wtadd, wtreplace, wmax, wmaxReplace, wmin, wminReplace
    gc.collect()
    return [accum, N, mini, maxi]
# ___________________________
def do_regrid(variable, lstInFile, outdir, stringBefore, yearStart, yearEnd, topLevel=0, bottomLevel=1000):

    createdFiles=[]
    nodata=1.e20
    
    if lstInFile is None:
        thisLogger.info( 'No file to process. Return' )
        return None

    if len(lstInFile)==0:
        thisLogger.info('Found no file to process, consider revising search pattern. Return.')
        return None

    (newGrid, latAxis, lonAxis, lat_bnds, lon_bnds) = makeGrid()

    for fileName in lstInFile:
        thisLogger.info('Regriding file: {0}'.format(fileName))

        thisFile = cdms2.open(fileName)
        # to reduce output file size and memory use, collect start/end times according to internal file encoding
        startTime = [t for t in thisFile[variable].getTime().asComponentTime() if (t.year==startYear)]
        endTime = [t for t in thisFile[variable].getTime().asComponentTime() if (t.year==endYear)]
        if len(startTime)==0 and len(endTime)==0: # this file does not contain useful data, next iteration
            continue
        if len(startTime)==0: # the first date is not in this file, process from the start
            startTime = thisFile[variable].getTime().asComponentTime()
        if len(endTime)==0: # the last date is not in this file, process up to the end
            endTime = thisFile[variable].getTime().asComponentTime()

        thisLogger.info('start time = {0}-{1:02}'.format(startTime[0].year, startTime[0].month) )
        thisLogger.info('end time = {0}-{1:02}'.format(endTime[-1].year, endTime[-1].month))

        if thisFile[variable].getLevel() is None:
            # some files do not have nodata set to 1.e20 (EC-EARTH), some have masked values set to something else (0 and 1.e20, for MRI): let's process our mask by identifying unchanged values
#            data = cdms2.createVariable(thisFile[variable].subRegion(time=(startTime[0], endTime[-1], 'cc')))
#            if (thisFile[variable].getMissing() < 1.e20): # some files do not have a mask set to 1.e20 or no mask at all (eg. EC-Earth and MRI with 0 in lands instead of 1.e20 
            tmp = cdms2.createVariable(thisFile[variable].subRegion( time=(startTime[0], endTime[-1], 'cc'), level=(topLevel, bottomLevel,'cc') ))
            data = autoMask(tmp, nodata)
            del tmp
            gc.collect()
#            else:
#                data = cdms2.createVariable(thisFile[variable].subRegion( time=(startTime[0], endTime[-1], 'cc'), level=(topLevel, bottomLevel,'cc') )) 
        else:
            verticalGrid = make_levels()
            topLevel = levelAxis.min()
            bottomLevel = levelAxis.max()
            if thisFile[variable].getMissing() is None:
                tmp = cdms2.createVariable(thisFile[variable].subRegion( time=(startTime[0], endTime[-1], 'cc'), level=(topLevel, bottomLevel,'cc') ))
                data = autoMask(tmp, nodata)
                del tmp
                gc.collect()
            else:
                data = cdms2.createVariable(thisFile[variable].subRegion( time=(startTime[0], endTime[-1], 'cc'), level=(topLevel, bottomLevel,'cc') )) 

        mask = numpy.array(data) < nodata
        if thisFile[variable].getLevel() is None:
            regrided = data.regrid(newGrid, missing=nodata, order=thisFile[variable].getOrder(), mask=mask)
        else:
            tmp = data.regrid(newGrid, missing=nodata, order=thisFile[variable].getOrder(), mask=mask)
            regrided = tmp.pressureRegrid( verticalGrid, method='linear')

        regrided.id=variable

        outfilename = '{0}/{1}{2}'.format(outdir, stringBefore, os.path.basename(fileName))
        createdFiles.append(outfilename )
        if os.path.exists(outfilename): os.remove(outfilename)
        outfile = cdms2.open(outfilename, 'w')
        outfile.write(regrided)
        outfile.close()
        thisFile.close()

        del mask, regrided
        gc.collect()

    del newGrid, latAxis, lonAxis, lat_bnds, lon_bnds
    gc.collect()
    return createdFiles
# ___________________________
# for a list of files: open all files, go from date 1 to date 2, compute avg for thisdate, save thisdate
# if a new grid is passed: regrid
def do_stats(variable, validYearList, monthList, lstInFile, outdir, stringBefore, outnameBase, minVar=-1.e20, maxVar=1.e20, doSTD=False):
    
    if validYearList is None:
        exitMessage('List of years to process is undefined, edit code. Exit 5.',5)

    createdFiles={}   
    nodata=1.e20

    if lstInFile is None:
        thisLogger.info('No file to process. Return.')
        return

    if len(lstInFile)==0:
        thisLogger.info('Found no file to process, consider revising search pattern.')
        return

    # open all files
    listFID=[]
    thisLogger.debug('Averaging with files')
    for ifile in lstInFile:
        thisLogger.debug(ifile)
        listFID.append(cdms2.open(ifile, 'r'))
        
    
    # go through the list of dates, compute ensemble average
    for iyear in validYearList:
        thisLogger.info('Processing year {0}'.format(iyear))
        for imonth in monthList:
            accumVar=None
            accumN=None
            mini=None
            maxi=None
            refGrid=None
            dims=None
            units=None
            for ifile in listFID:
                # get the actual file-date for this year/month (because day numbering varies, sometimes 1, 15 or 16...)
                if ifile[variable].getTime() is None:
                    if refGrid is None:
                        refGrid = ifile[variable].getGrid()
                        # axis=ifile[variable].getAxisList(omit='time')
                        dims = numpy.squeeze(ifile[variable]).shape
                        units= ifile[variable].units
                    [accumVar, accumN, mini, maxi]= updateCounters(accumVar, accumN, mini, maxi,
                                                                   numpy.array( ifile[variable]).ravel(),
                                                                   minVar, maxVar, nodata )
                else:
                    thisTime = [ii for ii in ifile[variable].getTime().asComponentTime() if (ii.year==iyear and ii.month==imonth)] 
                    if len(thisTime)==1:
                        if refGrid is None:
                            refGrid = ifile[variable].getGrid()
                            dims = numpy.squeeze(ifile[variable].subRegion(time=thisTime[0])).shape
                            units= ifile[variable].units

                        [accumVar, accumN, mini, maxi]= updateCounters(accumVar, accumN, mini, maxi,
                                                                       numpy.array( ifile[variable].subRegion(time=thisTime[0])).ravel(),
                                                                       minVar, maxVar, nodata )

            # compute average
            # it can happen that there is no data to process: if the input files for the current model has an ending date before the current date
            # in this case, accumN is None: do not save stats, and do not add a file name in createdFiles
            if accumN is not None:
                wtdivide = (accumN < nodata) * (accumN > 0)
            
            
                if wtdivide.any():
                    accumVar[wtdivide] = accumVar[wtdivide] / accumN[wtdivide]

                # compute std
                if doSTD:
                    thisLogger.info('Computing std: to be implemented')
                # create and save variables
                meanVar = cdms2.createVariable( accumVar.reshape(dims), typecode='f', id='mean_{0}'.format(variable),  fill_value=nodata, attributes=dict(long_name='mean', units=units) )
                meanVar.setGrid(refGrid)

                counter = cdms2.createVariable(accumN.reshape(dims), typecode='i', id='count', fill_value=nodata, attributes=dict(long_name='count', units='None') )
                counter.setGrid(refGrid)
                miniVar = cdms2.createVariable(mini.reshape(dims), typecode='f', id='minimum', fill_value=nodata, attributes=dict(long_name='minimum', units=units) )
                miniVar.setGrid(refGrid)
                maxiVar = cdms2.createVariable(maxi.reshape(dims), typecode='f', id='maximum', fill_value=nodata, attributes=dict(long_name='maximum', units=units) )
                maxiVar.setGrid(refGrid)

                outfilename = '{0}/{1}_{2}_{3}{4:02}.nc'.format(outdir, stringBefore, outnameBase, iyear, imonth )
                if os.path.exists(outfilename): os.remove(outfilename)
                thisLogger.debug('Saving stats to file {0}'.format(outfilename))
                outfile = cdms2.open(outfilename, 'w')
                outfile.write(meanVar)
                outfile.write(counter)
                outfile.write(miniVar)
                outfile.write(maxiVar)
                outfile.close()

                createdFiles['{0}{1:02}'.format(iyear,imonth)] = outfilename

    # close input files
    for ii in listFID: ii.close()

    del wtdivide, accumVar, mini, maxi, accumN
    gc.collect()
    return(createdFiles)
#___________________________
if __name__=="__main__":

    variable = None
    indir = None
    tmpdir = None
    outdir = None
    modelListFile=None
    startYear=None
    endYear=None
    monthList=range(1,13)
    regridFirst = True
    deleteRegrid = False
    modelStat = True
    rcp=None
    logFile='{0}.log'.format(__file__)
    minVar=-1.e20
    maxVar=1.e20
    topLevel=0
    bottomLevel=300

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
        elif arg=='-monthlist':
            ii = ii + 1
            monthList=decodeMonthList(sys.argv[ii])
        elif arg=='-regridfirst':
            ii=ii+1
            regridFirst=boolConvert(sys.argv[ii])
        elif arg=='-deleteregrid':
            ii = ii + 1
            deleteRegrid = boolConvert(sys.argv[ii])
        elif arg=='-rcp':
            ii=ii+1
            rcp=sys.argv[ii]
        elif arg=='-log':
            ii = ii + 1
            logFile = sys.argv[ii]
        ii = ii + 1

    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    thisLogger = logging.getLogger('MyLogger')
    thisLogger.setLevel(logging.DEBUG)
    handler = logging.handlers.RotatingFileHandler(logFile, maxBytes=1024*500, backupCount=5)
    thisLogger.addHandler(handler)

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
    if rcp is None:
        exitMessage('Please define an rcp, use option -rcp. Exit(15).',15)

    if tmpdir is None:
        tmpdir = '{0}/tmp_{1}'.format(outdir, id_generator() )

    if not os.path.exists(outdir): os.makedirs(outdir)
    if not os.path.exists(tmpdir): os.makedirs(tmpdir)

    # for netcdf3: set flag to 0
    cdms2.setNetcdfShuffleFlag(1)
    cdms2.setNetcdfDeflateFlag(1)
    cdms2.setNetcdfDeflateLevelFlag(3)

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

    processedFiles=None

    # process models individually
    for thisModel in modelList:
        thisLogger.info('Model {0}'.format(thisModel))
        pattern=re.compile('{0}_{1}_{2}_{3}_{4}_{5}.nc'.format(variable, 'Omon', thisModel, rcp, 'r.*i.*p.*', '.*') )
        lstInFile=[f for f in glob.glob('{0}/*.nc'.format(indir)) if (os.stat(f).st_size and pattern.match(os.path.basename(f) ) ) ]

        if regridFirst:
            regridedFiles = do_regrid(variable, lstInFile, tmpdir, 'regrid_', startYear, endYear, topLevel, bottomLevel)
        else:
            regridedFiles = lstInFile

        thisModelFiles = do_stats(variable, validYearList, monthList, regridedFiles, tmpdir, 'stats', '{0}_{1}_{2}'.format(variable,thisModel, rcp), minVar, maxVar )
        if deleteRegrid:
            for ii in regridedFiles: os.remove(ii)

        processedFiles = agregateDict(processedFiles, thisModelFiles)
        gc.collect()

    # if needed process ensembles
    if len(modelList)==1:
        thisLogger.info('>>> 1 model in input: job finished after first averaging round.')
    else:
        thisLogger.info( '>>> Averaging models averages, for each date')
        for idate in processedFiles: # iteration over keys
            thisYear = int(idate[0:4])
            thisMonth= int(idate[4:6])
            thisLogger.info('>> Averaging date {0}'.format(idate))
            listFiles = [x for x in flatten(processedFiles[idate])]
            thisLogger.info('>> averaging files '.format(listFiles))
            returnedList = do_stats('mean_{0}'.format(variable), [thisYear], [thisMonth], listFiles, outdir, 'ensemble', '{0}_{1}'.format(variable, rcp) , minVar, maxVar)

            gc.collect()


# end of file
