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
# ___________________________
def list_or_tuple(x):
    return isinstance(x, (list, tuple))
# ____________________________
def flatten(sequence, to_expand=list_or_tuple):
    for item in sequence:
        if to_expand(item):
            for subitem in flatten(item, to_expand):
                yield subitem
        else:
            yield item
# ____________________________
# dict{date:[filename]}
def agregateDict(refDict, newDict):
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
        if ikey in refDict.keys(): val.append(refDict[ikey])
        if ikey in newDict.keys(): val.append(newDict[ikey])
        result[ikey] = val

    return result
# ____________________________
def dateTime2Year(datetime):
    result=[]
    for ii in datetime:
        result.append(ii.year)
    return(numpy.array(result), sorted(set(result)))

def makeOutfileName(infile, outdir, prefix, year):
    return('{0}/{1}_{2}_{3}.nc'.format(outdir,prefix,os.path.basename(infile)[:-17], year))

# ____________________________
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
    latAxis.id='latitude'
    latAxis.long_name='Latitude'

    lonAxis = cdms2.createAxis(lon, lon_bnds)
    lonAxis.designateLongitude(True, 360.0)
    lonAxis.units='degrees_east'
    lonAxis.id='longitude'
    lonAxis.long_name='Longitude'

    return((cdms2.createGenericGrid(latAxis, lonAxis, lat_bnds, lon_bnds), latAxis, lonAxis, lat_bnds, lon_bnds))
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

# ___________________________
def do_regrid(variable, indir, lstInFile, outdir, stringBefore):

    createdFiles=[]
    nodata=1.e20
    
    if lstInFile is None:
        print 'No file to process. Return'
        return None

    if len(lstInFile)==0:
        print 'Found no file to process, consider revising search pattern. Return.'
        return None

    (newGrid, latAxis, lonAxis, lat_bnds, lon_bnds) = makeGrid()

    for fileName in lstInFile:
        print 'processing file ', fileName

        thisFile = cdms2.open(fileName)
        data = cdms2.createVariable(thisFile[variable])

        mask = numpy.array(data) < nodata
        regrided = data.regrid(newGrid, missing=nodata, order=thisFile[variable].getOrder(), mask=mask)
        regrided.id=variable

        outfilename = '{0}/{1}{2}'.format(outdir, stringBefore, os.path.basename(fileName))
        createdFiles.append(os.path.basename(outfilename) )
        if os.path.exists(outfilename): os.remove(outfilename)
        outfile = cdms2.open(outfilename, 'w')
        outfile.write(regrided)
        outfile.close()
        thisFile.close()

    return createdFiles
# ___________________________
# for a list of files: open all files, go from date 1 to date 2, compute avg for thisdate, save thisdate
# if a new grid is passed: regrid
def do_stats(variable, indir, lstInFile, outdir, stringBefore, outnameBase, minVar=-1.e20, maxVar=1.e20):
    
    if validYearList is None:
        exitMessage('List of years to process is undefined, edit code. Exit 5.',5)

    createdFiles={}   
    nodata=1.e20

    if lstInFile is None:
        print 'No file to process. Return.'
        return

    if len(lstInFile)==0:
        print 'Found no file to process, consider revising search pattern.'
        return

    print lstInFile

    # open all files
    listFID=[]
    for ifile in lstInFile: listFID.append(cdms2.open('{0}/{1}'.format(indir,ifile), 'r'))
    
    # go through the list of dates, compute ensemble average
    for iyear in validYearList:
        print 'Processing year {0}'.format(iyear)
        for imonth in range(1,13):
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
                if ifile[variable].getTime() is None:
                    if refGrid is None:
                        refGrid = ifile[variable].getGrid()
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
            wtdivide = (accumN < nodata) * (accumN > 0)

            if wtdivide.any():
                accumVar[wtdivide] = accumVar[wtdivide] / accumN[wtdivide]
            accumVar = accumVar/accumN

            # compute std

            # save variables
#            newGrid = cdms2.createGenericGrid(refGrid.getLatitude(), refGrid.getLongitude(), latBounds=refGrid.getLatitude().getBounds(), lonBounds=refGrid.getLongitude().getBounds())
            meanVar = cdms2.createVariable( accumVar.reshape(dims), typecode='f', id='mean_{0}'.format(variable),  fill_value=nodata, attributes=dict(long_name='mean', units=units) )

            meanVar.setGrid(refGrid)

            counter = cdms2.createVariable(accumN.reshape(dims), typecode='i', id='count', fill_value=nodata, attributes=dict(long_name='count', units='None') )
            miniVar = cdms2.createVariable(mini.reshape(dims), typecode='f', id='minimum', fill_value=nodata, attributes=dict(long_name='minimum', units=units) )
            maxiVar = cdms2.createVariable(maxi.reshape(dims), typecode='f', id='maximum', fill_value=nodata, attributes=dict(long_name='maximum', units=units) )
            outfilename = '{0}/{1}_{2}_{3}{4:02}.nc'.format(outdir, stringBefore, outnameBase, iyear, imonth )
            if os.path.exists(outfilename): os.remove(outfilename)
            outfile = cdms2.open(outfilename, 'w')
            outfile.write(meanVar)
            outfile.write(counter)
            outfile.write(miniVar)
            outfile.write(maxiVar)
            outfile.close()

            createdFiles['{0}{1:02}'.format(iyear,imonth)] = os.path.basename(outfilename)

    # close input files
    for ii in listFID: ii.close()

    return(createdFiles)
#___________________________
if __name__=="__main__":

    print 'To make this script properly work, ensure to source cdat setup file first (source /usr/local/uvcdat/VERSION/bin/setup_runtime.sh)'
    print

    variable = None
    indir = None
    tmpdir = None
    outdir = None
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
    for thisModel in modelList:

        pattern=re.compile('{0}_{1}_{2}_{3}_{4}_{5}.nc'.format(variable, 'Omon', thisModel, 'rcp85', 'r.*i.*p.*', '.*') )
        lstInFile=[f for f in glob.glob('{0}/*.nc'.format(indir)) if (os.stat(f).st_size and pattern.match(os.path.basename(f) ) ) ]

        regridedFiles = do_regrid(variable, indir, lstInFile, tmpdir, 'regrid_')

        thisModelFiles = do_stats(variable, tmpdir, regridedFiles, tmpdir, 'stats', '{0}_{1}_{2}'.format(variable,thisModel, 'rcp85') )

        for ii in regridedFiles: os.remove('{0}/{1}'.format(tmpdir, ii))

        processedFiles = agregateDict(processedFiles, thisModelFiles)
        

    print 'Averaging models, for each date:'
    for idate in processedFiles:
        print 'Averaging date ',idate
        print 'Averaging files ',processedFiles[idate]
        listFiles = [ x for x in flatten(processedFiles[idate]) ]
        print do_stats('mean_{0}'.format(variable), tmpdir, listFiles, tmpdir, 'ensemble', '{0}_{1}'.format(variable, 'rcp85') )







