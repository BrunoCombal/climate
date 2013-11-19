#!/bin/bash
# \author Bruno Combal
# \date November 2013

function tosEM(){
    rcp=rcp85 #rcp45
    model=modellist_tos.txt # modellisttos45.txt
    var=tos
    indir=/data/cmip5/rcp/rcp8.5/tos/ # /databis/cmip5_bis/rcp/rcp4.5/tos/
    outdir=/data/tmp/new_algo/${var}_${rcp}
    tmpdir=/data/tmp/new_algo/tmp_${var}_${rcp}
    bindir='./'

    ${bindir}/make_ensembleMean_tyx.py -v ${var} -path ${indir} -outdir ${outdir} -minVar 1 -maxVar 400 -modelList ${model} -startYear 2050 -endYear 2053 -rcp ${rcp}

    ${bindir}/make_ensembleMean_tyx.py -v ${var} -path ${indir} -outdir ${outdir} -minVar 1 -maxVar 400 -modelList ${model} -startYear 2054 -endYear 2057 -rcp ${rcp}

    ${bindir}/make_ensembleMean_tyx.py -v ${var} -path ${indir} -outdir ${outdir} -minVar 1 -maxVar 400 -modelList ${model} -startYear 2058 -endYear 2059 -rcp ${rcp}
    # produce some extra month for DHM 4 months rolling window
#    ${bindir}/make_ensembleMean_tyx.py -v ${var} -path ${indir} -outdir ${outdir} -minVar 1 -maxVar 400 -modelList ${model} -startYear 2029 -endYear 2029 -monthlist '10,11,12' -rcp ${rcp}
}

function dhm(){
    outdir='/data/tmp/new_algo/dhm/'
    indir='/data/tmp/new_algo/tos_rcp85'
    inpref='ensemble_tos_rcp85_'
    variable='mean_mean_tos'
    decad=$1
    climDir='/data/sst/reynolds_climatology/noaa_oist_v2/resized_fitted/sst.ltm.1971-2000_resized.nc'
    maxRealClim='/data/sst/reynolds_climatology/noaa_oist_v2/resized_fitted/max_sst.ltm.1971-2000_resized.nc'
    rmsAtMaxClim='/data/sst/reynolds_climatology/noaa_oist_v2/resized_fitted/rms_at_maxsst_resized.nc'
    modelClimDir='/data/cmip5/rcp/rcp8.5/toshist_ensemble/'
    modelClim='climato_tos_1971_2000_'
    bindir='./'

    ${bindir}/make_dhm.py -o ${outdir} -in ${indir} ${inpref} -var ${variable} -decad ${decad} -clim ${climDir} ${maxRealClim} ${rmsAtMaxClim} -modelClim ${modelClimDir} ${modelClim}

}

tosEM
dhm 2050