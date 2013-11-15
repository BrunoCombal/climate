#!/bin/bash
# \author Bruno Combal
# \date November 2013

function tosEM(){
    rcp=rcp45
    model=modellisttos45.txt
    var=tos
    indir=/databis/cmip5_bis/rcp/rcp4.5/tos/
    outdir=/data/tmp/new_algo/${var}_${rcp}
    tmpdir=/data/tmp/new_algo/tmp_${var}_${rcp}
    bindir='./'

    ${bindir}/make_ensembleMean_tyx.py -v ${var} -path ${indir} -outdir ${outdir} -minVar 1 -maxVar 400 -modelList ${model} -startYear 2030 -endYear 2033 -rcp ${rcp}

    ${bindir}/make_ensembleMean_tyx.py -v ${var} -path ${indir} -outdir ${outdir} -minVar 1 -maxVar 400 -modelList ${model} -startYear 2034 -endYear 2037 -rcp ${rcp}

    ${bindir}/make_ensembleMean_tyx.py -v ${var} -path ${indir} -outdir ${outdir} -minVar 1 -maxVar 400 -modelList ${model} -startYear 2038 -endYear 2039 -rcp ${rcp}
}

function dhm(){
    outdir='/data/tmp/new_algo/dhm/'
    indir='/data/tmp/new_algo/tos_rcp85'
    decad=2030
    climDir='/data/sst/reynolds_climatology/noaa_oist_v2/resized_fitted/sst.ltm.1971-2000_resized.nc'
    maxRealClim='/data/sst/reynolds_climatology/noaa_oist_v2/resized_fitted/max_sst.ltm.1971-2000_resized.nc'
    rmsAtMaxClim='/data/sst/reynolds_climatology/noaa_oist_v2/resized_fitted/rms_at_maxsst_resized.nc'
    modelClimDir='/data/cmip5/rcp/rcp8.5/toshist_ensemble/'
    modelClim='climato_tos_1971_2000_'
    bindir='./'

    ${bindir}/make_dhm.py -o ${outdir} -in ${indir} -decad ${decad} -clim ${climDir} ${maxRealClim} ${rmsAtMaxClim} -modelClim ${modelClimDir} ${modelClim}

}


dhm()