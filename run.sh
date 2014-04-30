#!/bin/bash
# \author Bruno Combal
# \date November 2013

# a series of function to call the different cliamte processing
# caution: this script is giving as an example, there is no input parameters tests.
# Read carefuly before using.

# _________________________________
# \brief Computes tos ensemble means
# $1 decade
# $2 rcp85, rcp45
function omlmaxEM(){
    decade=$1
    rcp=$2
    var=omlmax

    if [ ${rcp} = 'rcp85' ]; then
	rcp=rcp85
	indir=/datater/cmip5/rcp/rcp8.5/${var}/
	model=${indir}/modellist_omlmax.txt
    elif [ ${rcp} = 'rcp45' ]; then
	rcp=rcp45
	echo "rcp 4.5 has not been expected. Exit."
	exit 1
    else
	echo "Function omlmaxEM: wrong input parameter '$2'="${rcp}". Exit(1)."
	exit 1
    fi

    if [ ! -e ${model} ]; then
	ls ${model}
	echo "List of models does not exist ($model). Exit 2."
	exit 2
    fi

    outdir=/datater/tmp/new_algo/${var}_${rcp}
    tmpdir=/datater/tmp/new_algo/tmp_${var}_${rcp}
    mkdir -p ${outdir}
    mkdir -p ${tmpdir}
    bindir=${HOME}/github/climate/

    dStart=${decade}
    dEnd=$((decade + 3))
    echo "Processing ensembleMean from "${dStart}" to "${dEnd}
    ${bindir}/make_ensembleMean_tyx.py -v ${var} -path ${indir} -outdir ${outdir} -minVar 263 -maxVar 320 -modelList ${model} -startYear ${dStart} -endYear ${dEnd} -rcp ${rcp}

    dStart=$((dEnd + 1))
    dEnd=$((dStart + 3))
    echo "Processing ensembleMean from "${dStart}" to "${dEnd}
    ${bindir}/make_ensembleMean_tyx.py -v ${var} -path ${indir} -outdir ${outdir} -minVar 263 -maxVar 320 -modelList ${model} -startYear ${dStart} -endYear ${dEnd} -rcp ${rcp}

    dStart=$((dEnd + 1))
    dEnd=$((decade + 9))
    echo "Processing ensembleMean from "${dStart}" to "${dEnd}
    ${bindir}/make_ensembleMean_tyx.py -v ${var} -path ${indir} -outdir ${outdir} -minVar 263 -maxVar 320 -modelList ${model} -startYear ${dStart} -endYear ${dEnd} -rcp ${rcp}

}
# ________________________________
# \brief compute TOS ensemble mean  
function tosEM(){
    decade=$1
    rcp=$2
    var=tos

    if [ ${rcp} = 'rcp85' ]; then
	rcp=rcp85
	model=modellist_tos.txt
	indir=/data/cmip5/rcp/rcp8.5/tos/
    elif [ ${rcp} = 'rcp45' ]; then
	rcp=rcp45
	model=modellist_tos45.txt
	indir=/databis/cmip5_bis/rcp/rcp4.5/tos/
    else
	echo "Function tosEM: wrong input parameter '$2'="${rcp}". Exit(1)."
	exit 1
    fi

    if [ ! -e "${model}" ] ; then
	echo "Missing model list file "${model}" Exit(1)."
	exit 1
    fi

    outdir=/data/tmp/new_algo/${var}_${rcp}
    tmpdir=/data/tmp/new_algo/tmp_${var}_${rcp}
    mkdir -p ${outdir}
    mkdir -p ${tmpdir}
    bindir='./'

    dStart=${decade}
    dEnd=$((decade + 3))
    echo "Processing ensembleMean from "${dStart}" to "${dEnd}
    ${bindir}/make_ensembleMean_tyx.py -v ${var} -path ${indir} -outdir ${outdir} -minVar 263 -maxVar 320 -modelList ${model} -startYear ${dStart} -endYear ${dEnd} -rcp ${rcp}

    dStart=$((dEnd + 1))
    dEnd=$((dStart + 3))
    echo "Processing ensembleMean from "${dStart}" to "${dEnd}
    ${bindir}/make_ensembleMean_tyx.py -v ${var} -path ${indir} -outdir ${outdir} -minVar 263 -maxVar 320 -modelList ${model} -startYear ${dStart} -endYear ${dEnd} -rcp ${rcp}

    dStart=$((dEnd + 1))
    dEnd=$((decade + 9))
    echo "Processing ensembleMean from "${dStart}" to "${dEnd}
    ${bindir}/make_ensembleMean_tyx.py -v ${var} -path ${indir} -outdir ${outdir} -minVar 263 -maxVar 320 -modelList ${model} -startYear ${dStart} -endYear ${dEnd} -rcp ${rcp}
    
    # produce some extra month for DHM 4 months rolling window
    dStart=$((decade - 1))
    dEnd=$((decade - 1))
    ${bindir}/make_ensembleMean_tyx.py -v ${var} -path ${indir} -outdir ${outdir} -minVar 263 -maxVar 320 -modelList ${model} -startYear ${dStart} -endYear ${dEnd} -monthlist '10,11,12' -rcp ${rcp}
}

# ____________________________________
# \brief Computes thetao ensemble means
function thetaoEM(){
    rcp=rcp85 #rcp45
    model=modellist_thetao_rcp85.txt
    var=thetao
    indir=/databis/cmip5_bis/rcp/rcp8.5/thetao
    outdir=/data/tmp/new_algo/${var}_${rcp}
    tmpdir=/data/tmp/new_algo/tmp_${var}_${rcp}
    bindir='./'

    ${bindir}/make_ensembleMean_tzyx.py -v ${var} -path ${indir} -outdir ${outdir} -minVar 263 -maxVar 320 -modelList ${model} -startYear 2030 -endYear 2031 -rcp ${rcp} -resolution 1

}
# _____________________________________
# \brief Computes DHM (degree heating month)
function dhm(){
    outdir='/data/tmp/new_algo/dhm_rcp'${2}'/'
    indir='/data/tmp/new_algo/tos_'${2}
    inpref='ensemble_tos_'${2}'_'
    variable='mean_mean_tos'
    decad=$1
    climDir='/data/sst/reynolds_climatology/noaa_oist_v2/resized_fitted/sst.ltm.1971-2000_resized.nc'
    maxRealClim='/data/sst/reynolds_climatology/noaa_oist_v2/resized_fitted/max_sst.ltm.1971-2000_resized.nc'
    rmsAtMaxClim='/data/sst/reynolds_climatology/noaa_oist_v2/resized_fitted/rms_at_maxsst_resized.nc'
    modelClimDir='/data/cmip5/rcp/toshist_ensemble/'
    modelClim='climato_tos_1971_2000_'
    bindir='./'

    ${bindir}/make_dhm.py -o ${outdir} -in ${indir} ${inpref} -var ${variable} -decad ${decad} -clim ${climDir} ${maxRealClim} ${rmsAtMaxClim} -modelClim ${modelClimDir} ${modelClim}

}

# ___________________________________________
# \brief computes model trends
# $1: variable, eg. 'tos'
# $2: trendType, eg. 'esm'
function trends(){
    variable=$1
    trendType=$2
    indir='/data/cmip5/rcp/tos_controlrun'
    modelList=${indir}/lst_${variable}_ctrlrun_${trendType}_edited.txt
    outdir=/data/tmp/new_algo/control_run/${tos}_${trendType}
    degree=1

    ./make_modelTrend.py -deg 1 -o ${outdir} -p $indir -v ${variable} -trendType ${trendType} -modellist ${modelList} -annualAVG yes
}

# ___________________________________________
# \brief Main
# use thise section to call the functions

# uv-cdat library must first be sourced
source /usr/local/uvcdat/1.2.0/bin/setup_cdat.sh

#for ii in $(seq 2010 10 2080 )
#do
#    omlmaxEM $ii rcp85
#done
#trends 'tos' 'esm'

#tosEM 2010 rcp85

#tosEM 2040 rcp85
#dhm 2040 rcp85

#tosEM 2030 rcp45
#tosEM 2040 rcp45
#tosEM 2050 rcp45
dhm 2030 rcp85
dhm 2040 rcp85
dhm 2050 rcp85
#thetaoEM