#!/bin/bash
# \author Bruno Combal
# \date November 2013

rcp=rcp85
model=modellist_tos.txt
var=tos
indir=/data/cmip5/rcp/rcp8.5/tos
outdir=/data/tmp/new_algo/${var}_${rcp}
tmpdir=/data/tmp/new_algo/tmp_${var}_${rcp}

./make_ensembleMean_tyx.py -v ${var} -path ${indir} -outdir ${outdir} -minVar 1 -maxVar 400 -modelList ${model} -startYear 2030 -endYear 2033 -rcp ${rcp}

./make_ensembleMean_tyx.py -v ${var} -path ${indir} -outdir ${outdir} -minVar 1 -maxVar 400 -modelList ${model} -startYear 2034 -endYear 2037 -rcp ${rcp}

./make_ensembleMean_tyx.py -v ${var} -path ${indir} -outdir ${outdir} -minVar 1 -maxVar 400 -modelList ${model} -startYear 2038 -endYear 2039 -rcp ${rcp}

