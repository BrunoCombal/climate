#!/bin/bash
# \author Bruno Combal
# \date August 2014


source /usr/local/uvcdat/1.2.0/bin/setup_cdat.sh
indir=/data/tmp/francesco/input/new/tmp
outdir=/data/tmp/francesco/out
modelList=/data/tmp/francesco/models.txt

./make_ensembleMean_tzyx.py -v thetao -minVar 200 -maxVar 330 -monthlist 2 -startYear 2030 -endYear 2030 -path ${indir} -outdir ${outdir} -rcp rcp85 -log /data/tmp/francesco/log/ensemble.log -modellist ${modelList} -keeptmp

