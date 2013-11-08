#!/bin/bash
# \author Bruno Combal
# \date November 2013

function usage(){
    echo "Usage: ${0##*/} frequency dataPath"
    exit 1
}


if [ $# -ne 2 ]; then
    usage
fi

frequency=$1
indir=$2

ls ${indir}/*.nc | sed 's/_rcp.._r.*i.*p.*_[0-9]*-[0-9]*.nc//' | sed 's/.*_'${frequency}'_//' | sort | uniq
