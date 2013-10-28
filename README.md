climate
=======

Python and bash scripts to process climate projections data.
Codes developped for UNESCO-IOC (International Oceanographic Commission), in the framework of European Commision/7th frameworj programme "GEOWOW" and GEF "Transfoundary Water Assessment Program"

The scripts process CMIP5 models outputs, as obtained from ESGF data portal.
Python scripts require cdms2 library to be installed, this library generaly come with its own python. A way to have the libraries working is to source the initialisation script before calling the python script itself:
$: source /usr/local/uvcdat/1.2.0/bin/setup_cdat.sh
$: ./make_ensembleMean_TOS.py


- make_ensembleMean_TOS.py
a script to compute an ensemble mean (and min, max, std) from CMIP5/ESGF datasets (Sea Surface Temperature, named TOS, Temperature Of Surface, in CMIP5 naming convention).

