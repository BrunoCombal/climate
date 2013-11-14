climate
=======

Python and bash scripts to process climate projections data.
Codes developped for UNESCO-IOC (International Oceanographic Commission), in the framework of European Commision/7th framework programme "GEOWOW" and GEF "Transboundary Water Assessment Program"

The scripts process CMIP5 models outputs, as obtained from ESGF data portal.
Python scripts require cdms2 library to be installed, this library generaly come with its own python. A way to have the libraries working is to source the initialisation script before calling the python script itself:

$: source /usr/local/uvcdat/1.2.0/bin/setup_cdat.sh

$: ./make_ensembleMean_TOS.py

see wiki page:
https://github.com/BrunoCombal/climate/wiki

latest release:
https://github.com/BrunoCombal/climate/releases/latest