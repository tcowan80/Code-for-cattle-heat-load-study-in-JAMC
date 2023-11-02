# Code-for-cattle-heat-load-study-in-JAMC

This repository contains code to calculate the cattle heat stress indices for AGCD observations over Australia, 
and BARRA-R (version 1) reanalysis. The source input files are needed to successfully run the code.

AGCD data is available at: https://dapds00.nci.org.au/thredds/catalog/zv2/agcd/v1/catalog.html
BARRA data is available at: https://dapds00.nci.org.au/thredds/catalogs/cj37/BARRA/BARRA_R/BARRA_R.html

## calc_barra_hli.py
This derives the Heat Load Index (HLI) from BARRA-R v1 reanalysis.

## calc_barra_ahlu.py
This derives the Accumulated Heat Load Unit (AHLU) from BARRA-R v1 reanalysis.

## calc_obs_agcd_cci_1yr.py
This derives the Cattle Comfort Index (CCI) from AGCD observations (0.05 grid).

## calc_obs_agcd_hli_1yr.py
This derives the Heat Load Index (HLI) from AGCD observations (0.05 grid).

## calc_obs_agcd_thi.py
This derives the Temperature Humidity Index (THI) from AGCD observations (0.05 grid).
