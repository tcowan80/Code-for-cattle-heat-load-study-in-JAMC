## This calculates accumulated heat load units (AHLU)
# from heat load index (HLI) in AGCD observations.

# import python libraries

import os
import sys
import math
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np
import datetime
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import netCDF4 as nc4
import calendar
import cartopy.io.shapereader as shpreader
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import string
import numpy.ma as ma
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import shapefile
import geopandas as gpd
import pathlib

yyyymm1=sys.argv[1]; # year-month 1 of analysis (e.g., 200007 for July 2000)
yyyymm2=sys.argv[2]; # year-month 2 of analysis (e.g., 200106 for June 2001)
typeX=str(sys.argv[3]) ; #either analysis or forecast
hli_lower=int(sys.argv[4]) ; # usually set to 77
hli_upper=int(sys.argv[5]) ; # usually set to 86 for Angus, higher for tropical breeds, no higher than 96.

## example: python calc_barra_ahlu.py 200007 200106 'forecast' 77 86 ->  This calculates the AHLU
# from BARRA from July 2000 to June 2001 for Angus breed of cattle. 

ii=0

# open up HLI input file, change according to saved output from calc_barra_hli.py
hli_infile = xr.open_dataset('heat_load_index_smooth_hour_BARRA_R-'+typeX+'.'+str(yyyymm1)+'-'+str(yyyymm2)+'_halfwind_Aust.nc',decode_times=True, chunks={"time": 100})

# call in heat load index
hli_smooth = hli_infile.hli_smooth
nt,ny,nx = hli_smooth.shape
    
if (ii==0):
    ahlu = np.zeros((nt,ny,nx))
    
for doy in np.arange(0,nt):
    hli_tmp = hli_smooth[doy,:,:]
    excess = np.where(hli_tmp < hli_lower,hli_tmp - hli_lower,0)
    excess = np.where(hli_tmp > hli_upper,hli_tmp - hli_upper,excess)
    excess = np.where(excess < 0,excess/2,excess)	
    
    if (typeX == "an"):
        excess = excess*6
        

    if (doy==0):
        ahlu[doy,:,:] = ahlu[doy,:,:] + excess
    else:
        ahlu[doy,:,:] = ahlu[doy-1,:,:] + excess		    
    ahlu[doy,:,:] = np.where(ahlu[doy,:,:] < 0, 0, ahlu[doy,:,:])

    print(hli_smooth.time.dt.date[doy])

dout_ahlu = xr.Dataset(
            {"ahlu": (("time","latitude","longitude"), ahlu)},
                coords={
                    "time": hli_smooth.time,
                    "latitude": hli_smooth.latitude,
                    "longitude": hli_smooth.longitude
                },
)
dout_ahlu['ahlu'].attrs = {"standard_name": 'accumulated heat load unit',"long_name": 'AHLU for heat load index threshold of '+str(hli_upper)}
dout_ahlu.attrs = {'info':'Accumulated heat load units calculated using method described in Gaughan et al. 2008: \
		     A New Heat Load Index for Feedlot Cattle. \
		     Faculty Papers and Publications in Animal Science.,613.\
		     0.5x conversion factor applied to winds',
		     'id':'BARRA_v1',
		     'history':str(datetime.datetime.now())}

dout_ahlu.to_netcdf('accum_heat_load_unit_'+str(hli_upper)+'_BARRA_R-'+typeX+'.'+str(yyyymm1)+'-'+str(yyyymm2)+'_2m_logwind_Aust.nc',
		      encoding={'ahlu': {'_FillValue':-9999}}, unlimited_dims='time')
