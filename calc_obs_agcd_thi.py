# This code calculates THI from gridded AGCD observations

# it calls in Tmax and 9am/3pm vapor pressure and derives
# daily dewpoint temperature and then calculates THI
# Equation is from https://research.csiro.au/climate/wp-content/uploads/sites/54/2016/03/10_CAF_WorkingPaper10_pdf.pdf


# Author: Tim Cowan (tim.cowan@bom.gov.au)

# Assumption: user has access to the source input variable files

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

for yyyy in np.arange(1990,2018,1):
    
    thi_file = pathlib.Path('temp_humdity_index_day_agcd_' + str(yyyy) + '.nc')
    
    if thi_file.exists ():
        print(str(yyyy) + ' exists')
    
    if not thi_file.exists ():
        	
	# maximum temperature input file
        tasmax_infile = xr.open_dataset('/tasmax_daily.' + str(yyyy) + '.nc',decode_times=True, chunks={"time": 100})
    
        # vapour pressure input file
        vp09_infile = xr.open_dataset('/vprp_09_daily.' + str(yyyy) + '.nc',decode_times=True, chunks={"time": 100})
        vp15_infile = xr.open_dataset('/vprp_15_daily.' + str(yyyy) + '.nc',decode_times=True, chunks={"time": 100}) 
	
        # call in Tx
        tasmax = tasmax_infile.tasmax
        tasmax = xr.DataArray(tasmax, coords=[tasmax.time, tasmax.lat, tasmax.lon], dims=['time', 'lat', 'lon'])

        # call in vapour pressure
        vprp_09 = vp09_infile.vprp_09
        vprp_15 = vp15_infile.vprp_15
        vprp_09 = xr.DataArray(vprp_09, coords=[tasmax.time, tasmax.lat, tasmax.lon], dims=['time', 'lat', 'lon'])
        vprp_15 = xr.DataArray(vprp_15, coords=[tasmax.time, tasmax.lat, tasmax.lon], dims=['time', 'lat', 'lon'])
        vp = (vprp_09+vprp_15)/2
        vp = xr.DataArray(vp, coords=[tasmax.time, tasmax.lat, tasmax.lon], dims=['time', 'lat', 'lon'])
   	
	# calculate f(VP)
        fvp = np.log(vp/6.112)
	
	# calculate dew point temp (Td)
        Td = (234.5*fvp)/(17.67-fvp)
	
	# calculate THI
        thi = tasmax + (0.36*Td) + 41.2
        thi = xr.DataArray(thi, coords=[tasmax.time, tasmax.lat, tasmax.lon], dims=['time', 'lat', 'lon'])

## save the output to netCDF-----------------------------------
        dout_thi = xr.Dataset(
                    {"thi": (("time","lat","lon"), thi.data)},
                        coords={
                            "time": thi.time,
                            "lat": thi.lat,
                            "lon": thi.lon
                        },
        )
	
        dout_thi.to_netcdf('temp_humidity_index_day_agcd_' + str(yyyy) + '.nc')
        print(str(yyyy) + ' done')
