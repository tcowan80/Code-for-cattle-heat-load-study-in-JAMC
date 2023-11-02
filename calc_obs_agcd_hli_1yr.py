# This code calculates Heat load index (HLI) from gridded AGCD observations

# it calls in Tmax/Tmin/Solar/Winds and 9am/3pm vapor pressure and derives
# HLI based on the method described in Gaughan et al. 2008: doi: 10.2527/jas.2007-0305

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
    
    hli_smooth_file = pathlib.Path('heat_load_index_smooth_day_agcd_' + str(yyyy) + '.nc')
    
    if hli_smooth_file.exists ():
        print(str(yyyy) + ' exists')
    
    if not hli_file.exists ():
        
	# solar radiation input file
        rsds_infile = xr.open_dataset('rsds_daily.' + str(yyyy) + '.nc',decode_times=True, chunks={"time": 100})

        # maximum temperature input file
        tasmax_infile = xr.open_dataset('tasmax_daily.' + str(yyyy) + '.nc',decode_times=True, chunks={"time": 100})

        # minimum temperature input file
        tasmin_infile = xr.open_dataset('tasmin_daily.' + str(yyyy) + '.nc',decode_times=True, chunks={"time": 100})
    
        # vapour pressure input file
        vp09_infile = xr.open_dataset('vprp_09_daily.' + str(yyyy) + '.nc',decode_times=True, chunks={"time": 100})
        vp15_infile = xr.open_dataset('vprp_15_daily.' + str(yyyy) + '.nc',decode_times=True, chunks={"time": 100}) 
        
	# wind input file
        wndsp_av_infile = xr.open_dataset('wind_'+str(yyyy)+'.nc', decode_times=True, chunks={"time": 100})
	
	# call in SR
        rsds = rsds_infile.rsds.sel(lat=slice(-44,-10),lon=slice(112,154))
        rsds = xr.DataArray(rsds, coords=[rsds.time, rsds.lat, rsds.lon], dims=['time', 'lat', 'lon'])

        # call in Tx
        tasmax = tasmax_infile.tasmax.sel(lat=slice(-44,-10),lon=slice(112,154))
        tasmax = xr.DataArray(tasmax, coords=[rsds.time, rsds.lat, rsds.lon], dims=['time', 'lat', 'lon'])	
        
	# call in Tn
        tasmin = tasmin_infile.tasmin.sel(lat=slice(-44,-10),lon=slice(112,154))
        tasmin = xr.DataArray(tasmin, coords=[rsds.time, rsds.lat, rsds.lon], dims=['time', 'lat', 'lon'])
	
        # calculate average temperature
        tavg = (tasmax+tasmin)/2
        tavg = xr.DataArray(tavg, coords=[rsds.time, rsds.lat, rsds.lon], dims=['time', 'lat', 'lon'])

        # calculate BGT 
        bgt = (1.33*tavg) - (2.65*(tavg**0.5)) + (3.21*np.log10(rsds+1)) + 3.5
        sbgt = 1/(1+np.exp(((25-bgt)/2.25)))  # smoothed bgt function

        # make BGT data array
        bgt = xr.DataArray(bgt,  coords=[rsds.time, rsds.lat, rsds.lon], dims=['time', 'lat', 'lon'])

        # calculate relative humidity (rh) from vapour pressure
        # need to calculate saturation vapour pressure (svp)
        svp = 6.11*(10**((7.5*tavg)/(237.3+tavg)))
        svp1 = np.exp(1.8096 + ((17.269425*tavg)/(237.3+tavg)))

        # call in vapour pressure
        vprp_09 = vp09_infile.vprp_09.sel(lat=slice(-44,-10),lon=slice(112,154))
        vprp_15 = vp15_infile.vprp_15.sel(lat=slice(-44,-10),lon=slice(112,154))
        vprp_09 = xr.DataArray(vprp_09, coords=[rsds.time, rsds.lat, rsds.lon], dims=['time', 'lat', 'lon'])
        vprp_15 = xr.DataArray(vprp_15, coords=[rsds.time, rsds.lat, rsds.lon], dims=['time', 'lat', 'lon'])
        vp = (vprp_09+vprp_15)/2
        vp = xr.DataArray(vp, coords=[rsds.time, rsds.lat, rsds.lon], dims=['time', 'lat', 'lon'])
        rh = (vp/svp1) * 100 # in %
        rh = xr.DataArray(rh, coords=[rsds.time, rsds.lat, rsds.lon], dims=['time', 'lat', 'lon'])

        # call in wind
        wndsp_av = wndsp_av_infile.wind.sel(latitude=slice(-10,-44),longitude=slice(112,154))   
        wndsp_av = wndsp_av[:,::-1,:] ; # reverse latitudes in wind
        wndsp_av = xr.DataArray(wndsp_av, coords=[rsds.time, rsds.lat, rsds.lon], dims=['time', 'lat', 'lon'])
	
	# calculate heat load index (HLI)
        hli_lt25 = 10.66 + (0.28*rh) + (1.3*bgt) - wndsp_av
        hli_gt25 = 8.62 + (0.38*rh) + (1.55*bgt) - (0.5*wndsp_av) + np.exp(2.4-wndsp_av)
        hli = np.where(bgt<=25, hli_lt25, hli_gt25)
        
	# set any HLI values less than 50 to 50
        hli = np.where(hli<50,50,hli)	
        hli = xr.DataArray(hli, coords=[rsds.time, rsds.lat, rsds.lon], dims=['time', 'lat', 'lon'])
        
	# blended HLI
        hli_smooth = (sbgt*hli_gt25)+((1-sbgt)*hli_lt25)
        
        # set any HLI values less than 50 to 50
        hli_smooth = np.where(hli_smooth<50,50,hli_smooth)
        hli_smooth = xr.DataArray(hli_smooth, coords=[rsds.time, rsds.lat, rsds.lon], dims=['time', 'lat', 'lon'])
        hli_smooth.fillna(-9999)

## save the output to netCDF-----------------------------------
        dout_hli_smooth = xr.Dataset(
                    {"hli_smooth": (("time","lat","lon"), hli_smooth.data)},
                        coords={
                            "time": hli_smooth.time,
                            "lat": hli_smooth.lat,
                            "lon": hli_smooth.lon
                        },
        )
        dout_hli_smooth['hli_smooth'].attrs = {"standard_name": 'smoothed heat load index',"long_name": 'CATTLE HEAT LOAD INDEX'}
        dout_hli_smooth.attrs = {'info':'Heat load index calculated using method described in Mader et al. 2010: \
                                  A comprehensive index for assessing environmental stress in animals. \
                                  J Anim Sci.,88(6):2153-65. doi: 10.2527/jas.2009-2586. \
                                  Epub 2010 Jan 29. Erratum in: J Anim Sci. 2011 Sep;89(9):2955. PMID: 20118427.',
                                  'id':'Australian Gridded Climate Data (AGCD)',
                                 'history':str(datetime.datetime.now())}

        dout_hli_smooth.to_netcdf('heat_load_index_smooth_day_agcd_' + str(yyyy) + '.nc',
                                  encoding={'hli_smooth': {'_FillValue':-9999}}, unlimited_dims='time')
        print(str(yyyy) + ' done')
