# Calculates heat load index (HLI) 
# from BARRA_v1 reanalysis

# import python libraries

import os
import sys
import math
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
from pandas.tseries.offsets import DateOffset
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


yyyy1 = sys.argv[1] ; # year of analysis start (e.g., 2000)
mm1 = sys.argv[2] ; # month 1 start period of analysis (e.g., July = 07)
yyyy2 = sys.argv[3] ; # year of analysis end (e.g., 2001)
mm2 = sys.argv[4] # month 2 end period of analysis (e.g., June = 06)
typeX = str(sys.argv[5]) ; #either analysis or forecasts

## example: python calc_barra_hli.py 2000 07 2001 06 'forecast'  ->  This calculates the HLI
# from BARRA from July 2000 to June 2001. 


# input file directories on the GADI (NCI)
solar_dir = '/g/data/cj37/BARRA/BARRA_R/v1/'+typeX+'/slv/av_swsfcdown'
solar_file_list=[]

temp_dir = '/g/data/cj37/BARRA/BARRA_R/v1/'+typeX+'/slv/av_temp_scrn'
temp_file_list=[]

dewpt_dir = '/g/data/cj37/BARRA/BARRA_R/v1/'+typeX+'/slv/dewpt_scrn'
dewpt_file_list=[]

uwnd_dir = '/g/data/cj37/BARRA/BARRA_R/v1/'+typeX+'/slv/av_uwnd10m'
uwnd_file_list=[]

vwnd_dir = '/g/data/cj37/BARRA/BARRA_R/v1/'+typeX+'/slv/av_vwnd10m'
vwnd_file_list=[]


if (typeX == 'forecast'):
    sh_type = 'fc'
elif (typeX == 'analysis'):
    sh_type = 'an'
else:
    print("have to enter either forecast or analysis")
    sys.exit()       

## Assumption: User has access to the source input variable files and directories


# period of interest
if (typeX == 'forecast'):
    pr = pd.date_range(start=str(yyyy1)+'-'+str(mm1).zfill(2)+'-01',end=str(yyyy2)+'-'+str(mm2).zfill(2)+'-01', freq='1H') + DateOffset(hours=0.5)
elif (typeX == 'analysis'):
    pr = pd.date_range(start=str(yyyy1)+'-'+str(mm1).zfill(2)+'-01',end=str(yyyy2)+'-'+str(mm2).zfill(2)+'-01', freq='6H')
else:
    print("have to enter either forecast or analysis")
    sys.exit()    

yyyymmdd = pr.strftime('%Y%m%d')
yyyymm = pr.strftime('%Y%m')
hhmm = pr.strftime('%H%M')


solar_file_list = os.popen('ls '+solar_dir+'/'+str(yyyy1)+'/'+str(mm1).zfill(2)+'/av_swsfcdown*').read()
solar_file_list = solar_file_list.split('\n')
solar_file_list_arry = np.asarray(solar_file_list)
        
temp_file_list = os.popen('ls '+temp_dir+'/'+str(yyyy1)+'/'+str(mm1).zfill(2)+'/av_temp_scrn*').read()
temp_file_list = temp_file_list.split('\n')
temp_file_list_arry = np.asarray(temp_file_list)

dewpt_file_list = os.popen('ls '+dewpt_dir+'/'+str(yyyy1)+'/'+str(mm1).zfill(2)+'/dewpt_scrn*').read()
dewpt_file_list = dewpt_file_list.split('\n')
dewpt_file_list_arry = np.asarray(dewpt_file_list)

uwnd_file_list = os.popen('ls '+uwnd_dir+'/'+str(yyyy1)+'/'+str(mm1).zfill(2)+'/av_uwnd10m*').read()
uwnd_file_list = uwnd_file_list.split('\n')
uwnd_file_list_arry = np.asarray(uwnd_file_list)

vwnd_file_list = os.popen('ls '+vwnd_dir+'/'+str(yyyy1)+'/'+str(mm1).zfill(2)+'/av_vwnd10m*').read()
vwnd_file_list = vwnd_file_list.split('\n')
vwnd_file_list_arry = np.asarray(vwnd_file_list)

solar_infile = xr.open_mfdataset(solar_file_list_arry[0:-1],concat_dim='time',combine="nested",chunks={"time": 100},decode_times=True)
temp_infile = xr.open_mfdataset(temp_file_list_arry[0:-1],concat_dim='time',combine="nested",chunks={"time": 100},decode_times=True)
dewpt_infile = xr.open_mfdataset(dewpt_file_list_arry[0:-1],concat_dim='time',combine="nested",chunks={"time": 100},decode_times=True)
uwnd_infile = xr.open_mfdataset(uwnd_file_list_arry[0:-1],concat_dim='time',combine="nested",chunks={"time": 100},decode_times=True)
vwnd_infile = xr.open_mfdataset(vwnd_file_list_arry[0:-1],concat_dim='time',combine="nested",chunks={"time": 100},decode_times=True)

# input solar
av_swsfcdown = solar_infile.av_swsfcdown.sel(longitude=slice(112,155),latitude=slice(-45,-9))
av_swsfcdown["time"] = pr[0:-1]
print("Solar loaded")

# input temperature
av_temp_scrn = temp_infile.av_temp_scrn.sel(longitude=slice(112,155),latitude=slice(-45,-9))
av_temp_scrn = av_temp_scrn-273.15 # convert K to degC
av_temp_scrn = xr.DataArray(av_temp_scrn, coords=[av_swsfcdown.time, av_swsfcdown.latitude, av_swsfcdown.longitude], dims=['time', 'latitude', 'longitude'])
print("Screen temperature loaded")

# input dewpt tempetaure
dewpt_scrn = dewpt_infile.dewpt_scrn.sel(longitude=slice(112,155),latitude=slice(-45,-9))
dewpt_scrn = dewpt_scrn-273.15 # convert K to degC
dewpt_scrn = xr.DataArray(dewpt_scrn, coords=[av_swsfcdown.time, av_swsfcdown.latitude, av_swsfcdown.longitude], dims=['time', 'latitude', 'longitude'])
print("Dew-pt temperature loaded")

# calculate RH
rh = ((np.exp(1.8096+((17.2694*dewpt_scrn)/(237.3+dewpt_scrn))))/(np.exp(1.8096+((17.2694*av_temp_scrn)/(237.3+av_temp_scrn))))) * 100
rh = xr.DataArray(rh, coords=[av_swsfcdown.time, av_swsfcdown.latitude, av_swsfcdown.longitude], dims=['time', 'latitude', 'longitude'])
print("Relative humidity calculated")

# input winds
#av_uwnd10m = uwnd_infile.av_uwnd10m.sel(longitude=slice(112,155),latitude=slice(-45,-9)) * (4.87/(np.log((67.8*uwnd_infile.av_uwnd10m.sel(longitude=slice(112,155),latitude=slice(-45,-9))) - 5.42))) # for cattle heights # for cattle heights
av_uwnd10m = uwnd_infile.av_uwnd10m.sel(longitude=slice(112,155),latitude=slice(-45,-9)) * 0.5 # for cattle heights
av_uwnd10m = xr.DataArray(av_uwnd10m, coords=[av_swsfcdown.time, av_swsfcdown.latitude, av_swsfcdown.longitude], dims=['time', 'latitude', 'longitude'])

#av_vwnd10m = vwnd_infile.av_vwnd10m.sel(longitude=slice(112,155),latitude=slice(-45,-9.1)) * (4.87/(np.log((67.8*vwnd_infile.av_vwnd10m.sel(longitude=slice(112,155),latitude=slice(-45,-9.1))) - 5.42))) # for cattle heights
av_vwnd10m = vwnd_infile.av_vwnd10m.sel(longitude=slice(112,155),latitude=slice(-45,-9.1)) * 0.5 # for cattle heights
av_vwnd10m = xr.DataArray(av_vwnd10m, coords=[av_swsfcdown.time, av_swsfcdown.latitude, av_swsfcdown.longitude], dims=['time', 'latitude', 'longitude'])

# calculate wind speed
ws = np.sqrt((av_uwnd10m**2) + (av_vwnd10m**2))
ws = xr.DataArray(ws, coords=[av_swsfcdown.time, av_swsfcdown.latitude, av_swsfcdown.longitude], dims=['time', 'latitude', 'longitude'])
print("Wind speed calculated")

# calculate black globe temperature
bgt = (1.33*av_temp_scrn) - (2.65*(av_temp_scrn**0.5)) + (3.21*np.log10(av_swsfcdown+1)) + 3.5
sbgt = 1/(1+np.exp(((25-bgt)/2.25)))  # smoothed bgt function
print("Black globe temperature calculated")

# make data array
bgt = xr.DataArray(bgt, coords=[av_swsfcdown.time, av_swsfcdown.latitude, av_swsfcdown.longitude], dims=['time', 'latitude', 'longitude'])

# calculate heat load index (HLI)
hli_lt25 = 10.66 + (0.28*rh) + (1.3*bgt) - ws
hli_gt25 = 8.62 + (0.38*rh) + (1.55*bgt) - (0.5*ws) + np.exp(2.4-ws)
hli = np.where(bgt<=25, hli_lt25, hli_gt25)

# set any HLI values less than 50 to 50
hli = np.where(hli<50,50,hli)	
hli = xr.DataArray(hli, coords=[av_swsfcdown.time, av_swsfcdown.latitude, av_swsfcdown.longitude], dims=['time', 'latitude', 'longitude'])
print("Heat load index calculated")

# blended HLI
hli_smooth = (sbgt*hli_gt25)+((1-sbgt)*hli_lt25)

# set any HLI values less than 50 to 50
hli_smooth = np.where(hli_smooth<50,50,hli_smooth)
hli_smooth = xr.DataArray(hli_smooth, coords=[av_swsfcdown.time, av_swsfcdown.latitude, av_swsfcdown.longitude], dims=['time', 'latitude', 'longitude'])
print("Smoothed heat load index calculated")


# save output to netCDF (option to save BGT and HLI unfiltered)
#dout_bgt = xr.Dataset(
#	     {"bgt": (("time","lat","lon"), bgt.data)},
#		 coords={
#		     "time": bgt.time,
#		     "lat": bgt.lat,
#		     "lon": bgt.lon
#		 },
#)

#dout_hli = xr.Dataset(
#	     {"hli": (("time","lat","lon"), hli.data)},
#		 coords={
#		     "time": hli.time,
#		     "lat": hli.lat,
#		     "lon": hli.lon
#		 },
#)

dout_hli_smooth = xr.Dataset(
            {"hli_smooth": (("time","latitude","longitude"), hli_smooth.data)},
                coords={
                    "time": hli_smooth.time,
                    "latitude": hli_smooth.latitude,
                    "longitude": hli_smooth.longitude
                },
)

dout_hli_smooth['hli_smooth'].attrs = {"standard_name": 'smoothed heat load index',"long_name": 'CATTLE HEAT LOAD INDEX'}
dout_hli_smooth.attrs = {'info':'Heat load index calculated using method described in Mader et al. 2010: \
                                  A comprehensive index for assessing environmental stress in animals. \
                                  J Anim Sci.,88(6):2153-65. doi: 10.2527/jas.2009-2586. \
                                  Epub 2010 Jan 29. Erratum in: J Anim Sci. 2011 Sep;89(9):2955. PMID: 20118427.',
                                  'id':'Australian Gridded Climate Data (AGCD)',
                                 'history':str(datetime.datetime.now())}

#dout_bgt.to_netcdf('black_globe_temp_day_agcd_' + str(yyyy) + '.nc')
#dout_hli.to_netcdf('heat_load_index_day_agcd_' + str(yyyy) + '.nc')
dout_hli_smooth.to_netcdf('heat_load_index_smooth_hour_BARRA_R-'+sh_type+'.'+yyyymm[0]+'_halfwind_Aust.nc')
print(yyyymm[0],"done")
