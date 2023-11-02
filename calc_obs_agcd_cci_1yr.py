# This code calculates Cattle Comfort Index (CCI) from gridded AGCD observations

# it calls in Tmax/Tmin/Solar/Winds and 9am/3pm vapor pressure and derives
# CCI based on the method described in Mader et al. 2010: doi: 10.2527/jas.2009-2586

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
    
    cci_file = pathlib.Path('cattle_comfort_index_day_agcd_' + str(yyyy) + '.nc')
    
    if cci_file.exists ():
        print(str(yyyy) + ' exists')
    
    if not cci_file.exists ():
        
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
        wndsp_av = wndsp_av[:,::-1,:]  
        wndsp_av = xr.DataArray(wndsp_av, coords=[rsds.time, rsds.lat, rsds.lon], dims=['time', 'lat', 'lon'])

        # calculate temperature and humidity correction factor 
        eq1_p1 = np.exp((0.00182*rh+1.8*(1e-5)*tavg*rh))
        eq1_p2 = (0.000054*(tavg**2)+0.00192*tavg-0.0246)
        eq1_p3 = rh-30
        tavg_rh_cf = eq1_p1*eq1_p2*eq1_p3

        # calculate wind speed correction factor
        xx = (2.26*wndsp_av + 0.33)**(-2)
        b = 0.3
        eq2_p1 = ((2.9+1.14*1e-6*(wndsp_av**2.5) - (np.log(xx)/np.log(b))))
        eq2_p2 = (1/((2.26*wndsp_av + 0.23)**0.45)) * (eq2_p1)
        eq2_p3 = np.exp(eq2_p2)
        wndsp_av_cf = (-6.56/eq2_p3) - 0.00566*(wndsp_av**2) + 3.33 

        # calculate radiation correction factor
        rsds_cf = (0.0076*rsds)-(0.00002*rsds*tavg)+(0.00005*(tavg**2)*(rsds**0.5))+(0.1*tavg)-2

	# calculate CCI
        cci = tavg + tavg_rh_cf + wndsp_av_cf + rsds_cf

        # make data arrays for CCI and correction factors
        cci = xr.DataArray(cci,  coords=[rsds.time, rsds.lat, rsds.lon], dims=['time', 'lat', 'lon'])
        tavg_rh_cf = xr.DataArray(tavg_rh_cf,  coords=[rsds.time, rsds.lat, rsds.lon], dims=['time', 'lat', 'lon'])
        wndsp_av_cf = xr.DataArray(wndsp_av_cf,  coords=[rsds.time, rsds.lat, rsds.lon], dims=['time', 'lat', 'lon'])
        rsds_cf = xr.DataArray(rsds_cf,  coords=[rsds.time, rsds.lat, rsds.lon], dims=['time', 'lat', 'lon'])
        
        cci.fillna(1.e+20)
        tavg.fillna(1.e+20)
        tavg_rh_cf.fillna(1.e+20)
        wndsp_av_cf.fillna(1.e+20)
        rsds_cf.fillna(1.e+20)
			        
## save the output, including the RH, wind, and solar correction factors, and the CCI

#-----Tavg & RH correction factor -------------------------------------------------
        dout_tavg_rh_cf = xr.Dataset(
		    {"tavg_rh_cf": (("time","lat","lon"), tavg_rh_cf.data)},
			coords={
			    "time": tavg_rh_cf.time,
			    "lat": tavg_rh_cf.lat,
			    "lon": tavg_rh_cf.lon
			},
	)
        dout_tavg_rh_cf['tavg_rh_cf'].attrs = {"long_name": 'temperature and rel. humidity correction factor to CCI', "units":'degrees_C'}
        dout_tavg_rh_cf.attrs = {'id':'Australian Gridded Climate Data (AGCD)',
				 'publisher_name':'Australian Bureau of Meteorology',
				 'copyright':'(C) Copyright Commonwealth of Australia (2017), Bureau of Meteorology (ABN 92637 533532), see http://www.bom.gov.au/other/copyright.shtml for terms and conditions of reuse',
				 'history':str(datetime.datetime.now())}

#-----Wind speed correction factor -------------------------------------------------
        dout_wndsp_av_cf = xr.Dataset(
                    {"wndsp_av_cf": (("time","lat","lon"), wndsp_av_cf.data)},
                        coords={
                            "time": wndsp_av_cf.time,
                            "lat": wndsp_av_cf.lat,
                            "lon": wndsp_av_cf.lon
                        },
        )
        dout_wndsp_av_cf['wndsp_av_cf'].attrs = {"long_name": 'wind speed (wendy winds) correction factor to CCI', "units":'degrees_C'}
        dout_wndsp_av_cf.attrs = {'id':'Australian Gridded Climate Data (AGCD)',
                                 'publisher_name':'Australian Bureau of Meteorology',
                                 'copyright':'(C) Copyright Commonwealth of Australia (2017), Bureau of Meteorology (ABN 92637 533532), see http://www.bom.gov.au/other/copyright.shtml for terms and conditions of reuse',
                                 'history':str(datetime.datetime.now())}

#-----Solar correction factor -------------------------------------------------
        dout_rsds_cf = xr.Dataset(
		    {"rsds_cf": (("time","lat","lon"), rsds_cf.data)},
			coords={
			    "time": rsds_cf.time,
			    "lat": rsds_cf.lat,
			    "lon": rsds_cf.lon
			},
	)
        dout_rsds_cf['rsds_cf'].attrs = {"long_name": 'total radiation correction factor to CCI', "units":'degrees_C'}
        dout_rsds_cf.attrs = {'id':'Australian Gridded Climate Data (AGCD)',
				 'publisher_name':'Australian Bureau of Meteorology',
				 'copyright':'(C) Copyright Commonwealth of Australia (2017), Bureau of Meteorology (ABN 92637 533532), see http://www.bom.gov.au/other/copyright.shtml for terms and conditions of reuse',
				 'history':str(datetime.datetime.now())}

#-----CCI -------------------------------------------------	
        dout_cci = xr.Dataset(
                    {"cci": (("time","lat","lon"), cci.data)},
                        coords={
                            "time": cci.time,
                            "lat": cci.lat,
                            "lon": cci.lon
                        },
        )
        dout_cci['cci'].attrs = {"standard_name": 'comprehensive climate index',"long_name": 'CATTLE COMFORT INDEX', "units":'degrees_C'}
        dout_cci.attrs = {'info':'Cattle comfort index calculated using method described in Mader et al. 2010: \
                                  A comprehensive index for assessing environmental stress in animals. \
                                  J Anim Sci.,88(6):2153-65. doi: 10.2527/jas.2009-2586. \
                                  Epub 2010 Jan 29. Erratum in: J Anim Sci. 2011 Sep;89(9):2955. PMID: 20118427.',
                                  'id':'Australian Gridded Climate Data (AGCD)',
                                 'publisher_name':'Australian Bureau of Meteorology',
                                 'copyright':'(C) Copyright Commonwealth of Australia (2017), Bureau of Meteorology (ABN 92637 533532), see http://www.bom.gov.au/other/copyright.shtml for terms and conditions of reuse',
                                 'history':str(datetime.datetime.now())}


#-----Save netCDFs -------------------------------------------------		
        dout_tavg_rh_cf.to_netcdf('temperature_humidity_correction_factor_day_agcd_' + str(yyyy) + '.nc',
				  encoding={'tavg_rh_cf': {'_FillValue':1.e+20}}, unlimited_dims='time')
        dout_wndsp_av_cf.to_netcdf('windspeed_correction_factor_day_agcd_' + str(yyyy) + '.nc',
                                  encoding={'wndsp_av_cf': {'_FillValue':1.e+20}}, unlimited_dims='time')
        dout_rsds_cf.to_netcdf('total_radiation_correction_factor_day_agcd_' + str(yyyy) + '.nc',
				  encoding={'rsds_cf': {'_FillValue':1.e+20}}, unlimited_dims='time')
        dout_cci.to_netcdf('cattle_comfort_index_day_agcd_' + str(yyyy) + '.nc',
                                  encoding={'cci': {'_FillValue':1.e+20}}, unlimited_dims='time')
        print(str(yyyy) + ' done')
