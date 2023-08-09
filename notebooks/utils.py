import matplotlib.pyplot as plt
from matplotlib import ticker
import pandas as pd
import numpy as np
import glob
import metview as mv
import xarray as xr
import pytz
import datetime
import timezonefinder
from timezonefinder import TimezoneFinder
import warnings


def utc_offset(lat,lon,utc_timestring="2023-08-08 10:00:00"):
    """offset to UTC time in hours
    @param lat latitude of station/grid point
    @param lon longitude of station/grid point
    @param utc_timestring timestring in utc
    @return offset at (lon,lat) and utc_timestring to local time in hours"""
    tf=timezonefinder.TimezoneFinder()
    timezone_str=tf.certain_timezone_at(lat=lat,lng=lon)
    timezone=pytz.timezone(timezone_str)
    dt = datetime.datetime.strptime(utc_timestring, "%Y-%m-%d %H:%M:%S")
    deltat=timezone.utcoffset(dt) #in seconds
    return(deltat.total_seconds()/3600) #in hours

def calc_acc(mod,obs,res=1,nwindow=15):
    """calculates anomaly correlation coefficient (= correlation coefficient after removing seasonal cycle)
        based on moving window (default nwindow=15 [days]), as in LANDVER.py
        @param mod model data
        @param obs obs data
        @param: res resolution (1 h for ERA5 complete, 24 if only one hour is considered)
        @return acc"""
    n=np.shape(mod)[0]
    print(n)
    ap=24/res
    obs_anom=np.empty((n-2*ap*nwindow))
    mod_anom=np.empty((n-2*ap*nwindow))
    for i in range(n+nwindow*ap-1,n-nwindow*ap-1):
        dd=obs[i-nwindow*ap:i+nwindow*ap]
        zz=mod[i-nwindow*ap:i+nwindow*ap]
        obs_anom[i]=(obs[i]-np.nanmean(dd))/(np.nanstd(dd))
        mod_anom[i]=(mod[i]-np.nanmean(zz))/(np.nanstd(zz))
    return(np.corrcoef(mod_anom,obs_anom)[0,1])

def calc_acc_mon(mod,obs):
    """calculates anomaly correlation coefficient (= correlation coefficient after removing seasonal cycle)
       data already per month and hour (without any moving window)
       @param mod model data
       @param obs obs data
       @return acc"""
    mod_anom=(mod-np.nanmean(mod))/np.nanstd(mod)
    obs_anom=(obs-np.nanmean(obs))/np.nanstd(obs)
    return(np.corrcoef(mod_anom[~np.isnan(mod_anom)],obs_anom[~np.isnan(obs_anom)])[0,1])

def mean_diurnal_cycle(vec):
    """calculates mean diurnal cycle of vec"""
    cycle=np.empty((24))
    warnings.filterwarnings("ignore")
    for i in range(24):
        cycle[i]=np.nanmean(vec[i::24])
    return(cycle)

def mean_seasonal_cycle(vec,ndays=366):
    """calculates mean seasonal cycle of vec"""
    cycle=np.empty((12))
    res=int(len(vec)/ndays)
    if ndays==366:
        monlen=[0,31,29,31,30,31,30,31,31,30,31,30,31]
    else:
        monlen=[0,31,28,31,30,31,30,31,31,30,31,30,31]
    warnings.filterwarnings("ignore")
    for i in range(12):
        cycle[i]=np.nanmean(vec[sum(monlen[:(i)]*res):sum(monlen[:i+1]*res)])
    return(cycle)

def mean_seasonal_diurnal_cycle(vec,ndays=366,off=18):
    """calculates mean seasonal and diurnal cycle of vec as matrix with dimension (month,hour)=(12,24)"""
    cycle=np.empty((12,24))
    res=int(len(vec)/ndays)
    if ndays==366:
        monlen=[31,29,31,30,31,30,31,31,30,31,30,31]
    else:
        monlen=[31,28,31,30,31,30,31,31,30,31,30,31]
    warnings.filterwarnings("ignore")
    for i in range(12):
    	for j in range(24):
            tmp=vec[sum(monlen[:(i)]*res):sum(monlen[:i+1]*res)] #select month
            cycle[i,j]=np.nanmean(tmp[off+j::24]) #average per hour
    return(cycle)
    
def obs_daily_av(sh_obs,lh_obs,n=366):
    """calculates daily average flux for icos data (half-hourly, every second value considered)"""
    sh_av_obs=np.zeros(n)
    lh_av_obs=np.zeros(n)
    #model average
    for i in range(0,np.shape(sh_obs)[0]-48,48):
        #print(i)
        sh_av_obs[int((i+48)/48)]=np.nanmean(sh_obs[i:i+48:2])
        lh_av_obs[int((i+48)/48)]=np.nanmean(lh_obs[i:i+48:2])
    return sh_av_obs, lh_av_obs

def era_daily_av(sh_mod,lh_mod,n=366):
    """calculates daily average flux for era5 data (hourly, every value considered)"""
    sh_av_mod=np.zeros(n)
    lh_av_mod=np.zeros(n)
    #model average
    for i in range(0,np.shape(sh_mod)[0]-24,24):
        #print(i)
        sh_av_mod[int((i+24)/24)]=np.nanmean(sh_mod[i:i+24:])
        lh_av_mod[int((i+24)/24)]=np.nanmean(lh_mod[i:i+24:])
    return sh_av_mod, lh_av_mod

def deacc_24h(mod):
    """deaccumulates model data hourly"""
    for i in reversed(range(23)):
        mod[i::24]=mod[i::24]-mod[i-1::24]
    return(mod)
    
def acc_deacc(mod):
    """returns the 12 UTC value corresponding to accumulating 1h values and deacc on 6h"""
    mod12=mod[7::24]+mod[8::24]+mod[9::24]+mod[10::24]+mod[11::24]+mod[12::24]
    return(mod12/6)

def abline(slope, intercept,col="gray",lty="--"):
    """plots abline from slope and intercept as in R"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, lty,color=col)
