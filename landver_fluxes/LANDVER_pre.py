#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from importlib.metadata import files
import os
import numpy as np
import pandas as pd

# from Options_namelist import *
import pylab as pl
from Common_functions import showProgress, ensure_dir, get_parser
import glob, re
import metview as mv
import xarray as xr
from requests import get
from pandas import json_normalize

numbers = re.compile(r"(\d+)")


def get_elevation(lat=None, long=None):
    """
    script for returning elevation in m from lat, long
    """
    if lat is None or long is None:
        return None

    query = "https://api.open-elevation.com/api/v1/lookup" f"?locations={_lat},{_long}"

    # Request with a timeout for slow responses
    r = get(query, timeout=60)

    # Only get the json response in case of 200 or 201
    if r.status_code == 200 or r.status_code == 201:
        elevation = json_normalize(r.json(), "results")["elevation"].values[0]
    else:
        elevation = None
    return elevation


def numericalSort(value):
    parts = numbers.split(value)
    return parts


def preprocessData(cfg):
    """ preprocessing of analyses
    @param: cfg config class
    @return: no return, writes a csv file of preprocessed data for every experiment and every station (per network) in the pre_processed_data directory, 
            output file format: yyyymmddhh sh lh, output file name: fluxes_lat,lon [station]
    """
    # Ensure directories exist
    ensure_dir(cfg.output_dir)
    ensure_dir(cfg.output_dir + "/analysis_data")
    ensure_dir(cfg.output_dir + "/analysis_data/pre_processed_data")
    ensure_dir(cfg.output_dir + "/analysis_data/raw_data")

    grib_dir = cfg.output_dir + "/analysis_data/raw_data/"

    # Data times
    data_range = pd.date_range(cfg.analysis_period[0], cfg.analysis_period[1])
    years = pd.date_range(
        str(data_range[0].year), str(data_range[-1].year + 1), freq="A"
    )
    """
    if cfg.extract_SM or cfg.extract_ST: #extract from MARS - not yet adapted for fluxes

        # Extract grib files from MARS:
        for EXP in range(len(cfg.EXPVER)):
            n_tasks = data_range.size
            task = 0
            Time_freq = cfg.TIME[EXP].split("/")
            print("\nExtracting data for experiment " + cfg.EXPVER[EXP])
            print(data_range)
            for year in set(data_range.year):
                SM_date_retrieve = []
                SM_time_retrieve = []
                ST_date_retrieve = []
                ST_time_retrieve = []
                print(year)
                grib_dir = (
                    cfg.output_dir
                    + "/analysis_data/raw_data/"
                    + cfg.EXPVER[EXP]
                    + "_"
                    + cfg.CLASS[EXP]
                    + "/"
                    + str(year)
                    + "/"
                )
                ensure_dir(grib_dir)

                for d in data_range[data_range.year == year]:
                    YYYYMMDD = "%02d" % (d.year) + "%02d" % (d.month) + "%02d" % (d.day)

                    for t in Time_freq:
                        HH = t[:2]

                        if not os.path.isfile(
                            grib_dir + "/SM_" + YYYYMMDD + HH + ".grib"
                        ) and (YYYYMMDD not in SM_date_retrieve):
                            SM_date_retrieve.append(YYYYMMDD)

                        if HH not in SM_time_retrieve:
                            SM_time_retrieve.append(HH)

                        if not os.path.isfile(
                            grib_dir + "/ST_" + YYYYMMDD + HH + ".grib"
                        ) and (YYYYMMDD not in ST_date_retrieve):
                            ST_date_retrieve.append(YYYYMMDD)

                        if HH not in ST_time_retrieve:
                            ST_time_retrieve.append(HH)

                if cfg.extract_SM:

                    print(SM_date_retrieve)

                    if len(SM_date_retrieve) > 0:
                        # Retrieve all SM data for dates not yet retrieved

                        if (cfg.GRID[EXP]=="av") or (int(cfg.GRID[EXP][1:5])<2559):
                            SM_data = mv.retrieve(
                                class_=cfg.CLASS[EXP],
                                type=cfg.TYPE[EXP],
                                expver=cfg.EXPVER[EXP],
                                stream=cfg.STREAM[EXP],
                                anoffset=cfg.ANOFFSET[EXP],
                                levtype=cfg.LEVTYPE[EXP],
                                param=["39.128", "40.128", "41.128", "42.128"],
                                date=SM_date_retrieve,
                                time=SM_time_retrieve,
                                step=cfg.STEP[EXP],
                                grid=cfg.GRID[EXP],
                                domain=cfg.DOMAIN[EXP],
                            )
                            print(SM_data)

                            # Write files for each date
                            n_tasks = len(SM_date_retrieve)
                            task = 0

                            print(
                                "\nExtracting SM data for experiment "
                                + cfg.EXPVER[EXP]
                                + " and for class "
                                + cfg.CLASS[EXP]
                                + ", year "
                                + str(year)
                            )
                            for SM_day in SM_date_retrieve:
                                showProgress(task, n_tasks)
                                task += 1
                                for SM_time in SM_time_retrieve:
                                    SM_array = mv.read(
                                        data=SM_data, date=SM_day, time=SM_time
                                    )
                                    mv.write(
                                        grib_dir + "SM_" + SM_day + SM_time + ".grib",
                                        SM_array,
                                    )

                        else:
                            print(
                                "\nExtracting SM data for experiment "
                                + cfg.EXPVER[EXP]
                                + " and for class "
                                + cfg.CLASS[EXP]
                                + " from ecfs, year "
                                + str(year)
                            )
                            task = 0
                            for dates in SM_date_retrieve:

                                showProgress(task, n_tasks)
                                task += 1
                                SM_array=mv.ecfs(ecfs_domain = "ec:", file_name = "/RDX/prepIFS/"+cfg.EXPVER[EXP]+"/surface/control/output/swvl1_Param39_"+dates+"_"+dates+"_an.grb")

                                for layer in ("swvl2_Param40_","swvl3_Param41_","swvl4_Param42_"):
                                  SM_array.append(mv.ecfs(ecfs_domain = "ec:", file_name = "/RDX/prepIFS/"+cfg.EXPVER[EXP]+"/surface/control/output/"+layer+dates+"_"+dates+"_an.grb"))
    
                                for SM_time in SM_time_retrieve:
                                    SM_layer_time=mv.read(data=SM_array,time=SM_time)
                                    mv.write(
                                        grib_dir + "SM_" + dates + SM_time + ".grib",
                                        SM_layer_time,
                                    )                                    

                if cfg.extract_ST:

                    if len(ST_date_retrieve) > 0:

                        if (cfg.GRID[EXP]=="av") or (int(cfg.GRID[EXP][1:5])<2559):
                            n_tasks = len(ST_date_retrieve)
                            task = 0
                            print(
                                "\nExtracting ST data for experiment "
                                + cfg.EXPVER[EXP]
                                + " and for class "
                                + cfg.CLASS[EXP]
                                + ", year "
                                + str(year)
                            )
                            ST_data = mv.retrieve(
                                class_=cfg.CLASS[EXP],
                                type=cfg.TYPE[EXP],
                                expver=cfg.EXPVER[EXP],
                                stream=cfg.STREAM[EXP],
                                anoffset=cfg.ANOFFSET[EXP],
                                levtype=cfg.LEVTYPE[EXP],
                                param=["139.128", "170.128", "183.128", "236.128"],
                                date=ST_date_retrieve,
                                time=ST_time_retrieve,
                                step=cfg.STEP[EXP],
                                grid=cfg.GRID[EXP],
                                domain=cfg.DOMAIN[EXP],
                            )

                            for ST_day in ST_date_retrieve:
                                showProgress(task, n_tasks)
                                task += 1
                                for ST_time in ST_time_retrieve:
                                    ST_array = mv.read(
                                        data=ST_data, date=ST_day, time=ST_time
                                    )
                                    mv.write(
                                        grib_dir + "ST_" + ST_day + ST_time + ".grib",
                                        ST_array,
                                    )

                        else:
                            print(
                                "\nExtracting ST data for experiment "
                                + cfg.EXPVER[EXP]
                                + " and for class "
                                + cfg.CLASS[EXP]
                                + " from ecfs, year "
                                + str(year)
                            )
                            task = 0
                            n_tasks = len(ST_date_retrieve)

                            for dates in ST_date_retrieve:
                                showProgress(task, n_tasks)
                                task += 1
                                ST_array=mv.ecfs(ecfs_domain = "ec:", file_name = "/RDX/prepIFS/"+cfg.EXPVER[EXP]+"/surface/control/output/stl1_Param139_"+dates+"_"+dates+"_an.grb")

                                for layer in ("stl2_Param170_","stl3_Param183_","stl4_Param236_"):
                                  ST_array.append(mv.ecfs(ecfs_domain = "ec:", file_name = "/RDX/prepIFS/"+cfg.EXPVER[EXP]+"/surface/control/output/"+layer+dates+"_"+dates+"_an.grb"))
    
                                for ST_time in ST_time_retrieve:
                                    ST_layer_time=mv.read(data=ST_array,time=ST_time)
                                    mv.write(
                                        grib_dir + "ST_" + dates + ST_time + ".grib",
                                        ST_layer_time,
                                    )   ###################################################                                 
    """

    # Preprocess data:
    if cfg.pre_process_SH or cfg.pre_process_LH:

        # Data times

        for n, network in enumerate(cfg.Network): #for every network n, network

            for e, EXP in enumerate(cfg.EXPVER): #for every experiment e, EXP
                print("... pre-processing of experiment [index, name]: " + str(e) + ", " + str(EXP))
                Time_freq = cfg.TIME[e].split("/") #takes TIME from namelist -> Time_freq is list of selected times (hour)
                if len(Time_freq) > 1:
                    time_interval = str(int(Time_freq[1][:2]) - int(Time_freq[0][:2]))
                    data_range = pd.date_range(
                        cfg.analysis_period[0] + "-" + str(Time_freq[0]),
                        cfg.analysis_period[1] + "-" + str(Time_freq[-1]),
                        freq=time_interval + "H",
                    )
                elif len(Time_freq) == 1:
                    data_range = pd.date_range(
                        cfg.analysis_period[0] + "-" + str(Time_freq[0]),
                        cfg.analysis_period[1] + "-" + str(Time_freq[0]),
                    )

                preprocessed_in_situ_dir = cfg.in_situ_dir + "/" + network + "/"
                preprocessed_dir = (
                    cfg.output_dir
                    + "/analysis_data/pre_processed_data/"
                    + network
                    + "/"
                    + EXP
                    + "_"
                    + cfg.CLASS[e]
                    + "/"
                )
                ensure_dir(preprocessed_dir)
                grib_dir = (
                    cfg.output_dir
                    + "/analysis_data/raw_data/"
                    + EXP
                    + "_"
                    + cfg.CLASS[e]
                    + "/"
                )

                for yr in years.year:

                    try:
                        file_in_situ = sorted(
                            glob.glob(
                                preprocessed_in_situ_dir
                                + "/"
                                + str(yr)
                                + "/*"
                                + str(yr)
                            ),
                            key=numericalSort,
                        )

                    except:
                        import pdb

                        pdb.set_trace()
                        print(
                            "\nNo in situ data - skipping year "
                            + str(yr)
                            + " for network "
                            + network
                        )
                        continue

                    if len(file_in_situ) < 1:
                        print(
                            "\nNo in situ data - skipping year "
                            + str(yr)
                            + " for network "
                            + network
                        )
                        continue

                    if cfg.pre_process_SH or cfg.pre_process_LH:
                        try:
                            file = sorted( #list of all file names per year, per experiment (the two fluxes are in the same files)
                                glob.glob(grib_dir + str(yr) + "/fluxes*"),
                                key=numericalSort,
                            )
                            #print(file)
                            if len(file) < 1:
                                raise Exception(
                                    "No analysis fluxes files available in"
                                    + grib_dir
                                    + str(yr)
                                    + ". Check they are extracted."
                                )
                        except:
                            raise Exception(
                                "No analysis fluxes files available in "
                                + grib_dir
                                + str(yr)
                                + ". Check they are extracted."
                            )
                    print(
                        "\nPreprocessing data for experiment "
                        + EXP
                        + " and for class "
                        + cfg.CLASS[e]
                        + ", network "
                        + network
                        + ", year "
                        + str(yr)
                    )

                    ensure_dir(preprocessed_dir + "/" + str(yr))
                    ensure_dir(preprocessed_dir + "/" + str(yr))

                    #meta data stations
                    station_info = (
                        cfg.in_situ_dir
                        + "/station_info/"
                        + network
                        + ""
                        + "/station_coords_"
                        + str(yr)
                    )

                    stations = np.loadtxt(station_info, str, comments="%", skiprows=1) 
                    lat = stations[:, 1].astype("float") #list of lat
                    lon = stations[:, 2].astype("float") #list of lon

                    annual_data_range = data_range[data_range.year == yr]

                    if cfg.pre_process_SH or cfg.pre_procrocess_LH:
                        data_series_fluxes = xr.DataArray(
                            pl.empty((annual_data_range.size, len(lat), 2)), #2 instead of 4, since we have sh and lh (just sfc level)
                            coords=[annual_data_range, range(len(lat)), range(2)],
                            dims=["major_axis", "minor_axis", "items"],
                        )

                    n_tasks = annual_data_range.size
                    task = 0

                    nl=0
                    for idf in range(len(file)):
                        f=file[idf]
                        showProgress(task, n_tasks)
                        task += 1

                        #day_fc = d - pd.Timedelta(int(cfg.STEP[e]), "h")
                        #print(n,e,yr,d)
                        if cfg.pre_process_SH or cfg.pre_process_LH:
                            fset1 = mv.read(f)
                            #fset1 = mv.read( #this works as well
                            #    grib_dir
                            #    + str(day_fc.year)
                            #    + "/fluxes_"
                            #    + str(cfg.EXPVER[e]) #fluxes file names contain experiment
                            #    + "_"
                            #    + str(day_fc.year)
                            #    + "%02d" % (day_fc.month)
                                #+ "%02d" % (day_fc.day) #fluxes file names have no day
                                #+ "%02d" % (day_fc.hour) #fluxes file names have no hour
                                #+ ".grib" #fluxes files do not have the extension
                            #)
                            print("\n... read fluxes file: " + f)

                            file_sh=fset1["sshf"]
                            file_lh=fset1["slhf"]
                            print("... selected sshf and slhf from file")
                            coords=[lat,lon] #that's a list of pairs (lat,lon)
                            #print(coords)
                            #print(d)


                            data_series_SH = mv.nearest_gridpoint( #returns numpy array of shape (number of measurements,number of stations)
                                file_sh,lat,lon
                            )
                            data_series_LH = mv.nearest_gridpoint( #returns numpy array of shape (number of measurements,number of stations)
                                file_lh,lat,lon
                            )
                            print("... selected nearest gridpoints")
                            
                            ### experiment-dependent post processing ###
                            if EXP=="0001":
                                print("... post-processing specified for " + str(EXP))
                                hours_avail=["00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23"]
                                #hours_avail: hours for which the model output is provided
                                monlen=[31,29,31,30,31,30,31,31,30,31,30,31] #days per month, this is for 2020 -> needed for cutting
                                #throw overhang away (the files for shorter months are too long and contain already the beginning of the following month)
                                data_series_SH=data_series_SH[:monlen[idf]*26,:] #26 instead of 24 because of the two artificials zeros per day
                                data_series_LH=data_series_LH[:monlen[idf]*26,:]
                                print("after cutting: " + str(np.shape(data_series_SH)))
                                #convert unit and adapt sign convention:
                                data_series_SH=data_series_SH/3600*(-1)
                                data_series_LH=data_series_LH/3600*(-1)
                                #remove zero flux lines which resulted from former deaccumlation
                                idx=(data_series_SH[:,1]==0)
                                data_series_SH=data_series_SH[~idx,:]
                                data_series_LH=data_series_LH[~idx,:]
                                #prepare requested time resolution (the above data_series_XX contains all hours)
                                for t in Time_freq:
                                    if t not in hours_avail:
                                        print("WARNING: You request a time for which the model does not provide any output")
                                time_index_list=[i for i, item in enumerate(hours_avail) if item in Time_freq]
                                #accumulate on 6h for 00 and 12 utc
                                data_series_SH_new=data_series_SH[0::12]
                                data_series_SH_new[0::2]=1/6*(data_series_SH[19::24,:]+data_series_SH[20::24,:]+data_series_SH[21::24,:]+data_series_SH[22::24,:]+data_series_SH[23::24,:]+data_series_SH[0::24,:])
                                data_series_SH_new[1::2]=1/6*(data_series_SH[7::24,:]+data_series_SH[8::24,:]+data_series_SH[9::24,:]+data_series_SH[10::24,:]+data_series_SH[11::24,:]+data_series_SH[12::24,:])
                                data_series_SH=data_series_SH_new
                                data_series_LH_new=data_series_LH[0::12,:]
                                data_series_LH_new[0::2]=1/6*(data_series_LH[19::24,:]+data_series_LH[20::24,:]+data_series_LH[21::24,:]+data_series_LH[22::24,:]+data_series_LH[23::24,:]+data_series_LH[0::24,:])
                                data_series_LH_new[1::2]=1/6*(data_series_LH[7::24,:]+data_series_LH[8::24,:]+data_series_LH[9::24,:]+data_series_LH[10::24,:]+data_series_LH[11::24,:]+data_series_LH[12::24,:])
                                data_series_LH=data_series_LH_new                                
                                
                            elif EXP=="hyfs":
                                print("... post-processing specified for " + str(EXP))
                                hours_avail=["00","06","12","18"] #hours_avail: hours for which the model output is provided
                                #convert unit and adapt sign convention:
                                data_series_SH=data_series_SH/(3600*6)*-1
                                data_series_LH=data_series_LH/(3600*6)*-1
                                #deaccumulation
                                data_series_SH[3::4]=data_series_SH[3::4]-data_series_SH[2::4]
                                data_series_SH[2::4]=data_series_SH[2::4]-data_series_SH[1::4]
                                data_series_SH[1::4]=data_series_SH[1::4]-data_series_SH[0::4]
                                data_series_LH[3::4]=data_series_LH[3::4]-data_series_LH[2::4]
                                data_series_LH[2::4]=data_series_LH[2::4]-data_series_LH[1::4]
                                data_series_LH[1::4]=data_series_LH[1::4]-data_series_LH[0::4]
                                #prepare requested time resolution (the above data_series_XX contains all hours)
                                for t in Time_freq:
                                    if t not in hours_avail:
                                        print("WARNING: You request a time for which the model does not provide any output")
                                time_index_list=[i for i, item in enumerate(hours_avail) if item in Time_freq]
                                data_series_SH=data_series_SH[0::2,:] #quick fix
                                data_series_LH=data_series_LH[0::2,:]

                            else:
                                print("WARNING: you use an experiment for which no special post-processing is defined")


                            #concatenate of all files
                            nl2=np.shape(data_series_SH)[0]
                            data_series_fluxes[nl:(nl+nl2),:,0]=data_series_SH
                            data_series_fluxes[nl:(nl+nl2),:,1]=data_series_LH
                            nl=nl+nl2



                    ### write in new output file .dat ###
                    for lat_lon in range(len(lat)):
                        if cfg.pre_process_SH or cfg.pre_process_LH:
                            filename = (
                                preprocessed_dir
                                + "/"
                                + str(yr)
                                + "/"
                                + "fluxes"
                                + "_"
                                + "%.4f" % float(lat[lat_lon])
                                + ","
                                + "%.4f" % float(lon[lat_lon])
                                + ".dat"
                            )
                            data_series_fluxes.loc[:, lat_lon, :].to_pandas().round(
                                6
                            ).to_csv(
                                filename,
                                sep=" ",
                                header=False,
                                date_format="%Y%m%d%H",
                            )
                            


