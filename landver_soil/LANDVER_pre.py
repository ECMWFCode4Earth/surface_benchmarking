#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
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

    if cfg.extract_SM or cfg.extract_ST:

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
                                    )                                    


    # Preprocess data:
    if cfg.pre_process_SM or cfg.pre_process_ST:

        # Data times

        for n, network in enumerate(cfg.Network):

            for e, EXP in enumerate(cfg.EXPVER):
                Time_freq = cfg.TIME[e].split("/")
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

                preprocessed_in_situ_dir = cfg.in_situ_dir + "/" + network + "_ISMN/"
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

                    if cfg.pre_process_SM:
                        try:
                            file = sorted(
                                glob.glob(grib_dir + str(yr) + "/SM*"),
                                key=numericalSort,
                            )
                            if len(file) < 1:
                                raise Exception(
                                    "No analysis soil moisture files available in"
                                    + grib_dir
                                    + str(yr)
                                    + ". Check they are extracted."
                                )
                        except:
                            raise Exception(
                                "No analysis soil moisture files available in "
                                + grib_dir
                                + str(yr)
                                + ". Check they are extracted."
                            )

                    if cfg.pre_process_ST:
                        try:
                            file_ST = sorted(
                                glob.glob(grib_dir + str(yr) + "/ST*"),
                                key=numericalSort,
                            )
                            if len(file_ST) < 1:
                                raise Exception(
                                    "No analysis soil temperature files available in"
                                    + grib_dir
                                    + str(yr)
                                    + ". Check they are extracted."
                                )
                        except:
                            raise Exception(
                                "No analysis soil temperature files available in "
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

                    station_info = (
                        cfg.in_situ_dir
                        + "/station_info/"
                        + network
                        + "_ISMN"
                        + "/station_coords_"
                        + str(yr)
                    )

                    stations = np.loadtxt(station_info, str, comments="%", skiprows=1)
                    lat = stations[:, 1].astype("float")
                    lon = stations[:, 2].astype("float")

                    annual_data_range = data_range[data_range.year == yr]

                    if cfg.pre_process_SM:
                        data_series_SM = xr.DataArray(
                            pl.empty((annual_data_range.size, len(lat), 4)),
                            coords=[annual_data_range, range(len(lat)), range(4)],
                            dims=["major_axis", "minor_axis", "items"],
                        )

                    if cfg.pre_process_ST:
                        data_series_ST = xr.DataArray(
                            pl.empty((annual_data_range.size, len(lat), 4)),
                            coords=[annual_data_range, range(len(lat)), range(4)],
                            dims=["major_axis", "minor_axis", "items"],
                        )

                    n_tasks = annual_data_range.size
                    task = 0

                    for d in annual_data_range:
                        showProgress(task, n_tasks)
                        task += 1

                        day_fc = d - pd.Timedelta(int(cfg.STEP[e]), "h")

                        if cfg.pre_process_SM:
                            fset1 = mv.read(
                                grib_dir
                                + str(day_fc.year)
                                + "/SM_"
                                + str(day_fc.year)
                                + "%02d" % (day_fc.month)
                                + "%02d" % (day_fc.day)
                                + "%02d" % (day_fc.hour)
                                + ".grib"
                            )

                            data_series_SM.loc[d, :, :] = mv.nearest_gridpoint(
                                fset1, lat, lon
                            ).transpose()

                            data_series_SM.loc[d, :, :].fillna(1.7e38)

                        if cfg.pre_process_ST:
                            fset2 = mv.read(
                                grib_dir
                                + str(day_fc.year)
                                + "/ST_"
                                + str(day_fc.year)
                                + "%02d" % (day_fc.month)
                                + "%02d" % (day_fc.day)
                                + "%02d" % (day_fc.hour)
                                + ".grib"
                            )
                            data_series_ST.loc[d, :, :] = mv.nearest_gridpoint(
                                fset2, lat, lon
                            ).transpose()

                            data_series_ST.loc[d, :, :].fillna(1.7e38)
                            # Apply lapse rate correction to soil temperature
                            # for stations with large height difference # model-obs

                    for lat_lon in range(len(lat)):
                        if cfg.pre_process_SM:
                            filename = (
                                preprocessed_dir
                                + "/"
                                + str(yr)
                                + "/SM_"
                                + "%.4f" % float(lat[lat_lon])
                                + ","
                                + "%.4f" % float(lon[lat_lon])
                                + ".dat"
                            )

                            data_series_SM.loc[:, lat_lon, :].to_pandas().round(
                                6
                            ).to_csv(
                                filename, sep=" ", header=False, date_format="%Y%m%d%H"
                            )

                        if cfg.pre_process_ST:
                            filename_ST = (
                                preprocessed_dir
                                + "/"
                                + str(yr)
                                + "/ST_"
                                + "%.4f" % float(lat[lat_lon])
                                + ","
                                + "%.4f" % float(lon[lat_lon])
                                + ".dat"
                            )
                            data_series_ST.loc[:, lat_lon, :].to_pandas().round(
                                6
                            ).to_csv(
                                filename_ST,
                                sep=" ",
                                header=False,
                                date_format="%Y%m%d%H",
                            )
