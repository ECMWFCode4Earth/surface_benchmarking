#!/usr/bin/env python3
"""
Script to perform validation of two or more soil moisture or soil temperature analyses.
Last modified by David Fairbairn on 5/4/22. 
"""
import numpy as np
import glob
import pylab as pl
import pandas as pd
import re
import copy
from Common_functions import showProgress, ensure_dir, clean_dir
from scipy.stats import pearsonr
from warnings import filterwarnings
import metview as mv

filterwarnings("ignore")
import Validation_functions as ldasv

numbers = re.compile(r"(\d+)")


def in_situ_validation(cfg, var, times, land_type, df):

    # Initialize lat/lon/R dicts for station maps
    lat_networks = dict()
    lon_networks = dict()
    lat_networks_spring = dict()
    lon_networks_spring = dict()
    R_networks = dict()
    Bias_networks = dict()
    n_stations = dict()
    Ano_R_compare = dict()
    R = dict()
    Ano_R = dict()
    RMSD = dict()
    Ub_RMSD = dict()
    ME = dict()
    P_value = dict()
    ME = dict()
    P_value = dict()
    Ano_P_value = dict()
    w = dict()
    STD_an = dict()
    STD_obs = dict()
    Ano_R_conf_upper = dict()
    Ano_R_conf_lower = dict()

    if df is None:
        df = pd.DataFrame(
            data=None,
            columns=[
                "expver",
                "date",
                "type",
                "network",
                "variable",
                "metric",
                "link",
            ],
        )

    for e, EXP in enumerate([m + "_" + n for m, n in zip(cfg.EXPVER, cfg.CLASS)]):
        # Data times
        try:
            Time_freq = cfg.TIME[e].split("/")

            if times == "00_UTC":
                Time_freq = ["00"]
            elif times == "12_UTC":
                Time_freq = ["12"]

            if len(Time_freq) > 1:
                time_interval = str(int(Time_freq[1][:2]) - int(Time_freq[0][:2]))
                freq = time_interval + "H"
                data_range = pd.date_range(
                    cfg.analysis_period[0] + "-" + str(Time_freq[0]),
                    cfg.analysis_period[1] + "-" + str(Time_freq[-1]),
                    freq=freq,
                )
            elif len(Time_freq) == 1:
                data_range = pd.date_range(
                    cfg.analysis_period[0] + "-" + str(Time_freq[0]),
                    cfg.analysis_period[1] + "-" + str(Time_freq[0]),
                )
                time_interval = "24"
                freq = "24H"
        except:
            raise Exception(
                "Incorrect date in analysis period. Check Options namelist. "
            )

        years = pd.date_range(
            str(data_range[0].year), str(data_range[-1].year + 1), freq="A"
        )
        if land_type != "all_land":
            land_type_codes = ldasv.land_class_lookup_table(land_type)

        # For seasonal scores:
        mlist = [[12, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]]
        slist = ["winter", "spring", "summer", "autumn"]
        ylist = list()
        for yr in years.year:
            ylist.append(str(yr))
        sdict = {k: v for v, ks in zip(slist, mlist) for k in ks}

        start_period = (
            str(data_range[0].year)
            + "-%02d" % (data_range[0].month)
            + "-%02d" % (data_range[0].day)
        )
        end_period = (
            str(data_range[-1].year)
            + "-%02d" % (data_range[-1].month)
            + "-%02d" % (data_range[-1].day)
        )

        # Create the relevant post-processing directories
        Report_dir = (
            cfg.output_dir
            + "/"
            + var
            + "_validation_"
            + EXP
            + "_"
            + start_period
            + "_"
            + end_period
        )
        ensure_dir(Report_dir)
        tables_dir = Report_dir + "/tables/"
        ensure_dir(tables_dir)
        plots_dir = Report_dir + "/plots/"
        ensure_dir(plots_dir)
        time_series_dir = Report_dir + "/plots/time_series/"
        ensure_dir(time_series_dir)

        # Initialize lat/lon/R dicts for station maps
        lat_networks[EXP] = dict()
        lon_networks[EXP] = dict()
        lat_networks_spring[EXP] = dict()
        lon_networks_spring[EXP] = dict()
        R_networks[EXP] = dict()
        Bias_networks[EXP] = dict()
        Ano_R_compare[EXP] = dict()
        n_stations[EXP] = dict()
        R[EXP] = dict()
        Ano_R[EXP] = dict()
        RMSD[EXP] = dict()
        Ub_RMSD[EXP] = dict()
        ME[EXP] = dict()
        P_value[EXP] = dict()
        ME[EXP] = dict()
        Ano_R_conf_upper[EXP] = dict()
        Ano_R_conf_lower[EXP] = dict()
        P_value[EXP] = dict()
        Ano_P_value[EXP] = dict()
        w[EXP] = dict()
        STD_an[EXP] = dict()
        STD_obs[EXP] = dict()
        Data_available = dict()
        Analysis_available = dict()
        ST_data_available = dict()
        file = dict()
        file_ST = dict()
        file_in_situ = dict()

        # Check in situ data and analysis data are present:
        for n, network in enumerate(cfg.Network):
            file[network] = dict()
            file_ST[network] = dict()
            file_in_situ[network] = dict()
            frozen_conditions = dict()
            Data_available[network] = dict()
            ST_data_available[network] = dict()
            Analysis_available[network] = dict()
            preprocessed_analysis_dir = (
                cfg.output_dir
                + "/analysis_data/pre_processed_data/"
                + network
                + "/"
                + EXP
            )

            time_series_dir_net = time_series_dir + network
            ensure_dir(time_series_dir_net)

            preprocessed_in_situ_dir = cfg.in_situ_dir + "/" + network + "_ISMN/"
            depths, analysis_depths = ldasv.network_SM_depths(network)

            for yr in years.year:
                try:
                    file[network][str(yr)] = sorted(
                        glob.glob(
                            preprocessed_analysis_dir + "/" + str(yr) + "/" + var + "*"
                        ),
                        key=ldasv.numericalSort,
                    )
                    if len(file[network][str(yr)]) < 1:
                        Analysis_available[network][str(yr)] = False
                    else:
                        Analysis_available[network][str(yr)] = True
                except:
                    Data_available[network][str(yr)] = False
                try:
                    file_in_situ[network][str(yr)] = sorted(
                        glob.glob(
                            preprocessed_in_situ_dir
                            + "/"
                            + str(yr)
                            + "/"
                            + network
                            + "_"
                            + "%.2f" % (depths[0])
                            + "*"
                        ),
                        key=ldasv.numericalSort,
                    )

                    if len(file_in_situ[network][str(yr)]) < 1:
                        Data_available[network][str(yr)] = False
                    else:
                        Data_available[network][str(yr)] = True
                except:
                    Data_available[network][str(yr)] = False
        if var == "SM":
            # Check if analysis temperature files are available for filtering frozen soil moisture:
            if cfg.ST_quality_control:
                try:
                    file_ST[network][str(yr)] = sorted(
                        glob.glob(preprocessed_analysis_dir + "/" + str(yr) + "/ST*"),
                        key=ldasv.numericalSort,
                    )
                    if len(file_ST[network][str(yr)]) < 1:
                        ST_data_available[network][str(yr)] = False
                    else:
                        ST_data_available[network][str(yr)] = True
                except:
                    ST_data_available[network][str(yr)] = False

        # Loop over networks:
        for n, network in enumerate(cfg.Network):

            # Define directories:
            preprocessed_analysis_dir = (
                cfg.output_dir
                + "/analysis_data/pre_processed_data/"
                + network
                + "/"
                + EXP
            )
            preprocessed_in_situ_dir = cfg.in_situ_dir + "/" + network + "_ISMN/"
            depths, analysis_depths = ldasv.network_SM_depths(network)

            # Loop over layers:
            for l, layer in enumerate(cfg.validation_layer):

                time_series_dir_net_layer = time_series_dir + network + "/" + layer
                ensure_dir(time_series_dir_net_layer)
                time_series_dir_net_layer_times = (
                    time_series_dir + network + "/" + layer + "/" + times
                )
                ensure_dir(time_series_dir_net_layer_times)
                time_series_dir_net_layer_times_land = (
                    time_series_dir
                    + network
                    + "/"
                    + layer
                    + "/"
                    + times
                    + "/"
                    + land_type
                )
                ensure_dir(time_series_dir_net_layer_times_land)
                clean_dir(time_series_dir_net_layer_times_land)

                # Initialize dictionaries for each network/layer/period:
                R[EXP][network + "_" + layer] = dict()
                Ano_R[EXP][network + "_" + layer] = dict()
                RMSD[EXP][network + "_" + layer] = dict()
                Ub_RMSD[EXP][network + "_" + layer] = dict()
                Ano_R_conf_upper[EXP][network + "_" + layer] = dict()
                Ano_R_conf_lower[EXP][network + "_" + layer] = dict()
                ME[EXP][network + "_" + layer] = dict()
                P_value[EXP][network + "_" + layer] = dict()
                Ano_P_value[EXP][network + "_" + layer] = dict()
                w[EXP][network + "_" + layer] = dict()
                STD_an[EXP][network + "_" + layer] = dict()
                STD_obs[EXP][network + "_" + layer] = dict()
                lat_networks[EXP][network + "_" + layer] = dict()
                lon_networks[EXP][network + "_" + layer] = dict()
                lat_networks_spring[EXP][network + "_" + layer] = dict()
                lon_networks_spring[EXP][network + "_" + layer] = dict()
                R_networks[EXP][network + "_" + layer] = dict()
                Bias_networks[EXP][network + "_" + layer] = dict()
                Ano_R_compare[EXP][network + "_" + layer] = dict()

                R[EXP][network + "_" + layer]["period"] = list()
                Ano_R[EXP][network + "_" + layer]["period"] = list()
                RMSD[EXP][network + "_" + layer]["period"] = list()
                Ano_R_conf_upper[EXP][network + "_" + layer]["period"] = list()
                Ano_R_conf_lower[EXP][network + "_" + layer]["period"] = list()
                Ub_RMSD[EXP][network + "_" + layer]["period"] = list()
                ME[EXP][network + "_" + layer]["period"] = list()
                P_value[EXP][network + "_" + layer]["period"] = list()
                Ano_P_value[EXP][network + "_" + layer]["period"] = list()
                w[EXP][network + "_" + layer]["period"] = list()
                STD_an[EXP][network + "_" + layer]["period"] = list()
                STD_obs[EXP][network + "_" + layer]["period"] = list()

                for season in slist:
                    R[EXP][network + "_" + layer][season] = list()
                    Ano_R[EXP][network + "_" + layer][season] = list()
                    RMSD[EXP][network + "_" + layer][season] = list()
                    Ano_R_conf_upper[EXP][network + "_" + layer][season] = list()
                    Ano_R_conf_lower[EXP][network + "_" + layer][season] = list()
                    Ub_RMSD[EXP][network + "_" + layer][season] = list()
                    ME[EXP][network + "_" + layer][season] = list()
                    P_value[EXP][network + "_" + layer][season] = list()
                    Ano_P_value[EXP][network + "_" + layer][season] = list()
                    w[EXP][network + "_" + layer][season] = list()
                    STD_an[EXP][network + "_" + layer][season] = list()
                    STD_obs[EXP][network + "_" + layer][season] = list()

                for yr in years.year:
                    R[EXP][network + "_" + layer][str(yr)] = list()
                    Ano_R[EXP][network + "_" + layer][str(yr)] = list()
                    RMSD[EXP][network + "_" + layer][str(yr)] = list()
                    Ano_R_conf_upper[EXP][network + "_" + layer][str(yr)] = list()
                    Ano_R_conf_lower[EXP][network + "_" + layer][str(yr)] = list()
                    Ub_RMSD[EXP][network + "_" + layer][str(yr)] = list()
                    ME[EXP][network + "_" + layer][str(yr)] = list()
                    P_value[EXP][network + "_" + layer][str(yr)] = list()
                    Ano_P_value[EXP][network + "_" + layer][str(yr)] = list()
                    w[EXP][network + "_" + layer][str(yr)] = list()
                    STD_an[EXP][network + "_" + layer][str(yr)] = list()
                    STD_obs[EXP][network + "_" + layer][str(yr)] = list()

                # Initialize latitude and longitude lists for storing point locations:
                latitude = list()
                longitude = list()

                print(
                    "\n###############################################################################"
                )
                print(
                    "\nValidating "
                    + var
                    + " data for experiment "
                    + EXP[0:4]
                    + " and for class "
                    + cfg.CLASS[e]
                    + ", network "
                    + network
                    + ", layer "
                    + layer
                    + ", "
                    + times
                    + " and "
                    + land_type
                )
                print(
                    "\n###############################################################################"
                )

                # Filter out in situ point locations where both analysis and in situ files are available:
                lat_lon_list = dict()

                # Define lat/lon dicts
                for yr in years.year:
                    frozen_conditions[str(yr)] = list()
                    lat_lon_list[str(yr)] = dict()
                    lat_lon_list[str(yr)]["analysis"] = list()
                    lat_lon_list[str(yr)]["insitu"] = list()
                    lat_lon_list["merged"] = dict()
                    lat_lon_list["merged"]["analysis"] = list()
                    lat_lon_list["merged"]["insitu"] = list()
                    lat_lon_list["merged"]["combined_index"] = list()
                    lat_lon_list["uniform_prec"] = dict()
                    lat_lon_list["uniform_prec"]["analysis"] = list()
                    lat_lon_list["uniform_prec"]["insitu"] = list()

                    station_info = (
                        cfg.in_situ_dir
                        + "/station_info/"
                        + network
                        + "_ISMN"
                        + "/station_coords_"
                        + str(yr)
                    )

                    stations = np.loadtxt(station_info, str, comments="%", skiprows=1)

                    for i in range(len(file[network][str(yr)])):
                        lat, lon = ldasv.find_lat_lon(
                            file[network][str(yr)][i], "analysis"
                        )

                        if land_type == "all_land":
                            lat_lon_list[str(yr)]["analysis"].append(
                                "%.4f" % float(lat) + "_" + "%.4f" % float(lon)
                            )

                        elif land_type != "all_land":

                            station_index = np.where(
                                (stations[:, 1].astype(float) == float(lat))
                                & (stations[:, 2].astype(float) == float(lon))
                            )

                            if station_index[0].size > 0:

                                if (
                                    int(stations[station_index[0][0]][4])
                                    in land_type_codes
                                ):
                                    lat_lon_list[str(yr)]["analysis"].append(
                                        "%.4f" % float(lat) + "_" + "%.4f" % float(lon)
                                    )

                    for i in range(len(file_in_situ[network][str(yr)])):
                        lat, lon = ldasv.find_lat_lon(
                            file_in_situ[network][str(yr)][i], "insitu"
                        )

                        if land_type == "all_land":
                            lat_lon_list[str(yr)]["insitu"].append(
                                "%.4f" % float(lat) + "_" + "%.4f" % float(lon)
                            )

                        elif land_type != "all_land":

                            station_index = np.where(
                                (stations[:, 1].astype(float) == float(lat))
                                & (stations[:, 2].astype(float) == float(lon))
                            )

                            if station_index[0].size > 0:

                                if (
                                    int(stations[station_index[0][0]][4])
                                    in land_type_codes
                                ):
                                    lat_lon_list[str(yr)]["insitu"].append(
                                        "%.4f" % float(lat) + "_" + "%.4f" % float(lon)
                                    )
                year_count = 0

                for yr in range(len(years.year)):

                    if Data_available[network][str(years.year[yr])] == False:
                        print(
                            "\nNo in situ data - skipping year "
                            + str(years.year[yr])
                            + " for network "
                            + network
                        )
                        continue

                    if Analysis_available[network][str(years.year[yr])] == False:
                        print(
                            "\nNo analysis data - skipping year "
                            + str(years.year[yr])
                            + " for network "
                            + network
                        )
                        continue

                    if year_count == 0:
                        lat_lon_list["merged"]["analysis"] = lat_lon_list[
                            str(years.year[yr])
                        ]["analysis"]
                        lat_lon_list["merged"]["insitu"] = lat_lon_list[
                            str(years.year[yr])
                        ]["insitu"]
                    else:

                        # Accumulate all available stations over period:
                        if len(lat_lon_list["merged"]["analysis"]) > len(
                            lat_lon_list[str(years.year[yr])]["analysis"]
                        ):
                            new_stations_an = list(
                                set(lat_lon_list["merged"]["analysis"]).difference(
                                    lat_lon_list[str(years.year[yr])]["analysis"]
                                )
                            )
                        else:
                            new_stations_an = list(
                                set(
                                    lat_lon_list[str(years.year[yr])]["analysis"]
                                ).difference(lat_lon_list["merged"]["analysis"])
                            )

                        if len(lat_lon_list["merged"]["insitu"]) > len(
                            lat_lon_list[str(years.year[yr])]["insitu"]
                        ):
                            new_stations_in_situ = list(
                                set(lat_lon_list["merged"]["insitu"]).difference(
                                    lat_lon_list[str(years.year[yr])]["insitu"]
                                )
                            )
                        else:
                            new_stations_in_situ = list(
                                set(
                                    lat_lon_list[str(years.year[yr])]["insitu"]
                                ).difference(lat_lon_list["merged"]["insitu"])
                            )
                        if len(new_stations_an) > 0:
                            lat_lon_list["merged"]["analysis"].extend(new_stations_an)
                        if len(new_stations_in_situ) > 0:
                            lat_lon_list["merged"]["insitu"].extend(
                                new_stations_in_situ
                            )

                year_count = year_count + 1

                # Find lat/lon coordinates for common stations:
                for i in range(len(lat_lon_list["merged"]["analysis"])):
                    lat_analysis = re.split(
                        "[_,]", lat_lon_list["merged"]["analysis"][i]
                    )[0]
                    lon_analysis = re.split(
                        "[_,]", lat_lon_list["merged"]["analysis"][i]
                    )[1]
                    lat_lon_list["uniform_prec"]["analysis"].append(
                        "%.2f" % float(lat_analysis)
                        + "_"
                        + "%.2f" % float(lon_analysis)
                    )

                for i in range(len(lat_lon_list["merged"]["insitu"])):
                    lat_insitu = re.split("[_,]", lat_lon_list["merged"]["insitu"][i])[
                        0
                    ]
                    lon_insitu = re.split("[_,]", lat_lon_list["merged"]["insitu"][i])[
                        1
                    ]
                    lat_lon_list["uniform_prec"]["insitu"].append(
                        "%.2f" % float(lat_insitu) + "_" + "%.2f" % float(lon_insitu)
                    )

                lat_lon_list["merged"]["analysis_index"] = np.nonzero(
                    np.in1d(
                        lat_lon_list["uniform_prec"]["analysis"],
                        lat_lon_list["uniform_prec"]["insitu"],
                    )
                )[0].tolist()
                lat_lon_list["merged"]["insitu_index"] = np.nonzero(
                    np.in1d(
                        lat_lon_list["uniform_prec"]["insitu"],
                        lat_lon_list["uniform_prec"]["analysis"],
                    )
                )[0].tolist()

                n_tasks = len(lat_lon_list["merged"]["analysis_index"])
                task = 0

                # Loop over the stations:
                for i in range(len(lat_lon_list["merged"]["analysis_index"])):

                    if len(lat_lon_list["merged"]["analysis_index"]) == 0:
                        print("No " + land_type + " available for network " + network)
                        break

                    showProgress(task, n_tasks)
                    task += 1

                    lat_analysis = re.split(
                        "[_,]",
                        lat_lon_list["merged"]["analysis"][
                            lat_lon_list["merged"]["analysis_index"][i]
                        ],
                    )[0]
                    lon_analysis = re.split(
                        "[_,]",
                        lat_lon_list["merged"]["analysis"][
                            lat_lon_list["merged"]["analysis_index"][i]
                        ],
                    )[1]
                    lat_insitu = lat_analysis  # re.split("[_,]", lat_lon_list['merged']['insitu'][lat_lon_list['merged']['insitu_index'][i]])[0]
                    lon_insitu = lon_analysis  # re.split("[_,]", lat_lon_list['merged']['insitu'][lat_lon_list['merged']['insitu_index'][i]])[1]

                    latitude.append(lat_analysis)
                    longitude.append(lon_analysis)

                    insitu_df = pd.DataFrame(pl.empty((data_range.size)), data_range)
                    analysis_df = pd.DataFrame(pl.empty((data_range.size)), data_range)

                    if layer == "Surface":
                        insitu_depths = [depths[0]]
                        an_depths = [analysis_depths[0]]
                    elif layer == "Rootzone":
                        insitu_depths = depths
                        an_depths = analysis_depths

                    if (var == "SM") and (cfg.ST_quality_control):
                        ST_insitu_df = pd.DataFrame(
                            pl.empty((data_range.size)), data_range
                        )
                        ST_df = pd.DataFrame(pl.empty((data_range.size)), data_range)

                    (
                        insitu_df,
                        ST_insitu_df,
                        max_insitu,
                        min_insitu,
                    ) = ldasv.read_and_rescale_insitu_data(
                        var,
                        preprocessed_in_situ_dir,
                        years.year,
                        network,
                        lat_insitu,
                        lon_insitu,
                        insitu_depths,
                        data_range,
                        cfg.ST_quality_control,
                        cfg.daily_obs_time_average,
                        cfg.Rescale_data,
                        cfg.ST_QC_threshold - 273.15,
                    )

                    analysis_df, ST_df = ldasv.read_and_rescale_analysis_data(
                        var,
                        preprocessed_analysis_dir,
                        i,
                        years.year,
                        lat_analysis,
                        lon_analysis,
                        an_depths,
                        data_range,
                        freq,
                        cfg.ST_quality_control,
                        cfg.Rescale_data,
                        cfg.ST_ML10[e],
                        cfg.ST_QC_threshold,
                    )

                    insitu_df.loc[:][pl.isnan(analysis_df.loc[:])] = pl.nan

                    if cfg.filterOro:

                        for yr_count in years.year:
                            station_info = (
                                cfg.in_situ_dir
                                + "/station_info/"
                                + network
                                + "_ISMN"
                                + "/station_coords_"
                                + str(yr_count)
                            )
                            try:
                                stations = np.loadtxt(
                                    station_info, str, comments="%", skiprows=1
                                )
                            except:
                                continue

                            try:
                                st_index = np.where(
                                    (stations[:, 1] == lat_analysis)
                                    & (stations[:, 2] == lon_analysis)
                                )[0][0]
                            except:
                                continue

                            height = stations[st_index, 3]

                            if (cfg.GRID[e][0] == "o") or (cfg.GRID[e][0] == "O"):
                                grid_type = "_4"
                                clim_file = (
                                    "orog_"
                                    + str(np.int(cfg.GRID[e][1:]) + 1)
                                    + grid_type
                                )
                                foro = mv.read(cfg.clim_dir + "/" + clim_file + ".grb")
                            if (cfg.GRID[e][0] == "n") or (cfg.GRID[e][0] == "N"):
                                grid_type = "l_2"
                                clim_file = (
                                    "orog_"
                                    + str((np.int(cfg.GRID[e][1:]) * 2) - 1)
                                    + grid_type
                                )
                                foro = mv.read(cfg.clim_dir + "/" + clim_file + ".grb")
                            else:
                                foro = mv.read(
                                    cfg.clim_dir + "/" + cfg.clim_file + ".grb"
                                )

                            data_analysis_Oro = mv.nearest_gridpoint(
                                foro, float(lat_analysis), float(lon_analysis)
                            )
                            if (
                                np.abs(data_analysis_Oro - float(height))
                                > cfg.Oro_diff_threshold
                            ):
                                print(
                                    " Warning: screening where orog_diff>"
                                    + str(cfg.Oro_diff_threshold)
                                    + " for lat="
                                    + lat_analysis
                                    + ", lon="
                                    + lon_analysis
                                )
                                analysis_df.loc[str(yr_count)] = np.nan

                                if cfg.ST_quality_control and (var == "SM"):
                                    ST_df.loc[str(yr_count)] = np.nan

                    # ANOMALY time-series using 5 week moving average
                    insitu_ano = np.copy(insitu_df.values * np.nan)
                    analysis_ano = np.copy(analysis_df.values * np.nan)
                    ap = np.int(
                        24 / np.int(time_interval)
                    )  # first calculate frequency of data per day
                    for yy_len in range(14 * ap, len(insitu_df.loc[:]) - 14 * ap):
                        dd = insitu_df.loc[:].values[
                            yy_len - 15 * ap : yy_len + 15 * ap
                        ][
                            ~np.isnan(
                                insitu_df.loc[:].values[
                                    yy_len - 15 * ap : yy_len + 15 * ap
                                ]
                            )
                        ]
                        zz = analysis_df.values[yy_len - 15 * ap : yy_len + 15 * ap][
                            ~np.isnan(
                                analysis_df.values[yy_len - 15 * ap : yy_len + 15 * ap]
                            )
                        ]
                        if len(dd) > 5 * ap:
                            analysis_ano[yy_len] = (
                                analysis_df.values[yy_len] - np.mean(zz)
                            ) / np.std(zz)
                            insitu_ano[yy_len] = (
                                insitu_df.loc[:].values[yy_len] - np.mean(dd)
                            ) / np.std(dd)

                    # Filter NaN values where applicable
                    array = ldasv.filter_nan(analysis_df.values, insitu_df.values)
                    array_ano = ldasv.filter_nan(analysis_ano, insitu_ano)

                    

                    if (var == "SM") and cfg.ST_quality_control:
                        for yr in years.year:
                            if (
                                len(
                                    ST_insitu_df[str(yr)][
                                        ST_insitu_df[str(yr)]
                                        > (cfg.ST_QC_threshold - 273.15)
                                    ]
                                )
                                == 0
                            ) and (
                                len(ST_insitu_df[str(yr)][ST_insitu_df[str(yr)] > 0.0])
                                > 0
                            ):
                                frozen_conditions[str(yr)].append("0")
                            else:
                                frozen_conditions[str(yr)].append("1")

                    # Calculate stats where a reasonable amount of data is available:
                    if len(array[1]) > cfg.min_obs:

                        # Calculate scores over the entire period
                        bias_val, rmsd_val, corr_val, ano_corr_val = (
                            ldasv.bias(array[0], array[1]),
                            ldasv.rmse(array[0], array[1]),
                            ldasv.correlation(array[0], array[1]),
                            ldasv.correlation(array_ano[0], array_ano[1]),
                        )
                        std_an_val = np.std(array[0])
                        std_obs_val = np.std(array[1])
                        if (
                            (var == "SM")
                            and (cfg.SM_units == "Volumetric")
                            and cfg.Rescale_data
                        ):
                            rmsd_val = rmsd_val * (max_insitu - min_insitu)
                            bias_val = bias_val * (max_insitu - min_insitu)
                            std_an_val = std_an_val * (max_insitu - min_insitu)
                            std_obs_val = std_obs_val * (max_insitu - min_insitu)
                        ub_rmsd_val = np.sqrt(rmsd_val ** 2 - bias_val ** 2)
                        if (array[0].size > 2) and (array_ano[0].size > 2):
                            P_val = pearsonr(array[0], array[1])[1]
                            Ano_P_val = pearsonr(array_ano[0], array_ano[1])[1]
                        else:
                            P_val = pl.nan
                            Ano_P_val = pl.nan

                        if cfg.stat_sig:
                            insitu_df_ano = copy.deepcopy(insitu_df)
                            analysis_df_ano = copy.deepcopy(analysis_df)
                            insitu_df_ano.loc[:] = insitu_ano
                            analysis_df_ano.loc[:] = analysis_ano

                            if (
                                (Ano_P_val < 0.05)
                                and (P_val < 0.05)
                                and (ano_corr_val > -1.0)
                                and (ano_corr_val < 1.0)
                            ):
                                lower_conf, upper_conf = ldasv.Confidence_intervals(
                                    analysis_df_ano,
                                    insitu_df_ano,
                                    ano_corr_val,
                                    np.size(array_ano[0]),
                                )
                                Ano_R_conf_upper[EXP][network + "_" + layer][
                                    "period"
                                ].append(upper_conf)
                                Ano_R_conf_lower[EXP][network + "_" + layer][
                                    "period"
                                ].append(lower_conf)

                                if np.isnan(
                                    lower_conf
                                ):  # Set scores to nan if statistical significance sampling size insufficient:
                                    bias_val, rmsd_val, corr_val, ano_corr_val = (
                                        pl.nan,
                                        pl.nan,
                                        pl.nan,
                                        pl.nan,
                                    )
                                    std_an_val = pl.nan
                                    std_obs_val = pl.nan
                            else:
                                Ano_R_conf_upper[EXP][network + "_" + layer][
                                    "period"
                                ].append(pl.nan)
                                Ano_R_conf_lower[EXP][network + "_" + layer][
                                    "period"
                                ].append(pl.nan)

                        # Update period scores
                        R[EXP][network + "_" + layer]["period"].append(corr_val)
                        Ano_R[EXP][network + "_" + layer]["period"].append(ano_corr_val)
                        RMSD[EXP][network + "_" + layer]["period"].append(rmsd_val)
                        ME[EXP][network + "_" + layer]["period"].append(bias_val)
                        Ub_RMSD[EXP][network + "_" + layer]["period"].append(
                            ub_rmsd_val
                        )
                        P_value[EXP][network + "_" + layer]["period"].append(P_val)
                        Ano_P_value[EXP][network + "_" + layer]["period"].append(
                            Ano_P_val
                        )
                        STD_an[EXP][network + "_" + layer]["period"].append(std_an_val)
                        STD_obs[EXP][network + "_" + layer]["period"].append(
                            std_obs_val
                        )

                        s_index = analysis_df.index.to_series().dt.month.map(sdict)

                        # Calculate seasonal scores
                        for season in slist:

                            analysis_df_s = analysis_df[s_index == season]
                            insitu_df_s = insitu_df[s_index == season]
                            analysis_ano_s = analysis_ano[s_index == season]
                            insitu_ano_s = insitu_ano[s_index == season]

                            array = ldasv.filter_nan(
                                analysis_df_s.values, insitu_df_s.values
                            )
                            array_ano = ldasv.filter_nan(analysis_ano_s, insitu_ano_s)

                            bias_val, rmsd_val, corr_val, ano_corr_val = (
                                ldasv.bias(array[0], array[1]),
                                ldasv.rmse(array[0], array[1]),
                                ldasv.correlation(array[0], array[1]),
                                ldasv.correlation(array_ano[0], array_ano[1]),
                            )
                            std_an_val = np.std(array[0])
                            std_obs_val = np.std(array[1])
                            if (
                                (var == "SM")
                                and (cfg.SM_units == "Volumetric")
                                and cfg.Rescale_data
                            ):
                                rmsd_val = rmsd_val * (max_insitu - min_insitu)
                                bias_val = bias_val * (max_insitu - min_insitu)
                                std_an_val = std_an_val * (max_insitu - min_insitu)
                                std_obs_val = std_obs_val * (max_insitu - min_insitu)
                            ub_rmsd_val = np.sqrt(rmsd_val ** 2 - bias_val ** 2)

                            R[EXP][network + "_" + layer][season].append(corr_val)
                            Ano_R[EXP][network + "_" + layer][season].append(
                                ano_corr_val
                            )
                            RMSD[EXP][network + "_" + layer][season].append(rmsd_val)
                            ME[EXP][network + "_" + layer][season].append(bias_val)
                            Ub_RMSD[EXP][network + "_" + layer][season].append(
                                ub_rmsd_val
                            )

                            if (array[0].size > 2) and (array_ano[0].size > 2):
                                P_val = pearsonr(array[0], array[1])[1]
                                Ano_P_val = pearsonr(array_ano[0], array_ano[1])[1]
                            else:
                                P_val = pl.nan
                                Ano_P_val = pl.nan

                            P_value[EXP][network + "_" + layer][season].append(P_val)

                            Ano_P_value[EXP][network + "_" + layer][season].append(
                                Ano_P_val
                            )
                            STD_an[EXP][network + "_" + layer][season].append(
                                std_an_val
                            )
                            STD_obs[EXP][network + "_" + layer][season].append(
                                std_obs_val
                            )

                        # Calculate annual scores
                        for yr in years.year:

                            analysis_df_a = analysis_df[str(yr)]
                            insitu_df_a = insitu_df[str(yr)]
                            analysis_ano_a = analysis_ano[
                                np.where(analysis_df.index.to_series().dt.year == yr)[:]
                            ]
                            insitu_ano_a = insitu_ano[
                                np.where(analysis_df.index.to_series().dt.year == yr)[:]
                            ]

                            array = ldasv.filter_nan(
                                analysis_df_a.values, insitu_df_a.values
                            )
                            array_ano = ldasv.filter_nan(analysis_ano_a, insitu_ano_a)

                            bias_val, rmsd_val, corr_val, ano_corr_val = (
                                ldasv.bias(array[0], array[1]),
                                ldasv.rmse(array[0], array[1]),
                                ldasv.correlation(array[0], array[1]),
                                ldasv.correlation(array_ano[0], array_ano[1]),
                            )
                            std_an_val = np.std(array[0])
                            std_obs_val = np.std(array[1])
                            if (
                                (var == "SM")
                                and (cfg.SM_units == "Volumetric")
                                and cfg.Rescale_data
                            ):
                                rmsd_val = rmsd_val * (max_insitu - min_insitu)
                                bias_val = bias_val * (max_insitu - min_insitu)
                                std_an_val = std_an_val * (max_insitu - min_insitu)
                                std_obs_val = std_obs_val * (max_insitu - min_insitu)
                            ub_rmsd_val = np.sqrt(rmsd_val ** 2 - bias_val ** 2)

                            R[EXP][network + "_" + layer][str(yr)].append(corr_val)
                            Ano_R[EXP][network + "_" + layer][str(yr)].append(
                                ano_corr_val
                            )
                            RMSD[EXP][network + "_" + layer][str(yr)].append(rmsd_val)
                            ME[EXP][network + "_" + layer][str(yr)].append(bias_val)
                            Ub_RMSD[EXP][network + "_" + layer][str(yr)].append(
                                ub_rmsd_val
                            )
                            if (array[0].size > 2) and (array_ano[0].size > 2):
                                P_val = pearsonr(array[0], array[1])[1]
                                Ano_P_val = pearsonr(array_ano[0], array_ano[1])[1]
                            else:
                                P_val = pl.nan
                                Ano_P_val = pl.nan

                            P_value[EXP][network + "_" + layer][str(yr)].append(P_val)

                            Ano_P_value[EXP][network + "_" + layer][str(yr)].append(
                                Ano_P_val
                            )
                            STD_an[EXP][network + "_" + layer][str(yr)].append(
                                std_an_val
                            )
                            STD_obs[EXP][network + "_" + layer][str(yr)].append(
                                std_obs_val
                            )

                        if cfg.plot_time_series:

                            if (len(cfg.EXPVER) > 1) and (
                                EXP != (cfg.EXPVER[-1] + "_" + cfg.CLASS[-1])
                            ):

                                analysis_df.to_pickle(
                                    time_series_dir_net_layer_times_land
                                    + "/"
                                    + EXP
                                    + "_"
                                    + network
                                    + "_"
                                    + str(layer)
                                    + "_"
                                    + lat_analysis
                                    + "_"
                                    + lon_analysis
                                    + "_"
                                    + var
                                    + "_"
                                    + str(data_range[0].year)
                                    + "-%02d" % (data_range[0].month)
                                    + "-%02d" % (data_range[0].day)
                                    + "_"
                                    + str(data_range[-1].year)
                                    + "-%02d" % (data_range[-1].month)
                                    + "-%02d" % (data_range[-1].day)
                                    + "_"
                                    + times
                                    + "_"
                                    + land_type
                                    + "_R_anom_"
                                    + str(
                                        np.round(ano_corr_val,
                                            2
                                        )
                                    )
                                    )
                                    

                            if (var == "SM") and cfg.ST_quality_control:

                                if (len(cfg.EXPVER) > 1) and (
                                    EXP != (cfg.EXPVER[-1] + "_" + cfg.CLASS[-1])
                                ):
                                    ST_df.to_pickle(
                                        time_series_dir_net_layer_times_land
                                        + "/"
                                        + EXP
                                        + "_"
                                        + network
                                        + "_"
                                        + str(layer)
                                        + "_"
                                        + lat_analysis
                                        + "_"
                                        + lon_analysis
                                        + "_ST_"
                                        + str(data_range[0].year)
                                        + "-%02d" % (data_range[0].month)
                                        + "-%02d" % (data_range[0].day)
                                        + "_"
                                        + str(data_range[-1].year)
                                        + "-%02d" % (data_range[-1].month)
                                        + "-%02d" % (data_range[-1].day)
                                        + "_"
                                        + times
                                        + "_"
                                        + land_type
                                    )

                                else:


                                        ldasv.plot_SM_ST_time_series(
                                            analysis_df,
                                            insitu_df,
                                            ST_df,
                                            ST_insitu_df,
                                            network,
                                            lat_analysis,
                                            lon_analysis,
                                            layer,
                                            time_series_dir_net_layer_times_land,
                                            cfg.Rescale_data,
                                            data_range,
                                            EXP,
                                            [
                                                m + "_" + n
                                                for m, n in zip(cfg.EXPVER, cfg.CLASS)
                                            ],
                                            times,
                                            land_type,
                                            ano_corr_val,
                                        )
                            elif (EXP == (cfg.EXPVER[-1] + "_" + cfg.CLASS[-1])):

                                        ldasv.plot_time_series(
                                            var,
                                            analysis_df,
                                            insitu_df,
                                            network,
                                            lat_analysis,
                                            lon_analysis,
                                            layer,
                                            time_series_dir_net_layer_times_land,
                                            cfg.Rescale_data,
                                            data_range,
                                            EXP,
                                            [
                                                m + "_" + n
                                                for m, n in zip(cfg.EXPVER, cfg.CLASS)
                                            ],
                                            times,
                                            land_type,
                                            ano_corr_val,
                                        )

                    else:
                        # Remove point if not enough data available:
                        del latitude[-1]
                        del longitude[-1]

                n_stations[EXP][network + "_" + layer] = dict()

                # Warn if stations do not meet quality control

                if (var == "SM") and cfg.ST_quality_control:
                    for yr in years.year:
                        if (
                            Data_available[network][str(yr)] == True
                            and (Analysis_available[network][str(yr)] == True)
                            and (
                                sum(np.array(frozen_conditions[str(yr)][:]) == "0") > 0
                            )
                        ):
                            print(
                                "\n"
                                + str(
                                    sum(np.array(frozen_conditions[str(yr)][:]) == "0")
                                )
                                + " of "
                                + str(len(np.array(frozen_conditions[str(yr)][:])))
                                + " stations"
                                + " discarded for network "
                                + network
                                + " for year "
                                + str(yr)
                                + " due to observed soil temperature consistently below threshold of "
                                + str(cfg.ST_QC_threshold)
                                + " K."
                            )

                # Choose specific criteria for filtering scores i.e. P-Value<0.05:
                w[EXP][network + "_" + layer]["period"] = np.where(
                    (np.array(P_value[EXP][network + "_" + layer]["period"]) < 0.05)
                    & (
                        np.array(Ano_P_value[EXP][network + "_" + layer]["period"])
                        < 0.05
                    )
                    & (np.array(Ano_R[EXP][network + "_" + layer]["period"]) > -1.0)
                    & (np.array(Ano_R[EXP][network + "_" + layer]["period"]) < 1.0)
                )
                # Choose specific criteria for filtering seasonal scores i.e. P-Value<0.05:
                for season in slist:
                    w[EXP][network + "_" + layer][season] = np.where(
                        (np.array(P_value[EXP][network + "_" + layer][season]) < 0.05)
                        & (
                            np.array(Ano_P_value[EXP][network + "_" + layer][season])
                            < 0.05
                        )
                        & (np.array(Ano_R[EXP][network + "_" + layer][season]) > -1.0)
                        & (np.array(Ano_R[EXP][network + "_" + layer][season]) < 1.0)
                    )
                    lat_season = np.array(latitude)[
                        w[EXP][network + "_" + layer][season]
                    ]
                    lon_season = np.array(longitude)[
                        w[EXP][network + "_" + layer][season]
                    ]
                    lat_lon_networks = map("xx".join, zip(lat_season, lon_season))
                    n_stations[EXP][network + "_" + layer][season] = list(
                        lat_lon_networks
                    )  # removeDuplicates(lat_lon_networks)

                #    (np.array(R[network+'_'+layer][season])>0.0)&(np.array(Ano_R[network+'_'+layer][season]) > 0.0))

                # Choose specific criteria for filtering annual scores i.e. P-Value<0.05:
                for yr in years.year:
                    w[EXP][network + "_" + layer][str(yr)] = np.where(
                        (np.array(P_value[EXP][network + "_" + layer][str(yr)]) < 0.05)
                        & (
                            np.array(Ano_P_value[EXP][network + "_" + layer][str(yr)])
                            < 0.05
                        )
                        & (np.array(Ano_R[EXP][network + "_" + layer][str(yr)]) > -1.0)
                        & (np.array(Ano_R[EXP][network + "_" + layer][str(yr)]) < 1.0)
                    )

                    lat_year = np.array(latitude)[
                        w[EXP][network + "_" + layer][str(yr)]
                    ]
                    lon_year = np.array(longitude)[
                        w[EXP][network + "_" + layer][str(yr)]
                    ]
                    lat_lon_networks = map("xx".join, zip(lat_year, lon_year))
                    n_stations[EXP][network + "_" + layer][
                        str(yr)
                    ] = ldasv.removeDuplicates(list(lat_lon_networks))

                    if (
                        (len(w[EXP][network + "_" + layer][str(yr)][0]) == 0)
                        and (Data_available[network][str(yr)] == True)
                        and (Analysis_available[network][str(yr)] == True)
                    ):
                        print(
                            "\n Warning! Insufficient "
                            + var
                            + " observations available for the year "
                            + str(yr)
                            + ". Either try using "
                            + "different networks or extend the validation period."
                        )

                # Update latitude, longitude and R dicts for station map
                lat_period = np.array(latitude)[w[EXP][network + "_" + layer]["period"]]
                lon_period = np.array(longitude)[
                    w[EXP][network + "_" + layer]["period"]
                ]
                R_period = np.array(R[EXP][network + "_" + layer]["period"])[
                    w[EXP][network + "_" + layer]["period"]
                ]
                lat_lon_networks = map("xx".join, zip(lat_period, lon_period))
                xx = list(lat_lon_networks)
                yy = ldasv.removeDuplicates(xx)
                if len(xx) != len(yy):
                    print("duplicates - contact David F!!!")
                n_stations[EXP][network + "_" + layer]["period"] = yy

    # Make sure all experiments are validated for the same stations:

    n_stations["all_exp"] = dict()
    for n, network in enumerate(cfg.Network):
        for l, layer in enumerate(cfg.validation_layer):
            n_stations["all_exp"][network + "_" + layer] = dict()
            n_stations["all_exp"][network + "_" + layer]["period"] = list()

            for season in slist:
                n_stations["all_exp"][network + "_" + layer][season] = list()
            for yr in years.year:
                n_stations["all_exp"][network + "_" + layer][str(yr)] = list()

            for EXP in [m + "_" + n for m, n in zip(cfg.EXPVER, cfg.CLASS)]:
                if EXP == (cfg.EXPVER[0] + "_" + cfg.CLASS[0]):
                    n_stations["all_exp"][network + "_" + layer]["period"] = n_stations[
                        EXP
                    ][network + "_" + layer]["period"]
                    for season in slist:
                        n_stations["all_exp"][network + "_" + layer][
                            season
                        ] = n_stations[EXP][network + "_" + layer][season]
                    for yr in years.year:
                        n_stations["all_exp"][network + "_" + layer][
                            str(yr)
                        ] = n_stations[EXP][network + "_" + layer][str(yr)]
                else:
                    # Take common stations to enable period comparison between experiments:
                    if len(n_stations[EXP][network + "_" + layer]["period"]) > len(
                        n_stations["all_exp"][network + "_" + layer]["period"]
                    ):
                        n_stations["all_exp"][network + "_" + layer]["period"] = set(
                            n_stations[EXP][network + "_" + layer]["period"]
                        ).intersection(
                            n_stations["all_exp"][network + "_" + layer]["period"]
                        )
                    else:
                        n_stations["all_exp"][network + "_" + layer]["period"] = set(
                            n_stations["all_exp"][network + "_" + layer]["period"]
                        ).intersection(n_stations[EXP][network + "_" + layer]["period"])
                    for season in slist:
                        # Take common stations to enable annual comparison between experiments:
                        if len(n_stations[EXP][network + "_" + layer][season]) > len(
                            n_stations["all_exp"][network + "_" + layer][season]
                        ):
                            n_stations["all_exp"][network + "_" + layer][season] = set(
                                n_stations[EXP][network + "_" + layer][season]
                            ).intersection(
                                n_stations["all_exp"][network + "_" + layer][season]
                            )
                        else:
                            n_stations["all_exp"][network + "_" + layer][season] = set(
                                n_stations["all_exp"][network + "_" + layer][season]
                            ).intersection(
                                n_stations[EXP][network + "_" + layer][season]
                            )
                    for yr in years.year:
                        # Take common stations to enable seasonal comparison between experiments:
                        if len(n_stations[EXP][network + "_" + layer][str(yr)]) > len(
                            n_stations["all_exp"][network + "_" + layer][str(yr)]
                        ):
                            n_stations["all_exp"][network + "_" + layer][str(yr)] = set(
                                n_stations[EXP][network + "_" + layer][str(yr)]
                            ).intersection(
                                n_stations["all_exp"][network + "_" + layer][str(yr)]
                            )
                        else:
                            n_stations["all_exp"][network + "_" + layer][str(yr)] = set(
                                n_stations["all_exp"][network + "_" + layer][str(yr)]
                            ).intersection(
                                n_stations[EXP][network + "_" + layer][str(yr)]
                            )

    for n, network in enumerate(cfg.Network):
        for l, layer in enumerate(cfg.validation_layer):

            for EXP in [m + "_" + n for m, n in zip(cfg.EXPVER, cfg.CLASS)]:
                # Generate masks to filter common stations for each experiment validation (w[EXP][network+'_'+layer] represents
                # the mask of common stations for each experiment) over period/annual/seasonal dates:
                combined_stations = np.in1d(
                    n_stations[EXP][network + "_" + layer]["period"],
                    list(n_stations["all_exp"][network + "_" + layer]["period"]),
                )
                # Mask indices w[EXP][network+'_'+layer] for common stations over period:
                w[EXP][network + "_" + layer]["period"] = w[EXP][network + "_" + layer][
                    "period"
                ][0][np.where(combined_stations == True)]
                if list(n_stations[EXP][network + "_" + layer]["period"]):
                    # Redefine station locations in n_stations for each experiment to represent common stations over period:
                    n_stations[EXP][network + "_" + layer]["period"] = np.array(
                        n_stations[EXP][network + "_" + layer]["period"]
                    )[np.where(combined_stations == True)]
                for season in slist:
                    combined_stations = np.in1d(
                        n_stations[EXP][network + "_" + layer][season],
                        list(n_stations["all_exp"][network + "_" + layer][season]),
                    )
                    # Mask indices w[EXP][network+'_'+layer] for common stations over season:
                    w[EXP][network + "_" + layer][season] = w[EXP][
                        network + "_" + layer
                    ][season][0][np.where(combined_stations == True)]
                    if list(n_stations[EXP][network + "_" + layer][season]):
                        # Redefine station locations in n_stations for each experiment to represent common stations over season:
                        n_stations[EXP][network + "_" + layer][
                            season
                        ] = ldasv.removeDuplicates(
                            np.array(n_stations[EXP][network + "_" + layer][season])[
                                np.where(combined_stations == True)
                            ]
                        )
                for yr in years.year:
                    combined_stations = np.in1d(
                        n_stations[EXP][network + "_" + layer][str(yr)],
                        list(n_stations["all_exp"][network + "_" + layer][str(yr)]),
                    )
                    # Mask indices w[EXP][network+'_'+layer] for common stations over season:
                    w[EXP][network + "_" + layer][str(yr)] = w[EXP][
                        network + "_" + layer
                    ][str(yr)][0][np.where(combined_stations == True)]
                    # Redefine station locations in n_stations for each experiment to represent common stations over year:
                    n_stations[EXP][network + "_" + layer][str(yr)] = np.array(
                        n_stations[EXP][network + "_" + layer][str(yr)]
                    )[np.where(combined_stations == True)]

    # Generate tables and box plot scores for each network, each layer and each experiment:
    for n, network in enumerate(cfg.Network):
        for l, layer in enumerate(cfg.validation_layer):
            for exp, EXP in enumerate(
                [m + "_" + n for m, n in zip(cfg.EXPVER, cfg.CLASS)]
            ):
                # Generate station average R scores for maps
                lat_networks[EXP][network + "_" + layer]["period"] = [
                    x[: x.index("xx")] if "xx" in x else x
                    for x in n_stations["all_exp"][network + "_" + layer]["period"]
                ]
                lon_networks[EXP][network + "_" + layer]["period"] = [
                    x[x.index("xx") + 2 :] if "xx" in x else x
                    for x in n_stations["all_exp"][network + "_" + layer]["period"]
                ]
                R_period = np.array(R[EXP][network + "_" + layer]["period"])[
                    w[EXP][network + "_" + layer]["period"]
                ]
                R_networks[EXP][network + "_" + layer]["period"] = R_period
                Bias_period = np.array(ME[EXP][network + "_" + layer]["spring"])
                Bias_networks[EXP][network + "_" + layer]["period"] = Bias_period
                lat_networks_spring[EXP][network + "_" + layer]["period"] = [
                    x[: x.index("xx")] if "xx" in x else x
                    for x in n_stations["all_exp"][network + "_" + layer]["spring"]
                ]
                lon_networks_spring[EXP][network + "_" + layer]["period"] = [
                    x[x.index("xx") + 2 :] if "xx" in x else x
                    for x in n_stations["all_exp"][network + "_" + layer]["spring"]
                ]
                # Calculation of whether experiments perform significantly better than control (results in bar_chart_comparison plot):
                if exp > 0:
                    if cfg.stat_sig:
                        Ano_R_compare[EXP][network + "_" + layer]["period"] = (
                            np.array(
                                Ano_R[cfg.EXPVER[0] + "_" + cfg.CLASS[0]][
                                    network + "_" + layer
                                ]["period"]
                            )[
                                w[cfg.EXPVER[0] + "_" + cfg.CLASS[0]][
                                    network + "_" + layer
                                ]["period"]
                            ]
                            - np.array(Ano_R[EXP][network + "_" + layer]["period"])[
                                w[EXP][network + "_" + layer]["period"]
                            ]
                        )
                        Ano_R_compare[EXP][network + "_" + layer]["period"][
                            np.where(
                                Ano_R_compare[EXP][network + "_" + layer]["period"] > 0
                            )[:]
                        ] = 0
                        Ano_R_compare[EXP][network + "_" + layer]["period"][
                            np.where(
                                Ano_R_compare[EXP][network + "_" + layer]["period"] < 0
                            )[:]
                        ] = 1
                        Ano_R_conf_sig_1 = (
                            np.array(
                                Ano_R_conf_lower[cfg.EXPVER[0] + "_" + cfg.CLASS[0]][
                                    network + "_" + layer
                                ]["period"]
                            )[
                                w[cfg.EXPVER[0] + "_" + cfg.CLASS[0]][
                                    network + "_" + layer
                                ]["period"]
                            ]
                            - np.array(
                                Ano_R_conf_upper[EXP][network + "_" + layer]["period"]
                            )[w[EXP][network + "_" + layer]["period"]]
                        )
                        Ano_R_compare[EXP][network + "_" + layer]["period"][
                            np.where(Ano_R_conf_sig_1 > 0)[:]
                        ] = -1
                        Ano_R_conf_sig_2 = (
                            np.array(
                                Ano_R_conf_lower[EXP][network + "_" + layer]["period"]
                            )[w[EXP][network + "_" + layer]["period"]]
                            - np.array(
                                Ano_R_conf_upper[cfg.EXPVER[0] + "_" + cfg.CLASS[0]][
                                    network + "_" + layer
                                ]["period"]
                            )[
                                w[cfg.EXPVER[0] + "_" + cfg.CLASS[0]][
                                    network + "_" + layer
                                ]["period"]
                            ]
                        )
                        Ano_R_compare[EXP][network + "_" + layer]["period"][
                            np.where(Ano_R_conf_sig_2 > 0)[:]
                        ] = 2

                # Create the relevant post-processing directories
                Report_dir = (
                    cfg.output_dir
                    + "/"
                    + var
                    + "_validation_"
                    + EXP
                    + "_"
                    + start_period
                    + "_"
                    + end_period
                )
                tables_dir = Report_dir + "/tables/"
                plots_dir = Report_dir + "/plots/"
                time_series_dir = Report_dir + "/plots/time_series/"

                # Filter period and annual scores and put scores into table for each network and layer
                df = ldasv.draw_stats_table(
                    var,
                    R[EXP][network + "_" + layer],
                    Ano_R[EXP][network + "_" + layer],
                    RMSD[EXP][network + "_" + layer],
                    Ub_RMSD[EXP][network + "_" + layer],
                    ME[EXP][network + "_" + layer],
                    STD_an[EXP][network + "_" + layer],
                    STD_obs[EXP][network + "_" + layer],
                    w[EXP][network + "_" + layer],
                    ylist,
                    network,
                    cfg.Network,
                    tables_dir,
                    layer,
                    start_period,
                    end_period,
                    cfg.SM_units,
                    cfg.table_format,
                    EXP,
                    "annual",
                    n_stations[EXP][network + "_" + layer],
                    times,
                    land_type,
                    df,
                )

                # Draw annual and seasonal box plots:
                df = ldasv.draw_station_box_plots_sub_periods(
                    var,
                    EXP,
                    [R[EXP], Ano_R[EXP], RMSD[EXP], Ub_RMSD[EXP], ME[EXP]],
                    w[EXP],
                    plots_dir,
                    network,
                    cfg.Network,
                    layer,
                    years.year,
                    years.year,
                    "annual",
                    cfg.SM_units,
                    start_period,
                    end_period,
                    times,
                    land_type,
                    df,
                )
                df = ldasv.draw_station_box_plots_sub_periods(
                    var,
                    EXP,
                    [R[EXP], Ano_R[EXP], RMSD[EXP], Ub_RMSD[EXP], ME[EXP]],
                    w[EXP],
                    plots_dir,
                    network,
                    cfg.Network,
                    layer,
                    slist,
                    years.year,
                    "seasonal",
                    cfg.SM_units,
                    start_period,
                    end_period,
                    times,
                    land_type,
                    df,
                )

                # Filter seasonal scores and put scores into table for each network and layer
                df = ldasv.draw_stats_table(
                    var,
                    R[EXP][network + "_" + layer],
                    Ano_R[EXP][network + "_" + layer],
                    RMSD[EXP][network + "_" + layer],
                    Ub_RMSD[EXP][network + "_" + layer],
                    ME[EXP][network + "_" + layer],
                    STD_an[EXP][network + "_" + layer],
                    STD_obs[EXP][network + "_" + layer],
                    w[EXP][network + "_" + layer],
                    slist,
                    network,
                    cfg.Network,
                    tables_dir,
                    layer,
                    start_period,
                    end_period,
                    cfg.SM_units,
                    cfg.table_format,
                    EXP,
                    "seasonal",
                    n_stations[EXP][network + "_" + layer],
                    times,
                    land_type,
                    df,
                )

                # Draw box plot of CC and ACC over whole period (filtered scores for each network/layer):
                ldasv.draw_box_plot_period(
                    R[EXP][network + "_" + layer]["period"],
                    Ano_R[EXP][network + "_" + layer]["period"],
                    w[EXP][network + "_" + layer]["period"],
                    plots_dir,
                    network,
                    layer,
                    start_period,
                    end_period,
                    times,
                    land_type,
                )

    # Calculate scores (annual and seasonal) for all networks combined:
    for exp, EXP in enumerate([m + "_" + n for m, n in zip(cfg.EXPVER, cfg.CLASS)]):

        # Create the relevant post-processing directories
        Report_dir = (
            cfg.output_dir
            + "/"
            + var
            + "_validation_"
            + EXP
            + "_"
            + start_period
            + "_"
            + end_period
        )
        tables_dir = Report_dir + "/tables/"
        plots_dir = Report_dir + "/plots/"
        time_series_dir = Report_dir + "/plots/time_series/"

        # ldasv.draw_station_box_plots_period(R[EXP], Ano_R[EXP], w[EXP],plots_dir, cfg.Network, cfg.validation_layer, start_period, end_period)

        for layer in cfg.validation_layer:
            if var == "SM":
                units = cfg.SM_units
            elif var == "ST":
                units = cfg.ST_units

            df = ldasv.draw_station_box_plots_sub_periods(
                var,
                EXP,
                [R[EXP], Ano_R[EXP], RMSD[EXP], Ub_RMSD[EXP], ME[EXP]],
                w[EXP],
                plots_dir,
                "All_stations",
                cfg.Network,
                layer,
                years.year,
                years.year,
                "annual",
                units,
                start_period,
                end_period,
                times,
                land_type,
                df,
            )
            df = ldasv.draw_stats_table(
                var,
                R[EXP],
                Ano_R[EXP],
                RMSD[EXP],
                Ub_RMSD[EXP],
                ME[EXP],
                STD_an[EXP],
                STD_obs[EXP],
                w[EXP],
                ylist,
                "All_stations",
                cfg.Network,
                tables_dir,
                layer,
                start_period,
                end_period,
                units,
                cfg.table_format,
                EXP,
                "annual",
                n_stations[EXP],
                times,
                land_type,
                df,
            )
            df = ldasv.draw_stats_table(
                var,
                R[EXP],
                Ano_R[EXP],
                RMSD[EXP],
                Ub_RMSD[EXP],
                ME[EXP],
                STD_an[EXP],
                STD_obs[EXP],
                w[EXP],
                slist,
                "All_stations",
                cfg.Network,
                tables_dir,
                layer,
                start_period,
                end_period,
                units,
                cfg.table_format,
                EXP,
                "seasonal",
                n_stations[EXP],
                times,
                land_type,
                df,
            )

            if cfg.plot_time_series and (EXP == cfg.EXPVER[-1] + "_" + cfg.CLASS[-1]):

                for n, network in enumerate(cfg.Network):

                    time_series_dir_net_layer_land = (
                        time_series_dir
                        + network
                        + "/"
                        + layer
                        + "/"
                        + times
                        + "/"
                        + land_type
                    )

                    filename_html = "file://" + time_series_dir_net_layer_land

                    dic = {
                        "expver": [EXP],
                        "date": [start_period + " to " + end_period],
                        "type": ["Time series"],
                        "network": [network],
                        "variable": [layer + " " + var],
                        "metric": ["Time series plot"],
                        "link": filename_html,
                    }

                    df = df.append(pd.DataFrame(dic), ignore_index=True, sort=False)

            # Draw map of station scores

        df = ldasv.draw_station_map(
            df,
            start_period,
            end_period,
            lat_networks[EXP],
            lon_networks[EXP],
            R_networks[EXP],
            cfg.validation_layer,
            cfg.Network,
            plots_dir,
            var,
            EXP,
            times,
            land_type,
        )
        if cfg.stat_sig:
            if exp > 0:
                df = ldasv.draw_station_map(
                    df,
                    start_period,
                    end_period,
                    lat_networks[EXP],
                    lon_networks[EXP],
                    Ano_R_compare[EXP],
                    cfg.validation_layer,
                    cfg.Network,
                    plots_dir,
                    var,
                    EXP,
                    times,
                    land_type,
                    sig_comparisons=True,
                )

    if cfg.stat_sig:
        # Do statistical significance comparisons for Pearson R anomaly
        df = ldasv.error_bar_comparison(
            df,
            ([m + "_" + n for m, n in zip(cfg.EXPVER, cfg.CLASS)]),
            cfg.Network,
            Ano_R,
            w,
            Ano_R_conf_upper,
            Ano_R_conf_lower,
            plots_dir,
            tables_dir,
            cfg.validation_layer,
            "period",
            start_period,
            end_period,
            var,
            times,
            land_type,
        )

        if cfg.stat_sig:
            if exp > 0:

                df = ldasv.draw_station_box_plots_sub_periods_allExp(
                    ([m + "_" + n for m, n in zip(cfg.EXPVER, cfg.CLASS)]),
                    [R, Ano_R, RMSD, Ub_RMSD, ME],
                    w,
                    plots_dir,
                    cfg.Network,
                    cfg.validation_layer,
                    years.year,
                    years.year,
                    "annual",
                    units,
                    ([m + "_" + n for m, n in zip(cfg.EXPVER, cfg.CLASS)]),
                    var,
                    start_period,
                    end_period,
                    times,
                    land_type,
                    df,
                )

                df = ldasv.bar_chart_comparison(
                    df,
                    start_period,
                    end_period,
                    ([m + "_" + n for m, n in zip(cfg.EXPVER, cfg.CLASS)]),
                    cfg.Network,
                    Ano_R_compare,
                    plots_dir,
                    cfg.validation_layer,
                    "period",
                    times,
                    land_type,
                    var,
                )

    return df
