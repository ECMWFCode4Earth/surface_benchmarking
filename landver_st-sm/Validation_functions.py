"""
Functions for script LDAS_validate.py
Last modified by David Fairbairn on 12/8/19. 
"""

import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import pandas as pd
import re
import copy
from warnings import filterwarnings
import os
import cartopy.crs as ccrs
import cartopy as cart
import seaborn as sns
import glob

filterwarnings("ignore")
import math

numbers = re.compile(r"(\d+)")

###############################################################################
# Validation functions
###############################################################################


def Confidence_intervals(analysis_data, insitu_data, CC, sample_size):

    """
    Calculate lower and upper limit for confidence intervals using Fisher-Z transform with lag 1 autocorrelation
    (to remove the annual cycle correlations)
    """

    # Calculate lag-1 autocorrelation of time series
    insitu_autocorr = insitu_data.autocorr(lag=1)
    analysis_autocorr = analysis_data.autocorr(lag=1)
    # Calculate effective sample size from autocorrelation
    N_eff = (
        sample_size
        * (1.0 - (insitu_autocorr * analysis_autocorr))
        / (1.0 + (insitu_autocorr * analysis_autocorr))
    )

    Z_corr = math.log((1 + CC) / (1 - CC)) * 0.5

    std_eff = np.sqrt(1.0 / (N_eff - 3.0))
    Upper_z_limit = Z_corr + std_eff
    Lower_z_limit = Z_corr - std_eff
    Upper_r_limit = (math.exp(2.0 * Upper_z_limit) - 1.0) / (
        math.exp(2.0 * Upper_z_limit) + 1.0
    )
    Lower_r_limit = (math.exp(2.0 * Lower_z_limit) - 1.0) / (
        math.exp(2.0 * Lower_z_limit) + 1.0
    )

    return Lower_r_limit, Upper_r_limit


def removeDuplicates(listofElements):

    # Create an empty list to store unique elements
    uniqueList = []

    # Iterate over the original list and for each element
    # add it to uniqueList, if its not already there.
    for elem in listofElements:
        if elem not in uniqueList:
            uniqueList.append(elem)

    # Return the list of unique elements
    return uniqueList


def station_average_scores(listofElements, uniqueList, scores):

    average_scores = list()

    for elem in uniqueList:
        average_scores.append(
            np.nanmean(scores[np.where(np.array(listofElements) == elem)])
        )

    return average_scores


def filter_nan(s, o):
    """
    this functions removed the data  from simulated and observed data
    whereever the observed data contains nan
    this is used by all other functions, otherwise they will produce nan as
    output
    """
    data = np.array([s.flatten(), o.flatten()])
    data = np.transpose(data)
    data = data[~np.isnan(data).any(1)]
    return data[:, 0], data[:, 1]


def rmse(s, o):
    """
    Root Mean Squared Error
    input:
    s: simulated
    o: observed
    output:
    rmses: root mean squared error
    """
    s, o = filter_nan(s, o)
    return np.sqrt(np.mean((s - o) ** 2))


def correlation(s, o):
    """
    correlation coefficient
    input:
    s: simulated
    o: observed
    output:
    correlation: correlation coefficient
    """
    s, o = filter_nan(s, o)
    if s.size == 0:
        corr = np.NaN
    else:
        corr = np.corrcoef(o, s)[0, 1]
    return corr


def bias(s, o):
    """
    Bias
    input:
    s: simulated
    o: observed
    output:
    bias: bias
    """
    s, o = filter_nan(s, o)
    return np.mean(s - o)


def find_lat_lon(file_name, data_type):
    """
    Find latitude and longitude for file name
    input:
    s: file_name
    output:
    latitude,longitude
    """

    if data_type == "analysis":
        splitfile = re.split("[_,]", file_name)
        lat = splitfile[-2]
        lon = splitfile[-1][:-4]
    elif data_type == "insitu":
        splitfile = re.split("[_,]", file_name)
        lat = splitfile[-4]
        lon = splitfile[-3]
    return lat, lon


def network_SM_depths(network):
    if network == "USCRN":
        depths = [0.05, 0.10, 0.20, 0.50, 1.00]

    elif network == "SCAN":
        depths = [0.05, 0.10, 0.20, 0.50, 1.00]

    elif network == "SNOTEL":
        depths = [0.05, 0.20, 0.50]

    elif network == "SMOSMANIA":
        depths = [0.05, 0.10, 0.20, 0.30]

    elif network == "OZNET":
        depths = [0.04, 0.15, 0.45, 0.75]

    elif network == "REMEDHUS":
        depths = [0.05]

    elif network == "TERENO":
        depths = [0.05, 0.20, 0.50]

    # SM-DAS2 soil moisture depths
    analysis_depths = [0.07, 0.21, 0.72]

    return depths, analysis_depths


def land_class_lookup_table(land_type):
    land_dict = {"Croplands": [10, 12], "Grassland": [130]}

    return land_dict[land_type]


###############################################################################
# Functions to read analysis and observation time series
###############################################################################


def trapezoidal_method(depths, in_situ_obs):
    """
    Coversion of in situ observations to root-zone soil moisture
    input:
    obs: point observations at the sites
    depths: depths of the point observations
    output:
    RZ_obs: Root-zone integrated column
    """

    # Use trapezoidal method to calculated root-zone soil moisture from
    # point observations. Assume obs at 0.05 cm are also at surface
    RZ_obs = 0.5 * (in_situ_obs.loc[:, 0] + in_situ_obs.loc[:, 0]) * (depths[0] - 0.0)
    layer_total = depths[0] - 0.0

    for l in range(len(depths) - 1):
        RZ_obs = RZ_obs + 0.5 * (in_situ_obs.loc[:, l] + in_situ_obs.loc[:, l + 1]) * (
            depths[l + 1] - depths[l]
        )
        layer_total = layer_total + (depths[l + 1] - depths[l])
    if depths[-1] < 1.0:
        RZ_obs = RZ_obs + 0.5 * (
            in_situ_obs.loc[:, l + 1] + in_situ_obs.loc[:, l + 1]
        ) * (1.0 - depths[l + 1])
        layer_total = layer_total + (1.0 - depths[l + 1])

    RZ_obs = RZ_obs / layer_total

    return RZ_obs


def insitu_df_to_numeric(
    Data_dir, Yr, Network, Layer, Lat, Lon, Var_type, flag_type, ST_qc=False
):

    Insitu_yr = pd.read_pickle(
        Data_dir
        + "/"
        + str(Yr)
        + "/"
        + Network
        + "_"
        + "%.2f" % (Layer)
        + "_"
        + Lat
        + "_"
        + Lon
        + "_obs_"
        + str(Yr)
    )[Var_type]

    Insitu_yr_flag = pd.read_pickle(
        Data_dir
        + "/"
        + str(Yr)
        + "/"
        + Network
        + "_"
        + "%.2f" % (Layer)
        + "_"
        + Lat
        + "_"
        + Lon
        + "_obs_"
        + str(Yr)
    )[flag_type]

    return Insitu_yr, Insitu_yr_flag


def define_start_end(Yr, Data_range, Years):

    if (Yr == Data_range[0].year) and (len(Years) > 1):
        start = (
            str(Data_range[0].year)
            + "-%02d" % (Data_range[0].month)
            + "-%02d" % (Data_range[0].day)
            + "-%02d" % (Data_range[0].hour)
        )
        end = str(Data_range[0].year) + "-12-31" + "-%02d" % (Data_range[-1].hour)
    elif (Yr > Data_range[0].year) and (Yr < Data_range[-1].year) and (len(Years) > 1):
        start = str(Yr) + "-01-01" + "-%02d" % (Data_range[0].hour)
        end = str(Yr) + "-12-31" + "-%02d" % (Data_range[-1].hour)
    elif len(Years) == 1:
        start = (
            str(Data_range[0].year)
            + "-%02d" % (Data_range[0].month)
            + "-%02d" % (Data_range[0].day)
            + "-%02d" % (Data_range[0].hour)
        )
        end = (
            str(Data_range[-1].year)
            + "-%02d" % (Data_range[-1].month)
            + "-%02d" % (Data_range[-1].day)
            + "-%02d" % (Data_range[-1].hour)
        )
    else:
        start = str(Yr) + "-01-01" + "-%02d" % (Data_range[0].hour)
        end = (
            str(Data_range[-1].year)
            + "-%02d" % (Data_range[-1].month)
            + "-%02d" % (Data_range[-1].day)
            + "-%02d" % (Data_range[-1].hour)
        )
    return start, end


def extract_insitu_start_end(Daily_obs_time_average, Insitu_yr, Data_range, Start, End):

    if Daily_obs_time_average:
        insitu_df_layer = (
            Insitu_yr.asfreq("H")
            .loc[Start:End]
            .resample(Data_range.freq, base=12)
            .mean()
            .values
        )
    else:
        insitu_df_layer = Insitu_yr.asfreq("H").loc[Start:End].asfreq(Data_range.freq)
        if np.nansum(insitu_df_layer[:]) == 0.0:
            insitu_df_layer = (
                Insitu_yr.asfreq("H")
                .loc[Start:End]
                .resample("3H")
                .mean()
                .asfreq(Data_range.freq)
            )

    return insitu_df_layer


def read_and_rescale_insitu_data(
    var_type,
    data_dir,
    years,
    network,
    lat,
    lon,
    depths,
    data_range,
    ST_quality_control=False,
    daily_obs_time_average=False,
    Rescale_data=True,
    ST_QC_threshold_celcius=None,
):

    #    Optional outputs:
    ST_in_situ = None
    int_plus = None
    int_minus = None

    Data_df_layer = pd.DataFrame(pl.empty((data_range.size, len(depths))), data_range)
    if (var_type == "SM") and ST_quality_control:
        ST_df_layer = pd.DataFrame(pl.empty((data_range.size, len(depths))), data_range)

    for l, layer in enumerate(depths):
        for yr in years:
            if var_type == "SM":
                try:
                    insitu_yr, Insitu_yr_flag = insitu_df_to_numeric(
                        data_dir,
                        yr,
                        network,
                        layer,
                        lat,
                        lon,
                        var_type,
                        "SM_flag",
                        ST_quality_control,
                    )
                except:
                    Data_df_layer.loc[str(yr), l] = pl.nan
                    if ST_quality_control:
                        ST_df_layer.loc[str(yr), l] = pl.nan
                    continue
                if ST_quality_control:
                    ST_insitu_yr, ST_insitu_yr_flag = insitu_df_to_numeric(
                        data_dir,
                        yr,
                        network,
                        layer,
                        lat,
                        lon,
                        "TS",
                        "TS_flag",
                        ST_quality_control,
                    )
                    ST_insitu_yr[ST_insitu_yr_flag != "G"] = pl.nan
                    # Screen soil moisture where temperature is near-freezing
                    insitu_yr[ST_insitu_yr < ST_QC_threshold_celcius] = pl.nan
                    insitu_yr[pd.isnull(ST_insitu_yr)] = pl.nan

                # Screen erroneous SM observations
                insitu_yr[insitu_yr <= 0.0001] = pl.nan
                insitu_yr[insitu_yr > 1.0] = pl.nan

            elif var_type == "ST":
                try:
                    insitu_yr, insitu_yr_flag = insitu_df_to_numeric(
                        data_dir,
                        yr,
                        network,
                        layer,
                        lat,
                        lon,
                        "TS",
                        "TS_flag",
                        ST_quality_control,
                    )
                    insitu_yr[insitu_yr_flag != "G"] = pl.nan
                except:
                    Data_df_layer.loc[str(yr), l] = pl.nan
                    continue

            start, end = define_start_end(yr, data_range, years)
            try:
                Data_df_layer.loc[start:end, l] = extract_insitu_start_end(
                    daily_obs_time_average, insitu_yr, data_range, start, end
                )
            except:
                raise Exception(
                    "Problem with in situ data extraction for "
                    + lat
                    + "_"
                    + lon
                    + "_"
                    + " from "
                    + start
                    + " to "
                    + end
                    + " - contact David Fairbairn (david.fairbairn@ecmwf.int)"
                )

            if (var_type == "SM") and ST_quality_control:
                try:
                    ST_df_layer.loc[start:end, l] = (
                        ST_insitu_yr.asfreq("H")
                        .loc[start:end]
                        .resample("3H")
                        .mean()
                        .asfreq(data_range.freq)
                    )
                except:
                    ST_df_layer.loc[start:end, l] = (
                        ST_insitu_yr.asfreq("H")
                        .loc[start:end]
                        .resample("3H")
                        .mean()
                        .asfreq(data_range.freq)
                    )

    # Use trapezoidal method to integrate over root-zone:
    if len(depths) > 1:
        Data_in_situ = trapezoidal_method(depths, Data_df_layer)
        if (var_type == "SM") and ST_quality_control:
            ST_in_situ = trapezoidal_method(depths, ST_df_layer)
    else:
        Data_in_situ = Data_df_layer.loc[:, 0]
        if (var_type == "SM") and ST_quality_control:
            ST_in_situ = ST_df_layer.loc[:, 0]

    # Rescale if set for SM and return in situ dataframe:
    if (var_type == "SM") and ST_quality_control:
        int_plus = Data_in_situ.loc[:].max()
        int_minus = Data_in_situ.loc[:].min()
        if Rescale_data:
            Data_in_situ.loc[:] = (Data_in_situ.loc[:] - int_minus) / (
                int_plus - int_minus
            )

    elif (var_type == "SM") and not ST_quality_control:
        int_plus = Data_in_situ.loc[:].max()
        int_minus = Data_in_situ.loc[:].min()
        if Rescale_data:
            Data_in_situ.loc[:] = (Data_in_situ.loc[:] - int_minus) / (
                int_plus - int_minus
            )

    return Data_in_situ, ST_in_situ, int_plus, int_minus


def read_and_rescale_analysis_data(
    var_type,
    analysis_dir,
    i,
    years,
    lat,
    lon,
    depths,
    data_range,
    time_freq,
    ST_quality_control,
    Rescale_data,
    soilML10,
    ST_QC_threshold_kelvin=None,
):
    #   Optional output for soil temperature quality control
    ST_df_avg = None

    Data_df_layer = pd.DataFrame(pl.empty((data_range.size, len(depths))), data_range)
    if (var_type == "SM") and ST_quality_control:
        ST_df_layer = pd.DataFrame(pl.empty((data_range.size, len(depths))), data_range)

    for yr in years:

        try:
            analysis_yr = pd.read_csv(
                analysis_dir
                + "/"
                + str(yr)
                + "/"
                + var_type
                + "_"
                + lat
                + ","
                + lon
                + ".dat",
                header=None,
                delimiter=" ",
                index_col=0,
            )
            analysis_yr.index = analysis_yr.index.map(str)
            analysis_yr.index = pd.to_datetime(analysis_yr.index, format="%Y%m%d%H")
        except:
            Data_df_layer.loc[str(yr), :] = pl.nan
            if (var_type == "SM") and ST_quality_control:
                ST_df_layer.loc[str(yr), :] = pl.nan
            continue

        start, end = define_start_end(yr, data_range, years)
        try:
            if (var_type == "SM") or not soilML10:
                Data_df_layer.loc[start:end] = analysis_yr.loc[
                    start:end, 1 : len(depths)
                ].asfreq(time_freq)

            else:
                Data_df_layer.loc[start:end] = analysis_yr.loc[start:end, 4:4].asfreq(
                    time_freq
                )

        except:
            raise Exception(
                "Time period of preprocessed analysis "
                + var_type
                + " is not included in time period in Options namelist. Reprocess analysis "
                + var_type
                + " to include analysis time period."
            )

        if (var_type == "SM") and ST_quality_control:
            ST_ana_yr = pd.read_csv(
                analysis_dir + "/" + str(yr) + "/ST_" + lat + "," + lon + ".dat",
                header=None,
                delimiter=" ",
                index_col=0,
            )
            ST_ana_yr.index = ST_ana_yr.index.map(str)
            ST_ana_yr.index = pd.to_datetime(ST_ana_yr.index, format="%Y%m%d%H")
            try:
                ST_df_layer.loc[start:end] = ST_ana_yr.loc[
                    start:end, 1 : len(depths)
                ].asfreq(time_freq)
            except:
                raise Exception(
                    "Time period of preprocessed analysis ST is not included in time period in Options namelist. Reprocess analysis ST to include analysis time period."
                )

    for l, layer in enumerate(depths):

        if var_type == "SM":
            Data_df_layer.loc[:, :][Data_df_layer.loc[:, l] <= 0.0001] = pl.nan
            Data_df_layer.loc[:, :][Data_df_layer.loc[:, l] > 1.0] = pl.nan
        elif var_type == "ST":
            Data_df_layer.loc[:, :][Data_df_layer.loc[:, 0] > 1000.0] = pl.nan

        if (var_type == "SM") and ST_quality_control:
            Data_df_layer.loc[:, :][
                ST_df_layer.loc[:, l] < ST_QC_threshold_kelvin
            ] = pl.nan
            Data_df_layer.loc[:, :][pl.isnan(ST_df_layer.loc[:, l])] = pl.nan
            ST_df_layer.loc[:, l] = ST_df_layer.loc[:, l] * layer

        Data_df_layer.loc[:, :][pl.isnan(Data_df_layer.loc[:, 0])] = pl.nan
        Data_df_layer.loc[:, l] = Data_df_layer.loc[:, l] * layer

    if var_type == "SM":
        Data_analysis = Data_df_layer.sum(skipna=False, axis=1) / sum(depths)
    elif var_type == "ST":
        Data_analysis = Data_df_layer.sum(skipna=False, axis=1) / sum(depths) - 273.15

    # Rescale analysis SM to produce soil water index:
    if (var_type == "SM") and Rescale_data:
        Data_analysis.loc[:] = (Data_analysis.loc[:] - Data_analysis.loc[:].min()) / (
            Data_analysis.loc[:].max() - Data_analysis.loc[:].min()
        )

    if (var_type == "SM") and ST_quality_control:
        ST_df_avg = ST_df_layer.sum(skipna=False, axis=1) / sum(depths) - 273.15

    return Data_analysis, ST_df_avg


def numericalSort(value):
    parts = numbers.split(value)
    return parts


###############################################################################
# Plotting functions:
###############################################################################

""" 
 Plot single time series of soil moisture or soil temperature
"""

def plot_time_series(
    var_type,
    analysis,
    insitu,
    Network_name,
    lat_val,
    lon_val,
    layer_val,
    Plot_dir,
    SM_rescale,
    years,
    exp_cur,
    EXP_list,
    times,
    land_type,
    Anno_Corr_value,
):

    fig, axes = pl.subplots(nrows=1, ncols=1)
    fig.set_figheight(7)
    fig.set_figwidth(9)

    # If there is more than one experiment load pickle file for each previous experiment:
    if exp_cur != EXP_list[0]:
        var_exp = dict()
        Ranom_exp = list()
        # Find pickle file to plot for previous experiment:

        for exp_val, exp_name in enumerate(EXP_list[:-1]):

           # Find pickle files to plot for previous experiments:
            sample_colors = ["g", "k", "c", "m", "r", "b"]
            var_pickle = (
                Plot_dir
                + "/"
                + exp_name
                + "_"
                + Network_name
                + "_"
                + str(layer_val)
                + "_"
                + lat_val
                + "_"
                + lon_val
                + "_"
                + var_type
                + "_"
                + str(years[0].year)
                + "-%02d" % (years[0].month)
                + "-%02d" % (years[0].day)
                + "_"
                + str(years[-1].year)
                + "-%02d" % (years[-1].month)
                + "-%02d" % (years[-1].day)
                + "_"
                + times
                + "_"
                + land_type
                + "*")

            start_index = var_pickle.find("validation")

            var_pickle_end = var_pickle[start_index:].replace(EXP_list[-1], exp_name)
            var_pickle_start = var_pickle[:start_index]

            var_pickle_file = var_pickle_start + var_pickle_end

            if exp_val==0:
              var_pickle_control_file=var_pickle_file

            if glob.glob(var_pickle_file) and glob.glob(var_pickle_control_file):

                if exp_val==0:
                    var_pickle_control = glob.glob(var_pickle_file)[0]
                    var_pickle=copy.deepcopy(var_pickle_control)
                    R_index=var_pickle_control.find("R_anom_")
                    Ranom_control=float(var_pickle_control[R_index+7:])
                else:
                    var_pickle = glob.glob(var_pickle_file)[0]
                    R_index=var_pickle.find("R_anom_")
                    Ranom_exp.append(exp_name+"_"+str(np.round(float(var_pickle[R_index+7:]) - Ranom_control,2)))


                var_exp[exp_val] = pd.read_pickle(var_pickle)
                var_exp[exp_val].plot(
                    ax=axes,
                    linewidth=1,
                    linestyle=":",
                    color=sample_colors[exp_val],
                    label=exp_name,
                )
                if exp_val>0:
                  os.remove(var_pickle)

        if glob.glob(var_pickle_control_file):
          # Calculate Pearson R annom of current experiment vs control and append to list
          Ranom_exp.append(exp_cur+"_"+str(np.round(Anno_Corr_value - Ranom_control,2)))
          os.remove(var_pickle_control)
    else:
        # If only one experiment is validated, then store anomaly CC of this this experiment
        Ranom_exp.append(exp_cur+"_"+str(np.round(Anno_Corr_value,2)))

    analysis.plot(ax=axes, linewidth=1, linestyle="--", color="r", label=exp_cur)

    insitu.plot(ax=axes, linewidth=1, marker="o", color="b", label="In situ")

    axes.set_title(var_type + " at:" + lat_val + ", lon=" + lon_val, fontsize=10)
    if (var_type == "SM") and SM_rescale:
        axes.set_ylabel("SM (-)")
    elif (var_type == "SM") and not SM_rescale:
        axes.set_ylabel("SM ($m^3/m^3$)")
    elif var_type == "ST":
        axes.set_ylabel("Temperature  ($^{\circ}C$)")

    axes.set_xlabel("Time")

    axes.tick_params(axis="y", labelsize=9)
    axes.tick_params(axis="x", labelsize=9)
    axes.xaxis.label.set_fontsize(10)
    axes.yaxis.label.set_fontsize(10)

    handles, labels = axes.get_legend_handles_labels()
    plt.legend(handles, labels)

    pl.savefig(
        Plot_dir
        + "/"
        + "Rannom_vs_control_"
        + '_'.join(Ranom_exp)
        + "_"
        + "_lat_"
        + lat_val
        + "_lon_"
        + lon_val
        + "_"
        + var_type
        + "_"
        + Network_name
        + "_"
        + str(layer_val)
        + "_"
        + str(years[0].year)
        + "-%02d" % (years[0].month)
        + "-%02d" % (years[0].day)
        + "_"
        + str(years[-1].year)
        + "-%02d" % (years[-1].month)
        + "-%02d" % (years[-1].day)
        + "_"
        + times
        + "_"
        + land_type
        + ".pdf"
    )
    pl.close()

""" 
 Plot soil moisture time series with soil temperature quality control
"""

def plot_SM_ST_time_series(
    analysis_SM,
    in_situ_SM,
    analysis_ST,
    in_situ_ST,
    Network_name,
    lat_val,
    lon_val,
    layer_val,
    Plot_dir,
    SM_rescale,
    years,
    exp_cur,
    EXP_list,
    times,
    land_type,
    Anno_Corr_value,
):

    fig, axes = pl.subplots(nrows=2, ncols=1)
    fig.set_figheight(9)
    fig.set_figwidth(9)

    # If there is more than one experiment load pickle file for each previous experiment:
    if exp_cur != EXP_list[0]:
        SM_exp = dict()
        ST_exp = dict()
        Ranom_exp = list()

        # Find SM pickle files to plot for previous experiments:
        for exp_val, exp_name in enumerate(EXP_list[:-1]):
            sample_colors = ["g", "k", "c", "m", "r", "b"]

            SM_pickle_file=(Plot_dir
                + "/"
                + exp_name
                + "_"
                + Network_name
                + "_"
                + str(layer_val)
                + "_"
                + lat_val
                + "_"
                + lon_val
                + "_SM_"
                + str(years[0].year)
                + "-%02d" % (years[0].month)
                + "-%02d" % (years[0].day)
                + "_"
                + str(years[-1].year)
                + "-%02d" % (years[-1].month)
                + "-%02d" % (years[-1].day)
                + "_"
                + times
                + "_"
                + land_type
                + "*")
            
            start_index = SM_pickle_file.find("validation")

            SM_pickle_end = SM_pickle_file[start_index:].replace(EXP_list[-1], exp_name)

            SM_pickle_start = SM_pickle_file[:start_index]

            SM_pickle_file = SM_pickle_start + SM_pickle_end

            # Calculate Pearson R annom of previous experiment vs control and append to list
            if exp_val==0:
              SM_pickle_control_file=SM_pickle_file

            if glob.glob(SM_pickle_file) and glob.glob(SM_pickle_control_file):

                if exp_val==0:
                    SM_pickle_control = glob.glob(SM_pickle_file)[0]
                    SM_pickle=copy.deepcopy(SM_pickle_control)
                    R_index=SM_pickle_control.find("R_anom_")
                    try:
                      Ranom_control=float(SM_pickle_control[R_index+7:])
                    except:
                      import pdb; pdb.set_trace()
                else:
                    SM_pickle = glob.glob(SM_pickle_file)[0]
                    R_index=SM_pickle.find("R_anom_")
                    Ranom_exp.append(exp_name+"_"+str(np.round(float(SM_pickle[R_index+7:]) - Ranom_control,2)))

                # Plot soil moisture time series for previous experiments
                SM_exp[exp_val] = pd.read_pickle(SM_pickle)
                SM_exp[exp_val].plot(
                    ax=axes[0],
                    linewidth=1,
                    linestyle=":",
                    color=sample_colors[exp_val],
                    label=exp_name,
                )
                if exp_val>0:
                  os.remove(SM_pickle)


                # Find ST pickle files to plot for previous experiments:
                ST_pickle = (
                    Plot_dir
                    + "/"
                    + exp_name
                    + "_"
                    + Network_name
                    + "_"
                    + str(layer_val)
                    + "_"
                    + lat_val
                    + "_"
                    + lon_val
                    + "_ST_"
                    + str(years[0].year)
                    + "-%02d" % (years[0].month)
                    + "-%02d" % (years[0].day)
                    + "_"
                    + str(years[-1].year)
                    + "-%02d" % (years[-1].month)
                    + "-%02d" % (years[-1].day)
                    + "_"
                    + times
                    + "_"
                    + land_type
                )

                # Plot soil temperature time series for previous experiments
                ST_pickle_end = ST_pickle[start_index:].replace(EXP_list[-1], exp_name)
                ST_pickle_start = ST_pickle[:start_index]

                ST_pickle = ST_pickle_start + ST_pickle_end

                ST_exp[exp_val] = pd.read_pickle(ST_pickle)
                ST_exp[exp_val].plot(
                    ax=axes[1],
                    linewidth=1,
                    linestyle=":",
                    color=sample_colors[exp_val],
                    label=exp_name,
                )
                os.remove(ST_pickle)


        if glob.glob(SM_pickle_control_file):
          # Calculate Pearson R annom of current experiment vs control and append to list
          Ranom_exp.append(exp_cur+"_"+str(np.round(Anno_Corr_value - Ranom_control,2)))
          os.remove(SM_pickle_control)
    else:
        # If only one experiment is validated, then store anomaly CC of this this experiment
        Ranom_exp.append(exp_cur+"_"+str(np.round(Anno_Corr_value,2)))
        
    # Plot SM analysis and in situ data for latest experiment:
    analysis_SM.plot(ax=axes[0], linewidth=1, linestyle="--", color="r", label=exp_cur)
    in_situ_SM.plot(ax=axes[0], linewidth=1, marker="o", color="b", label="In situ")

    # Plot ST analysis and in situ data for latest experiment:
    analysis_ST.plot(ax=axes[1], linewidth=1, linestyle="--", color="r", label=exp_cur)
    in_situ_ST.plot(ax=axes[1], linewidth=1, marker="o", color="b", label="In situ")

    # Label SM plot axes
    axes[0].set_title("SM at:" + lat_val + ", lon=" + lon_val, fontsize=10)
    if SM_rescale:
        axes[0].set_ylabel("SM (-)")
    else:
        axes[0].set_ylabel("SM ($m^3/m^3$)")

    handles, labels = axes[0].get_legend_handles_labels()
    # legend = plt.legend(handles,labels)

    axes[0].legend(loc=0)
    axes[1].legend(loc=0)

    # handles, labels = axes[1].get_legend_handles_labels()
    # legend = plt.legend(handles,labels)

    # Label ST plot axes
    axes[1].set_title("ST at: lat=" + lat_val + ", lon=" + lon_val, fontsize=10)
    axes[1].set_ylabel("Temperature  ($^{\circ}C$)")
    axes[1].set_xlabel("Time")

    axes[0].tick_params(axis="y", labelsize=9)
    axes[0].tick_params(axis="x", labelsize=9)
    axes[0].xaxis.label.set_fontsize(10)
    axes[0].yaxis.label.set_fontsize(10)

    axes[1].tick_params(axis="y", labelsize=9)
    axes[1].tick_params(axis="x", labelsize=9)
    axes[1].xaxis.label.set_fontsize(10)
    axes[1].yaxis.label.set_fontsize(10)

    # Save figure to pdf with filename starting with Pearson R anom differences
    pl.savefig(
        Plot_dir
        + "/"
        + "Rannom_vs_control_"
        + '_'.join(Ranom_exp)
        + "_lat_"
        + lat_val
        + "_lon_"
        + lon_val
        + "_"
        + Network_name
        + "_"
        + str(layer_val)
        + "_SM_ST_"
        + str(years[0].year)
        + "-%02d" % (years[0].month)
        + "-%02d" % (years[0].day)
        + "_"
        + str(years[-1].year)
        + "-%02d" % (years[-1].month)
        + "-%02d" % (years[-1].day)
        + "_"
        + times
        + "_"
        + land_type
        + ".pdf"
    )
    pl.close()


def draw_stats_table(
    var_type,
    CC,
    ACC,
    RMSD,
    Ub_RMSD,
    BIAS,
    STD_an,
    STD_obs,
    w,
    time_periods,
    station_type,
    NET,
    table_dir,
    layer,
    start_year,
    end_year,
    SM_units,
    table_format,
    exp_name,
    plot_type,
    n_stations,
    times,
    land_type,
    df,
):

    nrows = 8
    ncols = len(time_periods) + 1

    fig = plt.figure(figsize=(1.0 * ncols, 0.4 * (nrows + 1)))
    ax = fig.add_subplot(111)

    ax.axis("off")

    CCs = dict()
    ACCs = dict()
    RMSDs = dict()
    Ub_RMSDs = dict()
    BIASs = dict()
    STDs_an = dict()
    STDs_obs = dict()
    col_labels = list()
    table_vals = list()
    n_stations_cum = dict()

    timespan = copy.deepcopy(time_periods)
    timespan.append("period")

    if var_type == "SM":
        if SM_units == "Volumetric":
            unit = "m$^3$/m$^3$"
        else:
            unit = "-"
    else:
        unit = "K"

    if station_type == "All_stations":
        # Calculate filtered scores for all networks

        for t in timespan:
            CCs[t] = list()
            ACCs[t] = list()
            RMSDs[t] = list()
            Ub_RMSDs[t] = list()
            BIASs[t] = list()
            STDs_an[t] = list()
            STDs_obs[t] = list()
            n_stations_cum[t] = list()
        # Concatenate the scores over all the networks:
        for n, network in enumerate(NET):
            for t in timespan:
                CCs[t] = np.append(
                    np.array(CCs[t]),
                    np.array(CC[NET[n] + "_" + layer][t])[w[NET[n] + "_" + layer][t]],
                )
                ACCs[t] = np.append(
                    np.array(ACCs[t]),
                    np.array(ACC[NET[n] + "_" + layer][t])[w[NET[n] + "_" + layer][t]],
                )
                RMSDs[t] = np.append(
                    np.array(RMSDs[t]),
                    np.array(RMSD[NET[n] + "_" + layer][t])[w[NET[n] + "_" + layer][t]],
                )
                Ub_RMSDs[t] = np.append(
                    np.array(Ub_RMSDs[t]),
                    np.array(Ub_RMSD[NET[n] + "_" + layer][t])[
                        w[NET[n] + "_" + layer][t]
                    ],
                )
                BIASs[t] = np.append(
                    np.array(BIASs[t]),
                    np.array(BIAS[NET[n] + "_" + layer][t])[w[NET[n] + "_" + layer][t]],
                )
                STDs_an[t] = np.append(
                    np.array(STDs_an[t]),
                    np.array(STD_an[NET[n] + "_" + layer][t])[
                        w[NET[n] + "_" + layer][t]
                    ],
                )
                STDs_obs[t] = np.append(
                    np.array(STDs_obs[t]),
                    np.array(STD_obs[NET[n] + "_" + layer][t])[
                        w[NET[n] + "_" + layer][t]
                    ],
                )
                n_stations_cum[t] = np.append(
                    np.array(n_stations_cum[t]),
                    np.array(n_stations[NET[n] + "_" + layer][t]),
                )
    else:
        for t in timespan:
            CCs[t] = np.array(CC[t])[w[t]]
            ACCs[t] = np.array(ACC[t])[w[t]]
            RMSDs[t] = np.array(RMSD[t])[w[t]]
            Ub_RMSDs[t] = np.array(Ub_RMSD[t])[w[t]]
            BIASs[t] = np.array(BIAS[t])[w[t]]
            STDs_an[t] = np.array(STD_an[t])[w[t]]
            STDs_obs[t] = np.array(STD_obs[t])[w[t]]
            n_stations_cum[t] = np.array(n_stations[t])
    for t in timespan:
        col_labels = np.append(
            col_labels, t + " (" + str(np.size(n_stations_cum[t])) + " stations)"
        )
        table_vals = np.append(
            table_vals,
            [
                np.array(
                    [
                        np.around(np.mean(CCs[t]), 5),
                        np.around(np.mean(ACCs[t]), 5),
                        np.around(np.mean(RMSDs[t]), 5),
                        np.around(np.mean(Ub_RMSDs[t]), 5),
                        np.around(np.mean(BIASs[t]), 5),
                        np.around(np.mean(STDs_an[t]), 5),
                        np.around(np.mean(STDs_obs[t]), 5),
                    ]
                )
            ],
        )

    table_vals = table_vals.reshape((len(timespan), 7)).transpose()

    row_labels = [
        "R (-)",
        "R anom (-)",
        "RMSD (" + unit + ")",
        "Ub RMSD (" + unit + ")",
        "Bias (" + unit + ")",
        "Std an (" + unit + ")",
        "Std obs (" + unit + ")",
    ]

    if table_format == "PDF":
        plt.table(
            cellText=table_vals,
            rowLabels=row_labels,
            colLabels=col_labels,
            loc="center right",
        )
        plt.tight_layout()
        fig.tight_layout()
        fig.savefig(
            table_dir
            + "/"
            + station_type
            + "_"
            + layer
            + "_"
            + str(start_year)
            + "_"
            + str(end_year)
            + "_"
            + times
            + "_"
            + land_type
            + "_"
            + plot_type
            + ".pdf",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()
    elif table_format == "ASCII":
        row_labels_ascii = [
            "R (-),",
            "R anom (-),",
            "RMSD (" + unit + "),",
            "Ub RMSD (" + unit + "),",
            "Bias (" + unit + "),",
            "Std an (" + unit + "),",
            "Std obs (" + unit + "),",
        ]

        for jj in range(len(row_labels_ascii[:])):
            row_labels_ascii[jj] = row_labels_ascii[jj].ljust(24)

        col_labels_ascii = [exp_name + ","]
        for jj in range(len(col_labels[:])):
            col_labels_ascii.append(col_labels[jj] + ",")
        for jj in range(len(col_labels[:])):
            col_labels_ascii[jj] = col_labels_ascii[jj].ljust(24)

        col_labels_ascii = "".join(col_labels_ascii)

        filename = (
            table_dir
            + "/"
            + station_type
            + "_"
            + layer
            + "_"
            + str(start_year)
            + "_"
            + str(end_year)
            + "_"
            + times
            + "_"
            + land_type
            + "_"
            + plot_type
            + ".txt"
        )

        f = open(
            filename,
            "w",
        )
        filename_html = "file://" + filename
        f.writelines(col_labels_ascii + "\n")
        for jj in range(len(row_labels_ascii[:])):
            table_vals_ascii = []
            for vals in table_vals[jj]:
                table_vals_ascii.append(str(vals) + ",")
            for kk in range(len(table_vals[jj])):
                table_vals_ascii[kk] = table_vals_ascii[kk].ljust(24)
            table_vals_ascii = "".join(table_vals_ascii)
            f.writelines(row_labels_ascii[jj] + table_vals_ascii + "\n")
        f.close()

        dic = {
            "expver": [exp_name],
            "date": [start_year + " to " + end_year],
            "type": ["Table"],
            "network": [station_type],
            "variable": [layer + " " + var_type],
            "metric": [plot_type + " scores "],
            "link": filename_html,
        }

        df = df.append(pd.DataFrame(dic), ignore_index=True, sort=False)

    return df


def draw_station_box_plots_period(
    CC, ACC, w, plots_dir, Networks, Layer_types, year_1, year_2
):

    #    ## BoxPlots

    fig, axes = plt.subplots(nrows=len(Networks), ncols=len(Layer_types))

    fig.set_figheight(len(Networks) * 7.5)
    fig.set_figwidth(14)

    alphabet = ["a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"]

    # Plot monthly stats---------------------------------------------------------------------------------------
    ii = 0
    for n, network in enumerate(Networks):

        for l, layer_type in enumerate(Layer_types):

            if layer_type == "Surface":
                data_to_plot = [
                    np.array(CC[network + "_" + layer_type]["period"])[
                        w[network + "_" + layer_type]["period"]
                    ],
                    np.array(ACC[network + "_" + layer_type]["period"])[
                        w[network + "_" + layer_type]["period"]
                    ],
                ]
            else:
                data_to_plot = np.array(CC[network + "_" + layer_type]["period"])[
                    w[network + "_" + layer_type]["period"]
                ]

            bp = axes[n, l].boxplot(data_to_plot, patch_artist=True)
            axes[n, l].axis([0.5, 2.5, 0.0, 1], fontsize=20)
            axes[n, l].set_title(
                alphabet[ii] + " Pearson CC " + network + " (" + str(layer_type) + ")",
                fontsize=20,
            )
            for box in bp["boxes"]:
                # change outline color
                box.set(color="#7570b3", linewidth=1.5)
            # change fill color
            #    box.set( facecolor = '#1b9e77' )
            # change color and linewidth of the whiskers
            for whisker in bp["whiskers"]:
                whisker.set(color="#7570b3", linewidth=2)
                # change color and linewidth of the caps
            for cap in bp["caps"]:
                cap.set(color="#7570b3", linewidth=2)
            # change color and linewidth of the medians
            for median in bp["medians"]:
                median.set(color="#b2df8a", linewidth=2)
            ## change the style of fliers and their fill
            for flier in bp["fliers"]:
                flier.set(marker="o", color="#e7298a", alpha=0.25)
            if layer_type == "Surface":
                axes[n, l].set_xticklabels(["CC", "CC_ano"], fontsize=20)
            else:
                axes[n, l].set_xticklabels(["CC"], fontsize=20)
            ## Remove top axes and right axes ticks
            axes[n, l].get_xaxis().tick_bottom()
            axes[n, l].get_yaxis().tick_left()
            axes[n, l].tick_params(axis="y", labelsize=20)
            axes[n, l].tick_params(axis="x", labelsize=20)
            ii = ii + 1
    plt.tick_params(labelsize=20)
    plt.savefig(
        plots_dir
        + "/box_plot_Pearson_R_"
        + network
        + "_"
        + str(year_1)
        + "_"
        + str(year_2)
        + "_period.pdf",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


def draw_station_map(
    df,
    start_period,
    end_period,
    lat_networks,
    lon_networks,
    data_dict,
    layers,
    Station_Networks,
    plots_dir,
    var_type,
    EXP,
    times,
    land_type,
    sig_comparisons=False,
):

    if sig_comparisons:
        vmin, vmax = -1.5, 2.5
        cm = plt.cm.get_cmap("jet", 4.0)
    else:
        vmin, vmax = 0.0, 1.0
        cm = plt.cm.get_cmap("viridis")

    map_labels = dict()
    countries = ["US", "France", "Aus", "Spain", "Germany"]
    for country in countries:
        map_labels[country] = list()
    for network in Station_Networks:
        if network in ["USCRN", "SCAN", "SNOTEL"]:
            map_labels["US"].append(network)
        if network in ["SMOSMANIA"]:
            map_labels["France"].append(network)
        if network in ["REMEDHUS"]:
            map_labels["Spain"].append(network)
        if network in ["OZNET"]:
            map_labels["Aus"].append(network)
        if network in ["TERENO"]:
            map_labels["Germany"].append(network)

    for country in countries:
        if len(map_labels[country]) == 0:
            continue

        for layer in layers:
            plt.figure(figsize=(12, 6))
            plt.subplot(1, 1, 1)

            if country == "US":
                ax = plt.axes(projection=ccrs.Robinson())
                markers = ("o", "^", "s")
                extent = [
                    -130.0,
                    -60.0,
                    25.0,
                    75.0,
                ]  # Define the extent (min lon, max lon, min lat, max lat)
                ax.set_extent(extent)
                ax.set_extent(extent)
                ax.coastlines()
                ax.add_feature(cart.feature.BORDERS)

            elif country == "France":
                ax = plt.axes(projection=ccrs.Robinson())
                extent = [
                    -4.1,
                    11,
                    41.3,
                    50.9,
                ]  # Define the extent (min lon, max lon, min lat, max lat)
                ax.set_extent(extent)
                ax.set_extent(extent)
                ax.coastlines()
                ax.add_feature(cart.feature.BORDERS)

            elif country == "Aus":
                ax = plt.axes(projection=ccrs.Robinson())
                extent = [
                    112.0,
                    154.0,
                    -44.0,
                    -9.0,
                ]  # Define the extent (min lon, max lon, min lat, max lat)
                ax.set_extent(extent)
                ax.set_extent(extent)
                ax.coastlines()
                ax.add_feature(cart.feature.BORDERS)

            elif country == "Spain":

                ax = plt.axes(projection=ccrs.Robinson())
                extent = [
                    -10.0,
                    4.0,
                    36.0,
                    44.0,
                ]  # Define the extent (min lon, max lon, min lat, max lat)
                ax.set_extent(extent)
                ax.set_extent(extent)
                ax.coastlines()
                ax.add_feature(cart.feature.BORDERS)

            elif country == "Germany":

                ax = plt.axes(projection=ccrs.Robinson())
                extent = [
                    6.0,
                    15.0,
                    47.0,
                    55.0,
                ]  # Define the extent (min lon, max lon, min lat, max lat)
                ax.set_extent(extent)
                ax.set_extent(extent)
                ax.coastlines()
                ax.add_feature(cart.feature.BORDERS)

            n = 0
            markers = ("o", "^", "s")
            for network in map_labels[country]:
                if len(lon_networks[network + "_" + layer]["period"]) == 0:
                    continue

                lon_scatter = np.array(
                    lon_networks[network + "_" + layer]["period"]
                ).astype("float")
                lat_scatter = np.array(
                    lat_networks[network + "_" + layer]["period"]
                ).astype("float")

                if country == "US" or country == "France" or country == "Germany":
                    scatter_size = 20
                elif country == "Spain" or country == "Aus":
                    scatter_size = 5

                im = ax.scatter(
                    x=lon_scatter,
                    y=lat_scatter,
                    marker=markers[n],
                    c=data_dict[network + "_" + layer]["period"],
                    vmin=vmin,
                    vmax=vmax,
                    cmap=cm,
                    transform=ccrs.PlateCarree(),
                    s=scatter_size,
                    zorder=2,
                    label=network,
                )

                n = n + 1

            if sig_comparisons:
                map_dir = (
                    plots_dir
                    + country
                    + "_map_Peason_R_vs_control_layer"
                    + layer
                    + "_"
                    + times
                    + "_"
                    + land_type
                    + ".png"
                )
                title = "Pearson R of " + EXP + " vs control in " + country
            else:
                map_dir = (
                    plots_dir
                    + country
                    + "_map_Pearson_R_layer_"
                    + layer
                    + "_"
                    + times
                    + "_"
                    + land_type
                    + ".png"
                )
                title = "Pearson R in " + country

            filename_html = "file://" + map_dir

            if sig_comparisons:
                expver = EXP + " vs control"
            else:
                expver = EXP

            dic = {
                "expver": [expver],
                "date": [start_period + " to " + end_period],
                "type": ["Map"],
                "network": ["inc"],
                "variable": [layer + " " + var_type],
                "metric": [country],
                "link": filename_html,
            }

            df = df.append(pd.DataFrame(dic), ignore_index=True, sort=False)

            if len(lon_networks[network + "_" + layer]["period"]) > 0:

                if sig_comparisons:
                    cbar = plt.gcf().colorbar(
                        im,
                        ax=ax,
                        ticks=np.arange(-1, 3),
                        orientation="horizontal",
                        fraction=0.04,
                    )
                    cbar.set_ticklabels(
                        ["Sig deg", "Non-sig deg", "Non-sig imp", "Sig imp"]
                    )
                else:
                    plt.gcf().colorbar(im, ax=ax, orientation="vertical", fraction=0.04)

                ax.legend()
                plt.title(title, fontsize=14)
                plt.savefig(map_dir, dpi=300, bbox_inches="tight")
                plt.close()
            else:
                print("No stations available for map")

    return df


def error_bar_comparison(
    df,
    EXPVER,
    NET,
    ACC,
    w,
    ACC_conf_upper,
    ACC_conf_lower,
    plot_dir,
    table_dir,
    layers,
    time_period,
    start_period,
    end_period,
    var_type,
    times,
    land_type,
):

    for layer_type in layers:
        # data to plot
        n_groups = len(NET)

        network_means = dict()
        Upper_conf = dict()
        Lower_conf = dict()
        rects = dict()
        network_labels = list()

        for EXP in EXPVER:
            rects[EXP] = list()
            network_means[EXP] = list()
            Upper_conf[EXP] = list()
            Lower_conf[EXP] = list()
            for n, network in enumerate(NET):
                network_means[EXP].append(
                    np.nanmean(
                        np.array(ACC[EXP][network + "_" + layer_type][time_period])[
                            w[EXP][network + "_" + layer_type][time_period]
                        ]
                    )
                )
                Upper_conf[EXP].append(
                    np.nanmean(
                        np.array(
                            ACC_conf_upper[EXP][network + "_" + layer_type][time_period]
                        )[w[EXP][network + "_" + layer_type][time_period]]
                    )
                )
                Lower_conf[EXP].append(
                    np.nanmean(
                        np.array(
                            ACC_conf_lower[EXP][network + "_" + layer_type][time_period]
                        )[w[EXP][network + "_" + layer_type][time_period]]
                    )
                )
                # if (EXP==EXPVER[0]):
                #  network_labels[n].append(network)
        # create plot
        fig, ax = plt.subplots()
        index = np.arange(n_groups)
        bar_width = 1.0 / (len(EXPVER))

        for exp, EXP in enumerate(EXPVER):

            # if (EXP=='0001'):
            #   label1='ERA-5'
            # elif (EXP=='guji'):
            #   label1='H27'
            # elif (EXP=='h5zn'):
            #   label1='H141'
            rects[EXP] = plt.errorbar(
                index + (bar_width * 0.5 * (exp + 1)),
                network_means[EXP],
                yerr=[
                    np.array(network_means[EXP]) - np.array(Lower_conf[EXP]),
                    np.array(Upper_conf[EXP]) - np.array(network_means[EXP]),
                ],
                fmt="o",
                elinewidth=2.0,
                label=EXP,
            )

        plt.xlabel("Network")
        plt.ylabel("R anomaly (-)")
        plt.title("R anomaly by network")
        plt.xlim(0.0, index[-1] + bar_width * (0.5 * (len(EXPVER) + 1)))
        plt.xticks(
            index + 0.5 * ((bar_width * 0.5) + (bar_width * 0.5 * (len(EXPVER)))), NET
        )
        plt.legend(loc=0)

        filename = (
            plot_dir
            + "/ACC_vs_all_exps_"
            + layer_type
            + "_"
            + times
            + "_"
            + land_type
            + ".png"
        )

        plt.savefig(
            filename,
            dpi=300,
            bbox_inches="tight",
        )
        plt.tight_layout()
        plt.close()

        filename_html = (
            "file://"
            + plot_dir
            + "/ACC_vs_all_exps_"
            + layer_type
            + "_"
            + times
            + "_"
            + land_type
            + ".png"
        )

        dic = {
            "expver": ["All"],
            "date": [start_period + " to " + end_period],
            "type": ["CI"],
            "network": ["inc"],
            "variable": [layer_type + " " + var_type],
            "metric": ["Pearson R anom (-)"],
            "link": filename_html,
        }

        #        dic[
        #            "link"
        #        ] += f'<a href="{filename2}">'

        df = df.append(pd.DataFrame(dic), ignore_index=True, sort=False)

        NET_2 = copy.deepcopy(NET)

        # Draw ASCII table of results

        ACC_combined = dict()
        Lower_conf_combined = dict()
        Upper_conf_combined = dict()

        row_labels = list()
        network_labels = list()
        NET_2.append("All")
        for exp, EXP in enumerate(EXPVER):
            ACC_combined[EXP] = list()
            Lower_conf_combined[EXP] = list()
            Upper_conf_combined[EXP] = list()

            for n, network in enumerate(NET_2):
                row_labels.append(EXP)
                network_labels.append(network)

                if network == "All":
                    network_means[EXP].append(np.nanmean(np.array(ACC_combined[EXP])))
                    Lower_conf[EXP].append(
                        np.nanmean(np.array(Lower_conf_combined[EXP]))
                    )
                    Upper_conf[EXP].append(
                        np.nanmean(np.array(Upper_conf_combined[EXP]))
                    )
                else:
                    ACC_combined[EXP].extend(
                        np.array(ACC[EXP][network + "_" + layer_type][time_period])[
                            w[EXP][network + "_" + layer_type][time_period]
                        ]
                    )
                    Lower_conf_combined[EXP].extend(
                        np.array(
                            ACC_conf_lower[EXP][network + "_" + layer_type][time_period]
                        )[w[EXP][network + "_" + layer_type][time_period]]
                    )
                    Upper_conf_combined[EXP].extend(
                        np.array(
                            ACC_conf_upper[EXP][network + "_" + layer_type][time_period]
                        )[w[EXP][network + "_" + layer_type][time_period]]
                    )

        for jj in range(len(row_labels[:])):
            row_labels[jj] = row_labels[jj].ljust(24)
            network_labels[jj] = network_labels[jj].ljust(24)

        col_labels = [
            "Exp",
            "Network",
            "Pearson R anom (-),",
            "Lower conf (95%),",
            "Upper conf (95%),",
        ]
        for jj in range(len(col_labels[:])):
            col_labels[jj] = col_labels[jj].ljust(24)

        f = open(
            table_dir
            + "/ACC_sig_scores_All_exps_"
            + layer_type
            + "_"
            + time_period
            + "_"
            + times
            + "_"
            + land_type
            + ".txt",
            "w",
        )
        col_labels = "".join(col_labels)
        f.writelines(col_labels + "\n")

        for exp, EXP in enumerate(EXPVER):
            for n, network in enumerate(NET_2):
                table_vals = [
                    EXP,
                    network,
                    str(network_means[EXP][n]),
                    str(Lower_conf[EXP][n]),
                    str(Upper_conf[EXP][n]),
                ]
                for kk in range(len(table_vals)):
                    table_vals[kk] = table_vals[kk].ljust(24)
                table_vals = "".join(table_vals)
                if (n == 0) and (exp == 0):
                    print(layer_type)
                    print("##################################################")
                    print(col_labels)
                print(table_vals)
                f.writelines(table_vals + "\n")
        f.close()

    return df


def bar_chart_comparison(
    df,
    start_period,
    end_period,
    EXPVER,
    NET,
    Ano_R_networks,
    plot_dir,
    layers,
    time_period,
    times,
    land_type,
    var_type,
):

    # data to plot
    n_groups = 4

    network_sum = dict()
    network_sum_sig = dict()
    Upper_conf = dict()
    Lower_conf = dict()
    rects = dict()
    colours = ["b", "r", "g", "k", "c", "m"]

    for layer_type in layers:

        for EXP in EXPVER[1:]:
            rects[EXP] = list()
            network_sum[EXP] = list()
            network_sum_sig[EXP] = list()
            Upper_conf[EXP] = list()
            Lower_conf[EXP] = list()
            for n, network in enumerate(NET):

                network_sum[EXP].extend(
                    Ano_R_networks[EXP][network + "_" + layer_type][
                        time_period
                    ].tolist()
                )

            network_sum_sig[EXP] = [
                np.sum(np.array(network_sum[EXP]) == -1),
                np.sum(np.array(network_sum[EXP]) == 0),
                np.sum(np.array(network_sum[EXP]) == 1),
                np.sum(np.array(network_sum[EXP]) == 2),
            ]
            # if (EXP==EXPVER[0]):
            #  network_labels[n].append(network)

        # create plot
        fig, ax = plt.subplots()
        index = np.arange(n_groups)
        bar_width = 0.35
        opacity = 1.0

        for exp, EXP in enumerate(EXPVER[1:]):

            rects[EXP] = plt.bar(
                index + (bar_width * (exp)),
                network_sum_sig[EXP],
                bar_width,
                alpha=opacity,
                color=colours[exp],
                label=EXP,
            )

            plt.text(
                0 + (bar_width * (exp)),
                network_sum_sig[EXP][0],
                str(
                    np.round(
                        network_sum_sig[EXP][0]
                        / len(np.array(network_sum[EXP]))
                        * 100.0,
                        1,
                    )
                )
                + "%",
            )
            plt.text(
                1 + (bar_width * (exp)),
                network_sum_sig[EXP][1],
                str(
                    np.round(
                        network_sum_sig[EXP][1]
                        / len(np.array(network_sum[EXP]))
                        * 100.0,
                        1,
                    )
                )
                + "%",
            )
            plt.text(
                2 + (bar_width * (exp)),
                network_sum_sig[EXP][2],
                str(
                    np.round(
                        network_sum_sig[EXP][2]
                        / len(np.array(network_sum[EXP]))
                        * 100.0,
                        1,
                    )
                )
                + "%",
            )
            plt.text(
                3 + (bar_width * (exp)),
                network_sum_sig[EXP][3],
                str(
                    np.round(
                        network_sum_sig[EXP][3]
                        / len(np.array(network_sum[EXP]))
                        * 100.0,
                        1,
                    )
                )
                + "%",
            )

        plt.xlabel("Pearson R anomaly of exp vs control")
        plt.ylabel("Number of stations")
        plt.title("Significance of correlation differences")
        plt.xticks(
            index + 0.5 * ((bar_width * 0.5) + (bar_width * 0.5 * (len(EXPVER[1:])))),
            ("Sig deg", "Non-sig deg", "Non-sig imp", "Sig imp"),
        )
        plt.legend()
        plt.tight_layout()
        plt.savefig(
            plot_dir
            + "/ACC_vs_control_bar_chart_"
            + layer_type
            + "_"
            + times
            + "_"
            + land_type
            + ".png",
            dpi=300,
            bbox_inches="tight",
        )

        plt.tight_layout()
        plt.close()

        filename_html = (
            "file://"
            + plot_dir
            + "/ACC_vs_control_bar_chart_"
            + layer_type
            + "_"
            + times
            + "_"
            + land_type
            + ".png"
        )

        dic = {
            "expver": ["All"],
            "date": [start_period + " to " + end_period],
            "type": ["CI"],
            "network": ["inc"],
            "variable": [layer_type + " " + var_type],
            "metric": ["Bar chart R anom (-)"],
            "link": filename_html,
        }

        #        dic[
        #            "link"
        #        ] += f'<a href="{filename2}">'

        df = df.append(pd.DataFrame(dic), ignore_index=True, sort=False)

    return df


def draw_box_plot_period(
    CC, ACC, w, plot_dir, network, layer_type, year_1, year_2, times, land_type
):

    #    ## BoxPlots
    fig = plt.figure()
    ax = fig.add_subplot(111)
    data_to_plot = [np.array(CC)[w], np.array(ACC)[w]]
    #
    bp = ax.boxplot(data_to_plot, patch_artist=True)
    ax.axis([0.5, 2.5, 0.0, 1])
    plt.title(
        "Pearson R " + network + " vs. SM-DAS-2 (" + str(layer_type) + ")", fontsize=14
    )
    for box in bp["boxes"]:
        # change outline color
        box.set(color="#7570b3", linewidth=1.5)
    # change fill color
    #    box.set( facecolor = '#1b9e77' )
    # change color and linewidth of the whiskers
    for whisker in bp["whiskers"]:
        whisker.set(color="#7570b3", linewidth=2)
        # change color and linewidth of the caps
    for cap in bp["caps"]:
        cap.set(color="#7570b3", linewidth=2)
    # change color and linewidth of the medians
    for median in bp["medians"]:
        median.set(color="#b2df8a", linewidth=2)
    ## change the style of fliers and their fill
    for flier in bp["fliers"]:
        flier.set(marker="o", color="#e7298a", alpha=0.25)
    ax.set_xticklabels(["R", "R_ano"], fontsize=15)
    ## Remove top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.tick_params(labelsize=15)
    plt.savefig(
        plot_dir
        + "/box_plot_Pearson_R_"
        + network
        + "_"
        + layer_type
        + "_"
        + str(year_1)
        + "_"
        + str(year_2)
        + "_"
        + times
        + "_"
        + land_type
        + "_period.pdf",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


def draw_station_box_plots_sub_periods(
    var_type,
    expver,
    scores,
    w,
    plot_dir,
    network,
    NET,
    layer_type,
    sub_periods,
    years,
    plot_type,
    SM_units,
    start_period,
    end_period,
    times,
    land_type,
    df,
):

    #    ## BoxPlots

    if var_type == "SM":
        if SM_units == "Volumetric":
            unit = "m3/m3"
        else:
            unit = "-"
    else:
        unit = "K"

    var_types = ["CC", "ACC", "RMSD", "Ub_RMSD", "Bias"]
    var_labels = [
        "R (-)",
        "R anom (-)",
        "RMSD (" + unit + ")",
        "Ub RMSD (" + unit + ")",
        "Bias (" + unit + ")",
    ]

    var_labels_no_units = ["R", "R anom", "RMSD", "UbRMSD", "Bias"]

    # Plot monthly stats---------------------------------------------------------------------------------------
    ii = 0

    for v, var in enumerate(scores):

        fig, axes = plt.subplots(nrows=1, ncols=1)

        fig.set_figheight(7.5)
        fig.set_figwidth(15)

        data_to_plot = list()
        col_labels = list()

        if network == "All_stations":
            var["All_stations" + "_" + layer_type] = dict()
            for p, period in enumerate(sub_periods):
                var["All_stations" + "_" + layer_type][str(period)] = list()
                for n, nets in enumerate(NET):
                    var["All_stations" + "_" + layer_type][str(period)].extend(
                        np.array(var[nets + "_" + layer_type][str(period)])[
                            w[nets + "_" + layer_type][str(period)]
                        ]
                    )
                if len(np.array(var[network + "_" + layer_type][str(period)])) > 1:
                    data_to_plot.append(
                        np.array(var[network + "_" + layer_type][str(period)])
                    )
                    col_labels.append(
                        str(period)
                        + "\n ("
                        + str(len(var[network + "_" + layer_type][str(period)][:]))
                        + ")"
                    )

        else:
            for p, period in enumerate(sub_periods):
                if (
                    len(
                        np.array(var[network + "_" + layer_type][str(period)])[
                            w[network + "_" + layer_type][str(period)]
                        ]
                    )
                    > 1
                ):
                    data_to_plot.append(
                        np.array(var[network + "_" + layer_type][str(period)])[
                            w[network + "_" + layer_type][str(period)]
                        ]
                    )
                    col_labels.append(
                        str(period)
                        + "\n ("
                        + str(len(w[network + "_" + layer_type][str(period)][:]))
                        + ")"
                    )

        if len(data_to_plot) > 0:

            data_to_plot_df = pd.DataFrame(data_to_plot).transpose()
            data_to_plot_df.columns = col_labels[:]
            data_to_plot_df = data_to_plot_df.assign(expver=expver)
            df_scores = data_to_plot_df

            df_scores_long_form = df_scores.melt(
                id_vars=["expver"],
                var_name="Year (n stations)",
                value_name=var_labels[v],
            )
            sns.boxplot(
                data=df_scores_long_form,
                x="Year (n stations)",
                y=var_labels[v],
                hue="expver",
                showfliers=False,
                ax=axes,
            )
            sns.pointplot(
                data=df_scores_long_form,
                x="Year (n stations)",
                y=var_labels[v],
                hue="expver",
                ci=None,
                ax=axes,
            )
            axes.yaxis.grid(True)

            plt.xlabel("Year (n stations)", fontsize=16)
            plt.ylabel(var_labels[v], fontsize=16)
            plt.tick_params(axis="both", which="major", labelsize=14)

            ii = ii + 1
            plt.tick_params(labelsize=20)
            filename = (
                plot_dir
                + "/box_plot_"
                + var_types[v]
                + "_"
                + network
                + "_"
                + layer_type
                + "_"
                + str(years[0])
                + "_"
                + str(years[-1])
                + "_"
                + times
                + "_"
                + land_type
                + "_"
                + plot_type
                + ".pdf"
            )
            plt.savefig(
                filename,
                dpi=300,
                bbox_inches="tight",
            )
            plt.close()

            #            if (var_types[v] == "CC") or (var_types[v] == "ACC"):

            filename_html = "file://" + filename

            dic = {
                "expver": [expver],
                "date": [start_period + " to " + end_period],
                "type": ["Box plot"],
                "network": [network],
                "variable": [layer_type + " " + var_type],
                "metric": [plot_type + " " + var_labels_no_units[v]],
                "link": filename_html,
            }

            df = df.append(pd.DataFrame(dic), ignore_index=True, sort=False)

    return df


def draw_station_box_plots_sub_periods_allExp(
    expvers,
    scoresAll,
    wAll,
    plot_dir,
    NET,
    layers,
    sub_periods,
    years,
    plot_type,
    SM_units,
    expl,
    var_type,
    start_period,
    end_period,
    times,
    land_type,
    df,
):

    #    ## BoxPlots

    if var_type == "SM":
        if SM_units == "Volumetric":
            unit = "m$^3$/m$^3$"
        else:
            unit = "-"
    else:
        unit = "K"

    var_types = ["CC", "ACC", "RMSD", "Ub_RMSD", "Bias"]
    var_labels = [
        "R (-)",
        "R anom (-)",
        "RMSD (" + unit + ")",
        "Ub RMSD (" + unit + ")",
        "Bias (" + unit + ")",
    ]
    var_labels_no_units = ["R", "R anom", "RMSD", "UbRMSD", "Bias"]
    colors = ["k", "r", "g", "b", "m", "y", "c"]
    expc = dict()
    for e, EXP in enumerate(expvers):
        expc[EXP] = colors[e]

    # Plot monthly stats---------------------------------------------------------------------------------------
    ii = 0
    dx = 0.15

    nexps = len(expvers)

    NETWORKS = copy.deepcopy(NET)
    NETWORKS.append("All_stations")

    for network in NETWORKS:

        for layer_type in layers:

            for v, cvars in enumerate(scoresAll):

                fig, axes = plt.subplots(nrows=1, ncols=1)

                fig.set_figheight(7.5)
                fig.set_figwidth(15)

                wp = dict()
                var = dict()
                bps = list()

                for e, EXP in enumerate(expvers):
                    xp = []
                    data_to_plot = list()
                    #            print(EXP)
                    col_labels = list()

                    wp = wAll[EXP]
                    var = cvars[EXP]

                    if network == "All_stations":
                        var["All_stations" + "_" + layer_type] = dict()
                        for p, period in enumerate(sub_periods):
                            var["All_stations" + "_" + layer_type][str(period)] = list()
                            for n, nets in enumerate(NET):

                                try:
                                    var["All_stations" + "_" + layer_type][
                                        str(period)
                                    ].extend(
                                        np.array(
                                            var[nets + "_" + layer_type][str(period)]
                                        )[wp[nets + "_" + layer_type][str(period)]]
                                    )
                                except:
                                    import pdb

                                    pdb.set_trace()
                            if (
                                len(
                                    np.array(
                                        var[network + "_" + layer_type][str(period)]
                                    )
                                )
                                > 1
                            ):
                                data_to_plot.append(
                                    np.array(
                                        var[network + "_" + layer_type][str(period)]
                                    )
                                )
                                col_labels.append(
                                    str(period)
                                    + "\n ("
                                    + str(
                                        len(
                                            var[network + "_" + layer_type][
                                                str(period)
                                            ][:]
                                        )
                                    )
                                    + ")"
                                )
                                xp.append(p * 1 + e * dx)

                    else:
                        for p, period in enumerate(sub_periods):

                            try:

                                if (
                                    len(
                                        np.array(
                                            var[network + "_" + layer_type][str(period)]
                                        )[wp[network + "_" + layer_type][str(period)]]
                                    )
                                    > 1
                                ):
                                    data_to_plot.append(
                                        np.array(
                                            var[network + "_" + layer_type][str(period)]
                                        )[wp[network + "_" + layer_type][str(period)]]
                                    )

                                    col_labels.append(
                                        str(period)
                                        + "\n ("
                                        + str(
                                            len(
                                                wp[network + "_" + layer_type][
                                                    str(period)
                                                ][:]
                                            )
                                        )
                                        + ")"
                                    )

                                    xp.append(p * 1 + e * dx)

                            except:
                                import pdb

                                pdb.set_trace()

                    if len(data_to_plot) > 0:
                        #                print('PRINTING ALL DATA')
                        #                print(len(data_to_plot), len(xp))

                        data_to_plot_df = pd.DataFrame(data_to_plot).transpose()
                        data_to_plot_df.columns = col_labels[:]
                        data_to_plot_df = data_to_plot_df.assign(expver=EXP)
                        if e == 0:
                            df_scores = data_to_plot_df
                        else:
                            df_scores = pd.concat([df_scores, data_to_plot_df])

                        ii = ii + 1
                        # store boxplot for legend

                try:

                    df_scores_long_form = df_scores.melt(
                        id_vars=["expver"],
                        var_name="Year (n stations)",
                        value_name=var_labels[v],
                    )
                    sns.boxplot(
                        data=df_scores_long_form,
                        x="Year (n stations)",
                        y=var_labels[v],
                        hue="expver",
                        showfliers=False,
                        ax=axes,
                    )
                    sns.pointplot(
                        data=df_scores_long_form,
                        x="Year (n stations)",
                        y=var_labels[v],
                        hue="expver",
                        ci=None,
                        ax=axes,
                    )
                    axes.yaxis.grid(True)
                    #                axes.legend([b for b in bps[:]], [l for l in expl], ncol=2, loc="lower right")

                    plt.xlabel("Year (n stations)", fontsize=16)
                    plt.ylabel(var_labels[v], fontsize=16)
                    plt.tick_params(axis="both", which="major", labelsize=14)

                    filename = (
                        plot_dir
                        + "/box_plot_all_"
                        + var_types[v]
                        + "_"
                        + network
                        + "_"
                        + layer_type
                        + "_"
                        + str(years[0])
                        + "_"
                        + str(years[-1])
                        + "_"
                        + times
                        + "_"
                        + land_type
                        + "_"
                        + plot_type
                        + ".pdf"
                    )

                    plt.savefig(
                        filename,
                        dpi=300,
                        bbox_inches="tight",
                    )
                    plt.close()

                    filename_html = "file://" + filename

                    dic = {
                        "expver": ["All"],
                        "date": [start_period + " to " + end_period],
                        "type": ["Box plot"],
                        "network": [network],
                        "variable": [layer_type + " " + var_type],
                        "metric": [plot_type + " " + var_labels_no_units[v]],
                        "link": filename_html,
                    }

                    df = df.append(pd.DataFrame(dic), ignore_index=True, sort=False)

                except:
                    continue

    return df
