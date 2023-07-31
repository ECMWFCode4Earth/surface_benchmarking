#!/usr/bin/env python2
# -*- coding: utf-8 -*-

################general options: 
#-------------------------------------------------------------------------------------
analysis_period = ['2020-01-01','2020-12-31'] #Analysis period.

#Output data directory where all analysis data and results are to be stored: 
#output_dir="/Users/dadf/landver_test_21_04_23/sfc_evaluation_out"#"/perm/rd/dadf/H14_vs_H26"
output_dir="/home/lauracma/Documents/ecmwf_proj/data/model_fluxes/" #adapted

#Experimental parameters
EXPVER =    ["0001","hyfs"] #Experiment names needed for extraction, preprocessing and validation 
#First experiment is considered the control 
CLASS =     ["ea","rd"]
TYPE =      ["AN","AN"]
ANOFFSET =  ["9","9"]
STREAM =    ["oper","oper"]
REPRES =    ["SH","SH"]
LEVTYPE =   ["SFC","SFC"]
TIME      = ["00/12","00/12"] 
STEP =      ["00","00"]
DOMAIN =    ["G","G"]
GRID =      ["av","av"]
EXPNAME =  ["ERA5","hyfs"]

"""
notes:
- SM, ST files had the format yyyymmddhh but the fluxes files have yyyymm, so we need to extract the rest from the file itself, but probably it's easier to use the resolution than
- era5 and hyfs contain differently accumulated data, so they have to be processed accordingly (in the loop over EXPVER, you have to know which experiment)
- add a namelist parameter for resolution? or still use time?
-> all in LANDVER_pre.py
"""