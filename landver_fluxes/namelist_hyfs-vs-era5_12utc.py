#!/usr/bin/env python2
# -*- coding: utf-8 -*-

######################################################################
# --- example namelist ---
# validation of hyfs vs era5
# 12 utc only (with time conversion)
# whole year 2020
# accumulation of era5 on 6h as well (to make it comparable with hyfs)
######################################################################
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
TIME = ["12","12"]
STEP =      ["00","00"]
DOMAIN =    ["G","G"]
GRID =      ["av","av"]
EXPNAME =  ["ERA5","hyfs"]

#Switches
LT2UTC = True #LT to UTC?
UTC2LT = False #UTC to LT? (both False means, that no conversion will be applies)
ACC6h = True #Accumulate on 6h?

