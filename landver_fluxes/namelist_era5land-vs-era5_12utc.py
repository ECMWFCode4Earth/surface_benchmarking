#!/usr/bin/env python2
# -*- coding: utf-8 -*-

######################################################################
# --- example namelist ---
# validation of era5land vs era5
# 12 utc (test run)
# oct-dec 2020
# no accumulation, use instant fluxes
######################################################################
#-------------------------------------------------------------------------------------
analysis_period = ['2020-10-01','2020-12-31'] #Analysis period.

#Output data directory where all analysis data and results are to be stored: 
#output_dir="/Users/dadf/landver_test_21_04_23/sfc_evaluation_out"#"/perm/rd/dadf/H14_vs_H26"
output_dir="/home/lauracma/Documents/ecmwf_proj/data/model_fluxes/" #adapted

#Experimental parameters
EXPVER =    ["0001","0001"] #Experiment names needed for extraction, preprocessing and validation 
#First experiment is considered the control 
CLASS =     ["ea","l5"]
TYPE =      ["AN","AN"]
ANOFFSET =  ["9","9"]
STREAM =    ["oper","oper"]
REPRES =    ["SH","SH"]
LEVTYPE =   ["SFC","SFC"]
TIME = ["12","12"]
STEP =      ["00","00"]
DOMAIN =    ["G","G"]
GRID =      ["av","av"]
EXPNAME =  ["ERA5","ERA5Land"]

#Switches
LT2UTC = True #LT to UTC? (shifts insitu (!) data)
UTC2LT = False #UTC to LT? (both False means, that no conversion will be applied)
ACC6h = False #Accumulate on 6h? (accumulates model/obs data)

#REMARKS:
#- used 3 months only
#- the provided era5land files had the wrong names (2022 instead of 2020) and sshf is not accumulated whereas slhf is