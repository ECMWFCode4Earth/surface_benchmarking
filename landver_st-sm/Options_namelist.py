#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 14:10:37 2017

@author: dadf
"""

################general options: 
#-------------------------------------------------------------------------------------
analysis_period = ['2020-01-01','2020-12-31'] #Analysis period.

#Output data directory where all analysis data and results are to be stored: 
#output_dir="/Users/dadf/landver_test_21_04_23/sfc_evaluation_out"#"/perm/rd/dadf/H14_vs_H26"
output_dir="/home/lauracma/Documents/ecmwf_proj/data/sfc_evaluation_out/" #adapted


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
