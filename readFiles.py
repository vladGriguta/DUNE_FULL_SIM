#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 21:23:12 2018

@author: vladgriguta
"""
import pandas as pd
import numpy as np

def Read_SN(location):
    # Read SN events decay time
    time_SN_vis = pd.read_csv(location+'/time_SN_vis.txt')
    time_SN_vuv = pd.read_csv(location+'/time_SN_vuv.txt')
    
    # Concatenate all into a new DataFrame object
    time_SN_vuv['vuv'] = 1
    time_SN_vis['vuv'] = 0
    frames = [time_SN_vuv,time_SN_vis]
    timeSN = pd.concat(frames)
    timeSN = timeSN.sort_values('time')
    timeSN = timeSN.reset_index(drop=True)
    del time_SN_vis, time_SN_vuv
    
    nr_events_SN = int(np.max(timeSN['event']))+1
    
    return timeSN,nr_events_SN

def Read_Ar39(location):
    # Read Ar39 events decay time
    time_Ar_vis = pd.read_csv(location+'/time_Ar_vis.txt')
    time_Ar_vuv = pd.read_csv(location+'/time_Ar_vuv.txt')
    
    # Concatenate all into a new DataFrame object
    time_Ar_vuv['vuv'] = 1
    time_Ar_vis['vuv'] = 0
    frames = [time_Ar_vuv,time_Ar_vis]
    timeAr = pd.concat(frames)
    timeAr = timeAr.sort_values('time')
    timeAr = timeAr.reset_index(drop=True)
    del time_Ar_vis, time_Ar_vuv
    
    # compute the number of Ar39 events
    nr_events_Ar = int(np.max(timeAr['event']))+1
        
    return timeAr,nr_events_Ar