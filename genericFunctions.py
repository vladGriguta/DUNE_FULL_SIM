#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 20:41:47 2018

@author: vladgriguta
"""
import matplotlib.pyplot as plt
import numpy as np;
import pandas as pd

def DivideDataByRes(timeSN,timeAr,simulationTime,resolution=50):
    """
    This function transforms the array containing the arival time of all photons
    to an array containing the number of photons detected in each resolution
    time
    
    Input:
    timeSN:     DataFrame with all arrival times from SN events
    timeAr:     DataFrame with all arrival times from Ar39 events
    resolution: Minimum time in which the detector can distinguish separate photon
                hits
    simulationTime:
                Total duration of the simulation
                
    Output:
    events:     array of length int(simulationTime/resolution) showing how many photons
                were detected in a given time resolution and which PDs have they excited
    time_and_PD:2D array with arrival time and PD of all photons sorted by arrival time
    """
    time_and_PD = np.concatenate((np.array(timeSN[['time','pmt']]),np.array(timeAr[['time','pmt']])))
    time_and_PD = np.array(time_and_PD[time_and_PD[:,0].argsort()])
    # Need to solve the issue with the memory usage.
    # But maybe not this way. Maybe should think of a way to compute which PMTs
    # are being hit for each time frame inside Candidates
    
    """
    events = []
    for i in range(0,int(simulationTime/resolution)):
        events.append([])
        events[i].append(0)
    """
    events = np.zeros((int(simulationTime/resolution)))    
    for i in range(0,len(time_and_PD)):
        events[int(time_and_PD[i][0]/resolution)] += 1
        #PDs[int(time_and_PD[i][0]/resolution)].append(time_and_PD[i][1])
        #events[int(time_and_PD[i][0]/resolution)].append(time_and_PD[i][1])
        
    return events, time_and_PD






def GridSearch(events, eventsSN, thresholdVals, integTimeVals, threshold2Vals, trigDur,resolution,noise,time_and_PD,simulationTime):

    """
    Function that searches through a grid of threshold and integTime vals and
    returs a Dataframe object with the appropriate values of SN detection 
    efficiency and fake trigger rate
    
    Input variables:
    events:        A dataframe with timings of all photons recorded in the simulation
    eventsSN:      A dataframe with the decay time of all SN events
    thresholdVals: A numpy ARRAY containing the values of the threshold
    integTimeVals: A numpy ARRAY with integration time windows
    trigDur:       A numpy array with the freeze-out times of the triggering algorithm
    threshold2Vals:    A numpy array with the maximum number of photons to hit at
                   least one Photon Detector in integTime
    resolution:    The resolution assumed for the photodetector system
    noise:         Mean number of Ar39 photons per microsecond
    time_and_PD:   Numpy array with the timing and PD hit by all photons recorded
    simulationTime:Obvious
    
    Output variables:
    df_eff, df_fake:
        Pandas Panels (3D Dataframes) storing the SN detection efficiency
        and the rate of fake triggers per second for each combination of the 
        variables in the grid search
    """
    
    def sendStaticData():
        return events, eventsSN, trigDur,resolution,noise,time_and_PD,simulationTime
    

    index1 = list(map(str, thresholdVals))
    index2 = list(map(str, integTimeVals))
    index3 = list(map(str, threshold2Vals))
    
    df_eff = pd.Panel(items = index1,major_axis = index2, minor_axis = index3)
    df_fake = pd.Panel(items = index1,major_axis = index2, minor_axis = index3)
    """
    progress = 0    
    for i in range(0,len(thresholdVals)):
        for j in range(0,len(integTimeVals)):
            for k in range(0,len(threshold2Vals)):
                
                # Print progress if applicable
                if(((i+1)*(j+1)*(k+1)) % int((len(thresholdVals)+1)*(len(integTimeVals)+1)
                    *(len(threshold2Vals)+1)/10) == 0):
                    progress += 10
                    print(str(progress) + ' % Completed!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                
                threshold_local = thresholdVals[i]
                integTime_local = integTimeVals[j]
                threshold2_local = threshold2Vals[k]
                
                [SNCandidates, fakeTrig, _] = Candidates(events, eventsSN, threshold_local, 
                    integTime_local, threshold2_local, trigDur,resolution,noise,time_and_PD)
                
                df_eff[str(threshold_local)].iloc[j][k] = (100 * np.sum(SNCandidates>0)/len(SNCandidates))
                df_fake[str(threshold_local)].iloc[j][k] = (len(fakeTrig)* 1000000 / simulationTime)
                print('This run was finished')
                
    """
    
    from pathos.multiprocessing import ProcessingPool as Pool

    x = thresholdVals
    y = integTimeVals
    z = threshold2Vals
    
    y_,x_,z_ = np.meshgrid(y,x,z)
    x = x_.flatten()
    y = y_.flatten()
    z = z_.flatten()
    
    results = []
    
    import CandidateSearch
    results.append(Pool().map(CandidateSearch.Candidates,x,y,z))
    

    
    return df_eff, df_fake, results
