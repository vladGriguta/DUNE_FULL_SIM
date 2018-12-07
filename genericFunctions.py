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

    events = np.zeros((int(simulationTime/resolution)))    
    for i in range(0,len(time_and_PD)):
        events[int(time_and_PD[i][0]/resolution)] += 1
        #PDs[int(time_and_PD[i][0]/resolution)].append(time_and_PD[i][1])
        #events[int(time_and_PD[i][0]/resolution)].append(time_and_PD[i][1])
        
    return events, time_and_PD






class GridSearch:

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
    def __init__(self, thresholdVals, integTimeVals, threshold2Vals, events, eventsSN,
                 trigDur,resolution,noise,time_and_PD,simulationTime):
        self.thresholdVals = thresholdVals
        self.integTimeVals = integTimeVals
        self.threshold2Vals = threshold2Vals
        self.events=events
        self.eventsSN = eventsSN
        self.trigDur = trigDur
        self.resolution = resolution
        self.noise = noise
        self.time_and_PD = time_and_PD
        self.simulationTime = simulationTime
    
    def ActualGridSearch(self):

        #from pathos.multiprocessing import ProcessingPool as Pool
        
        thresholdVals = self.thresholdVals
        integTimeVals = self.integTimeVals
        threshold2Vals = self.threshold2Vals
        events = self.events
        eventsSN = self.eventsSN
        trigDur = self.trigDur
        resolution = self.resolution 
        noise = self.noise
        time_and_PD = self.time_and_PD
        simulationTime = self.simulationTime
        
        x = []
        y = []
        z = []
        
        for i in range(len(thresholdVals)):
            for j in range(len(integTimeVals)):
                for k in range(len(threshold2Vals)):
                    x.append(thresholdVals[i])
                    y.append(integTimeVals[j])
                    z.append(threshold2Vals[k])
        
        varyingData = []
        for i in range(len(x)):
            varyingData.append([x[i],y[i],z[i]])
        
        constantData = [events,eventsSN,trigDur,resolution,noise,time_and_PD,simulationTime]
        
        import CandidateSearch
        #results.append(Pool().map(CandidateSearch.Candidates.CandidateSearch(),[x,y,z]))
        import multiprocessing
        import itertools
        
        freeProc = 2
        n_proc=multiprocessing.cpu_count()-freeProc
        print('Code will run using '+str(n_proc)+' processors. This is equal to the total number of processors available minus '+str(freeProc))
        with multiprocessing.Pool(processes=n_proc) as pool:
            result_list=pool.starmap(CandidateSearch.Candidates, zip(varyingData, itertools.repeat(constantData)))
            pool.close()

        return result_list

    
    
    
