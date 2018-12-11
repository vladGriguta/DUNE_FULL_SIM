#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 20:41:11 2018

@author: vladgriguta
"""


# This is the high level coding structure where all functions
# are called

import numpy as np
import pandas as pd
import genericFunctions
import readFiles



# Define constant parameters in the simulation
simulationTime = 100*1000000  # 10 seconds in microseconds
resolution = 0.05 # 50 nanoseconds


# Call functions to read from files
location = '../dataDUNE/10Seconds'
print('Reading data from files... Simulation time = '+str(simulationTime/1000000)+' seconds')
"""
timeSN,nr_events_SN = readFiles.Read_SN(location)
timeAr,nr_events_Ar = readFiles.Read_Ar39(location)
"""
timeAr,nr_events_Ar = readFiles.Read_Ar39_Oana()
timeSN,nr_events_SN = readFiles.Read_SN_Oana()

print('Finished reading from files.')

# Compute the noise (mean number of photons from Ar39 per microseconds)
noise = float(len(timeAr))/simulationTime

print('The mean number of Ar39 events per microsecond is '+str(nr_events_Ar/simulationTime)
      +'. Compare with expected number of (about) 0.8 photons per microsecond.')
print('The mean number of Ar39 photons per microsecond is '+str(len(timeAr)/simulationTime)
      +'. Compare with expected number of (about) 1.8 photons per microsecond.')

### Testing the code

events, time_and_PD = genericFunctions.DivideDataByRes(timeSN,timeAr,simulationTime,resolution=resolution)
print('The data from Ar and SN has been saved to the new array')

# For now, all events are generated in the centre of each time_frame (2.5 ms)
eventsSN = pd.DataFrame()
eventsSN['eventTime'] = np.linspace(0,int(simulationTime/2500)-1,int(simulationTime/2500))*(2500)+1250
eventsSN['event'] = np.linspace(0,int(simulationTime/2500)-1,int(simulationTime/2500))

# Delete unnecesary data 
import gc
del timeSN, timeAr
gc.collect()
"""
SNCandidates, fakeTrig, SNTrig = genericFunctions.Candidates(events=events, threshold=18,
                                            integTime=2.5, eventsSN=eventsSN, trigDur=2.5,
                                            resolution=resolution,noise=noise, threshold2=threshold2,
                                            time_and_PD=time_and_PD)
"""

thresholdVals =  [0.468 *200,0.468 *240,0.468 *280]
integTimeVals = [0.3,0.4,0.5,0.6,0.7]
threshold2Vals = [[0,0],[3,6],[4,6],[5,6],[6,6],[7,6],[8,6]] # [number_PDs,PhotonRateInEach]
trigDur = 0

# Create an instance of the Grid search class which saves all relevant attributes
# for the analysis

GS = genericFunctions.GridSearch(thresholdVals, integTimeVals, threshold2Vals,events,
                                 eventsSN, trigDur,resolution,noise,time_and_PD,simulationTime)

results = GS.ActualGridSearch()

print('The results of the full simulation are:\n'+str(results))

locationResults = '../dataDUNE/'
with open(locationResults+'resultsParallelImp_11.12.18_.txt', 'w') as f:
    for item in results:
        f.write("%s\n" % str(item))
