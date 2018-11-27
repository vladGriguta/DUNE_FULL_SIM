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
import CandidateSearch



# Define constant parameters in the simulation
simulationTime = 10*1000000  # 10 seconds in microseconds
resolution = 0.05 # 50 nanoseconds


# Call functions to read from files
location = '../dataDUNE/10Seconds'
print('Reading data from files... Simulation time = '+str(simulationTime/1000000)+' seconds')
timeSN,nr_events_SN = readFiles.Read_SN(location)
timeAr,nr_events_Ar = readFiles.Read_Ar39(location)
print('Finished reading from files.')

# Compute the noise (mean number of photons from Ar39 per microseconds)
noise = float(nr_events_Ar)/simulationTime

print('The mean number of Ar39 events per microsecond is '+str(nr_events_Ar/simulationTime)
      +'. Compare with expected number of (about) 0.8 photons per microsecond.')
print('The mean number of Ar39 photons per microsecond is '+str(len(timeAr)/simulationTime)
      +'. Compare with expected number of (about) 1.8 photons per microsecond.')

### Testing the code

events, time_and_PD = genericFunctions.DivideDataByRes(timeSN,timeAr,simulationTime,resolution=resolution)

# For now, all events are generated in the centre of each time_frame (2.5 ms)
eventsSN = pd.DataFrame()
eventsSN['eventTime'] = np.linspace(0,int(simulationTime/2500)-1,int(simulationTime/2500))*(2500)+1250
eventsSN['event'] = np.linspace(0,int(simulationTime/2500)-1,int(simulationTime/2500))

import gc
del timeSN, timeAr
gc.collect()
"""
# Create another array to be used for imposing threshold on the distribution
# of photons among Photon Detectors
threshold2 = [1,2] # number of PDs,number of Photons
SNCandidates, fakeTrig, SNTrig = genericFunctions.Candidates(events=events, threshold=18,
                                            integTime=2.5, eventsSN=eventsSN, trigDur=2.5,
                                            resolution=resolution,noise=noise, threshold2=threshold2,
                                            time_and_PD=time_and_PD)
"""

thresholdVals = [18,22]
integTimeVals = [2,3]
threshold2Vals = [[1,5],[2,3]]
trigDur = 2.5


CS = CandidateSearch.Candidates(events, eventsSN, trigDur,resolution,noise,
                                time_and_PD,simulationTime)
GS = genericFunctions.GridSearch(thresholdVals, integTimeVals, threshold2Vals)

results = GS.ActualGridSearch()
print(results)


"""
# Save the pandas pannels
df1 = pd.Panel.to_frame(df_eff)
df1.to_csv('../dataDUNE/resultsFullSim/efficiencies_'+str(trigDur)+'.csv', sep=',')
df2 = pd.Panel.to_frame(df_fake)
df2.to_csv('../dataDUNE/resultsFullSim/fakeRate_'+str(trigDur)+'.csv', sep=',')

import plottingFunctions
locationPlots = '../dataDUNE/10SecondsResults'
plottingFunctions.plotGridSearch(panel=df_eff,noise=noise, trigDur=trigDur,
                                 location=locationPlots, isEff=True)
plottingFunctions.plotGridSearch(panel=df_fake,noise=noise, trigDur=trigDur,
                                 location=locationPlots, isEff=False)

"""