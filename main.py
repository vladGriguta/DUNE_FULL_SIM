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
simulationTime = 2.5*1000000  # 2.5 seconds in microseconds
resolution = 0.05 # 50 nanoseconds


# Call functions to read from files
location = 'oldQuantumEff'
timeSN,nr_events_SN = readFiles.Read_SN(location)
timeAr,nr_events_Ar = readFiles.Read_Ar39(location)

# Compute the noise (mean number of photons from Ar39 per microseconds)
noise = float(nr_events_Ar)/simulationTime

print('The number of Ar39 events in '+str(simulationTime)+' seconds is '
      +str(nr_events_Ar)+ '. Compare with expected number of (about) 1.8 photons per us.')

### Testing the code

events, time_and_PD = genericFunctions.DivideDataByRes(timeSN,timeAr,simulationTime,resolution=resolution)

# For now, all events are generated in the centre of each time_frame (2.5 ms)
eventsSN = pd.DataFrame()
eventsSN['eventTime'] = np.linspace(0,int(simulationTime/2500)-1,int(simulationTime/2500))*(2500)+1250
eventsSN['event'] = np.linspace(0,int(simulationTime/2500)-1,int(simulationTime/2500))

# Create another array to be used for imposing threshold on the distribution
# of photons among Photon Detectors

threshold2 = [1,2] # number of PDs,number of Photons

SNCandidates, fakeTrig, SNTrig = genericFunctions.Candidates(events=events, threshold=13,
                                            integTime=2.5, eventsSN=eventsSN, trigDur=2.5,
                                            resolution=resolution,noise=noise, threshold2=threshold2,
                                            time_and_PD=time_and_PD)

SNEff = np.count_nonzero(SNCandidates)/len(SNCandidates)




