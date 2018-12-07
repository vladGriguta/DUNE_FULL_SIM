#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 21:04:34 2018

@author: vladgriguta
"""

import matplotlib.pyplot as plt
import numpy as np;
import pandas as pd

def checkRecentTrig(time_of_event,SNTrig,fakeTrig,trigDur):
    """
    Function used inside Candidates() to check whether the trigger was recently
    activated.
    Return True if trigger was recently activated and False otherwise
    """
    # Use try-except to avoid error when arrays are empty
    try:
        if(time_of_event-SNTrig[len(SNTrig)-1]<trigDur):
            return True
            #print(time_of_event-SNTrig[len(SNTrig)-1])
    except:
        pass
    try:
        if(time_of_event-fakeTrig[len(fakeTrig)-1]<trigDur):
            return True
            #print(time_of_event-fakeTrig[len(fakeTrig)-1])
    except:
        pass
    
    return False

def secondCondition(threshold2,time_and_PD,photons_in_region,photonsPassed):
    """
    Function that checks if at least threshold2[0] PMTs were hit by threshold2[1]
    photons.
    
    """
    PD_distrib = np.zeros(120)
    for i in range(0,photons_in_region):
        # Start at the position photonsPassed
        PD_distrib[int(time_and_PD[photonsPassed+i][1])] += 1
    
    PDs_meeting_requirement = 0
    for i in range(0,120):
        if(PD_distrib[i]>threshold2[1]):
            PDs_meeting_requirement += 1
            if(PDs_meeting_requirement >= threshold2[0]):
                return True
    return False



def Candidates(varyingData,constantData):
    """
    This function searches through the 
    
    The input variables are: 
    events:     A dataframe with timings of all photons recorded in the simulation
    eventsSN:   A dataframe with the decay time of all SN events
    noise:      The mean number of Ar39 events per microsecond
    threshold:  The threshold of photons that activates the trigger
    integTime:  The integration time window
    trigDur:    The freeze-out time of the triggering algorithm
    threshold2: The maximum number of photons to hit at least n Photon Detectors in integTime
                (n_PDs, n_photons)
    resolution: The resolution assumed for the photodetector system
    noise:      Mean number of Ar39 photons per microsecond
    time_and_PD:Numpy array with the timing and PD hit by all photons recorded
    
    The output is:
        
    SNCandidates: 1D Array containing the number of times each SN event is flagged
                  size = number of SN events
    fakeTrig:     1D Array containing the time when the triggering algorithm flagged
                  a SN event incorectly
    SNTrig:       1D Array containing the time when the triggering algorithm flagged
                  a SN event correctly
    """
    [threshold,integTime,threshold2] = varyingData
    
    [events,eventsSN,trigDur,resolution,noise,time_and_PD,simulationTime] = constantData
    

    
    # Declare output variables
    SNCandidates = np.zeros(len(eventsSN))
    SNTrig = []
    fakeTrig = []
    
    # Update threshold to account for different integTime
    threshold_in_photons = threshold * noise * integTime 
    
    # Compute the number of bins corresponding to the integTime
    integBins = int(integTime/resolution)
    
    # Start by computing how many events there are in the first integTime
    photons_in_region = int(np.sum(events[0:integBins]))
    # And which PDs are being hit
    photonsPassed = 0
    
    # Go through the array and get those points where the total number of events
    # counted is above the threshold_in_photons
    progress = 0
    for i in range(integBins+1,len(events)):
        # Start by printing the progress
        if(i % int(len(events)/10) == 0):
            progress += 10
            print(str(progress) + ' % Completed for threshold='+str(threshold)+
                  ', integrationTime = '+str(integTime)+
                  ', threshold2 = '+str(threshold2)+
                  ', trigDuration = '+str(trigDur))
        
        # update number of events by substracting the element furthest away from
        # current position and adding the element in current position
        photons_in_region = int(photons_in_region - events[i-1-integBins] + events[i])
        
        
        if(photons_in_region>threshold_in_photons):
            
            if(secondCondition(threshold2,time_and_PD,photons_in_region,photonsPassed) == True):
                # The time when the event is flagged as SN is defined as
                time_of_event = (i-float(integBins/2))*resolution
                
                # Check if trigger has not activated recently
                recentTrig = checkRecentTrig(time_of_event,SNTrig,trigDur,trigDur)
                
                if(not recentTrig):            
                    # Check if flag is within vicinity of an actual SN event
                    if(np.min(abs(time_of_event-eventsSN['eventTime'])) < integTime):
                        SNCandidates[np.argmin(abs(time_of_event-eventsSN['eventTime']))] +=1
                        SNTrig.append(time_of_event)
                    else:
                        fakeTrig.append(time_of_event)
                    
        photonsPassed += int(events[i-1-integBins])
    
    eff = (100 * np.sum(SNCandidates>0)/len(SNCandidates))
    fakeRate = (len(fakeTrig)* 1000000 / simulationTime)
    
    print('\n\n\n The results for the following inputs: \nthreshold='+str(threshold)+'\nintegrationTime = '+str(integTime)+
      '\n threshold2 = '+str(threshold2)+'\n trigDuration = '+str(trigDur)+ '\nare:\n'+
      'Efficiency of SN detection: '+ str(eff)+ '\nNumber of fake triggers: '+str(fakeRate)+'\n')
    
    locationResults = '../dataDUNE/10SecondsResults/'+'ParallelImp_trig'+str(trigDur)+'_v1_'+'.txt'
    f = open(locationResults, 'a+')
    f.write(str(threshold)+','+str(integTime)+','+str(threshold2)+','+str(eff)+','+str(fakeRate)+'\n')
    f.close()
    
    return eff, fakeRate