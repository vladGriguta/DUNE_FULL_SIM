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


def Candidates(events, eventsSN, threshold, integTime, threshold2, trigDur,resolution,noise,time_and_PD):
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
            """
            else:
                print('The event did not succeed due to second threshold')
            """
                    
        photonsPassed += int(events[i-1-integBins])
    
    return SNCandidates, fakeTrig, SNTrig


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

    index1 = list(map(str, thresholdVals))
    index2 = list(map(str, integTimeVals))
    index3 = list(map(str, threshold2Vals))
    
    df_eff = pd.Panel(items = index1,major_axis = index2, minor_axis = index3)
    df_fake = pd.Panel(items = index1,major_axis = index2, minor_axis = index3)

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
                
    return df_eff, df_fake





# Visualizing the results
def Plot_Trigger_Distrib(eventsSN,trigEf,fakeRate,threshold,SN_event_nr_bins,mean_events,SN_event_time):
    # Plot the SN events that did not pass trigger
    X1 = []
    Y1 = []
    c1 = []
    X2 = []
    Y2 = []
    c2 = []
    
    for i in range(0,len(eventsSN)):
        if(eventsSN['TriggerResponse'][i]):
            X1.append(eventsSN['energy'][i])
            Y1.append(eventsSN['distanceToAnode'][i])
            c1.append('blue')
        else:
            X2.append(eventsSN['energy'][i])
            Y2.append(eventsSN['distanceToAnode'][i])
            c2.append('red')
    plt.scatter(X1,Y1,c = c1,label = 'Detected',alpha=0.7)
    plt.scatter(X2,Y2,c = c2,label = 'Undetected',alpha=0.7)
    plt.legend()
    plt.title('Scatter plot of SN events')
    plt.xlabel('Energy / MeV')
    plt.ylabel('Distance to anode / cm')
    textstr = '\n'.join((
        r'$\mathrm{SNEfficiency}=%.2f $' % (trigEf, )+'%',
        r'$\mathrm{FakeEventsRate}=%d s^{-1} $' % (fakeRate, )))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    plt.text(25, 375, textstr, fontsize=12,
            verticalalignment='top', bbox=props)
    thr = int(threshold / (SN_event_nr_bins*mean_events))
    plt.savefig('week5/SNtime'+str(SN_event_time)+'.thr'+str(thr) +'.jpg', format='jpg')