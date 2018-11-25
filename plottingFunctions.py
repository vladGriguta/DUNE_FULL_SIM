#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 24 19:31:57 2018

@author: vladgriguta
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

def plotGridSearch(panel,noise,trigDur,location='../dataDUNE/10SecondsResults',isEff=True):
    """
    This is the function designed for plotting the grid search values. It takes
    as input the panel (3D DataFrame) containing SN efficiencies (or fake triggering
    rates) for each value along the three axis (items=threshold1;
    major_axis=integration time; minor_axis=threshold2) and it delivers 3D plots 
    of this data.
    
    
    """
    # In next part the panel characteristics are converted to 3D arrays that 
    # can be plotted using the matplotlib library    
    ###########################################################################
    
    temp = panel.items
    th1 = np.zeros(np.shape(temp))
    for i in range(len(temp)):
        th1[i] = int(temp[i])
    temp = panel.major_axis
    inTime = np.zeros(np.shape(temp))
    for i in range(len(temp)):
        inTime[i] = int(temp[i])
    th2 = panel.minor_axis
    
    possible_colors = ['red','orange','yellow','green','blue','purple']
    colors = possible_colors[0:len(th2)]
    #colors = cm.rainbow(np.linspace(start=0, stop=1, num=len(th2)))
    
    th2_to_color = {}
    for key, val in zip(th2,colors):
        th2_to_color[key] = val
    
    vals = panel.values
    ##########################################################################
    
    
    th1_, inTime_ = np.meshgrid(th1, inTime, indexing='ij')
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    
    # add each point individually
    groups = th2
    for i in range(len(groups)):
        if(isEff == True):
            ax.scatter(th1_,inTime_,vals[:][:][i],c=th2_to_color[groups[i]],label=groups[i],linewidth=0.5,alpha=0.8)
        else:
            # For fake trigger rate represent in logarithmic scale
            ax.scatter(th1_,inTime_,np.log10(vals[:][:][i]),c=th2_to_color[groups[i]],label=groups[i],linewidth=0.5,alpha=0.8)
        
    plt.legend(loc=2) 
    ax.set_xlabel('Threshold1 / \u03C3')
    ax.set_ylabel('Integration Time / seconds')

    
    
    if(isEff == True):
        ax.set_zlabel('SN Detection Efficiency')
        fig.savefig(location+'efficiency_trigDur='+str(trigDur)+'.png',format='png')
    else:
        ax.set_zlabel('Log(Rate of False Triggers)')
        fig.savefig(location+'fakeRate_trigDur='+str(trigDur)+'.png',format='png')
    
    
    
    
    
    
    
    
    
    