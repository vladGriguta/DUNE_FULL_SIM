#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 24 19:31:57 2018

@author: vladgriguta
"""
import pandas as pd
import numpy as np


def plotGridSearch(df_eff,df_fake):
    
    eff = pd.Panel.to_xarray(df_eff)
    #fake = pd.Panel.to_xarray(df_fake)
    
    for i in range(0,len(eff.items)):
        for j in range(0,len(eff.major_axis)):
            for k in range(0,len(eff.minor_axis)):
                eff[i][j][k] = i+j+k
    
    th1 = eff.items
    th1_ = np.zeros(len(th1))
    for i in range(len(th1)):
        th1_[i] = int(th1[i])
        
    inTime = eff.major_axis
    inTime_ = np.zeros(len(inTime))
    for i in range(len(inTime)):
        inTime_[i] = int(inTime[i])
        
    th2 = eff.minor_axis
    th2_ = []
    for i in range(len(th2)):
        th2_.append(str(th2[i]))
    
    import matplotlib.cm as cm
    colors = cm.rainbow(np.linspace(start=0, stop=1, num=len(th2)))
    
    keys = th2_
    vals = colors
    myDict = {}
    for key, val in zip(keys,vals):
        myDict[key] = val
    
    vals = eff.values
    vals_ = np.zeros(np.shape(vals))
    for i in range(np.shape(vals)[0]):
        for j in range(np.shape(vals)[1]):
            for k in range(np.shape(vals)[2]):
                vals_[i][j][k] = int(vals[i][j][k])
    
    
    
    th1, inTime, th2 = np.meshgrid(th1_, inTime_,th2_, indexing='ij')
    
    colors = np.zeros(np.concatenate((np.array(np.shape(th2)),np.array([4]))))
    for i in range(np.shape(th2)[0]):
        for j in range(np.shape(th2)[1]):
            for k in range(np.shape(th2)[2]):
                colors[i][j][k] = myDict[th2[i][j][k]]
    
    
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    
    ax.scatter(th1,inTime,vals,c=colors.reshape(8,4),linewidth=0.5,alpha=0.8)