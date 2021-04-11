#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 21:15:27 2021

@author: felipe
"""

import sys, os
sys.path.append(os.path.abspath('../spectrum'))
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import Data
import scipy.signal as signal

prop_cycle = plt.rcParams['axes.prop_cycle']
colorsProp = prop_cycle.by_key()['color']

filenames=['/home/felipe/Dropbox/UTFSM/NeuralField/timeseries/baselines/Sleep-caso-1-Phi.txt',
           '/home/felipe/Dropbox/UTFSM/NeuralField/timeseries/baselines/Sleep-caso-2-Phi.txt',
           '/home/felipe/Dropbox/UTFSM/NeuralField/timeseries/baselines/Sleep-caso-3-Phi.txt',
           '/home/felipe/Dropbox/UTFSM/NeuralField/timeseries/baselines/Sleep-caso-4-Phi.txt',
           '/home/felipe/Dropbox/UTFSM/NeuralField/timeseries/baselines/Sleep-caso-5-Phi.txt'
           ]


for num,filename in enumerate(filenames):
    data, datazScore,std=Data.loadStoredData(filename)
    data=data.to_numpy(dtype=float)
    data=np.mean(data,axis=1)-np.mean(data)
    lowFrequencySP=9
    highFrequencySP=16
    lowFrequencySO=0.5
    highFrequencySO=1.25
    fs=100
    bSP,aSP=signal.cheby1(4,1e-6,[lowFrequencySP/(fs/2), highFrequencySP/(fs/2)],'bandpass')
    SPData=signal.filtfilt(bSP,aSP,data,axis=0)
    bSO,aSO=signal.cheby1(4,1e-6,[lowFrequencySO/(fs/2), highFrequencySO/(fs/2)],'bandpass')
    SOData=signal.filtfilt(bSO,aSO,data,axis=0)
    hm=signal.hamming(100)
    SPrms=np.convolve(np.sqrt(SPData**2),hm/np.sum(hm**2))
    SOrms=np.convolve(np.sqrt(SOData**2),hm/np.sum(hm**2))
    # if num==3:
    #     print(1.25*np.mean(SOrms))
    #     print(1.25*np.mean(SPrms))
    # th_SO=0.2215394096435029
    # th_SP=0.15249695533902907
    
    
    #%%
    fig1=plt.figure(figsize=(2,1))
    ax=fig1.add_subplot(frame_on=False)
    ax.tick_params(axis='both',which='both',length=0)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if num<2:
        offset=1e-2
        th=5e-4
        ax.vlines(0.5,-1e-2,1e-2,colors='k')
        ax.hlines([-0.96e-2,0,0.96e-2],0,0.5,colors='k')
        SOrms1=np.ma.masked_less_equal(SOrms, th)
        SPrms1=np.ma.masked_less_equal(SPrms, th)
    else:
        offset=2e-2
        th=2e-3
        ax.vlines(0.5,-3e-2,3e-2,colors='k')
        ax.hlines([-2.88e-2,-2e-2,-1e-2,0,1e-2,2e-2,2.89e-2],0,0.5,colors='k')
        SOrms1=np.ma.masked_less_equal(SOrms, th)
        SPrms1=np.ma.masked_less_equal(SPrms, th)
        
    ax.plot(np.linspace(0.5,20,1869),data[130::],'k',linewidth=0.35)
    plt.plot(np.linspace(0.5,20,1968),SOrms[130::]+offset,linewidth=2,color=colorsProp[1])
    plt.plot(np.linspace(0.5,20,1968),SPrms[130::]+offset,linestyle=(2,(1.5,0.25)),linewidth=2,color=colorsProp[4])
    ax.hlines(0,0,20,colors='k',linewidths=1)
    
    fig1.savefig(filename+'-excerptPlot.pdf',dpi=300)
    
    
    
    
    
    
    
    
    
    
    