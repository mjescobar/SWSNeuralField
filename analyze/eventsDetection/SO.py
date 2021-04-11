#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:40:32 2019

@author: felipe
"""

import numpy as np
import scipy.signal as signal


def SODetection(data,neg_threshold=-40e-6,p2p_threshold=80e-6, SOp2p_threshold=0,fs=100):
    shapeData=np.shape(data)
    SOdetection=np.zeros_like(data)
    SOsegmentation=np.zeros_like(data)
    SOp2p=[]
    lenTime=shapeData[0]
    lowFreq=0.5
    highFreq=1.25
    b,a=signal.cheby1(4,1e-6,[lowFreq/(fs/2), highFreq/(fs/2)],'bandpass')
    filteredSignal=signal.filtfilt(b,a,data)
    zeros_down=[]
    SOs=[]
    zeros_prev=np.zeros((lenTime//100,),dtype=int)
    zeros_post=np.zeros((lenTime//100,),dtype=int)
    zero_post=0
    good_so=[]
    for n in range(lenTime):
        #Detect the zero-crosses from up to down in the filtered signal
        if n>0 and filteredSignal[n-1]>0 and filteredSignal[n]<0:
            zeros_down.append(n)
    len_zeros=len(zeros_down)
    zeros_down=np.array(zeros_down)
    nzero=0 #Counter of zeros
    
    if len(zeros_down)>0:
        n=zeros_down[nzero] #Start sample
    else:
        n=lenTime
        
    count=0 #Counter of SO
    while (n<lenTime and nzero<len_zeros-1):
        #Search for overpassing of the threshold
        if filteredSignal[n]<neg_threshold:
            if n>zeros_down[0]:
                #after the second zero-cross
                #select the previous zero-cross
                zero_prev=zeros_down[np.argwhere(zeros_down<n)[-1][0]]
                #Avoid overlaps
                if zero_prev<zero_post:
                    nzero+=1
                    n=zeros_down[nzero]
                    zero_prev=zeros_down[np.argwhere(zeros_down<n)[-1][0]]
            else:
                #Start at zero index time (this not guarantee a zero cross) 
                zero_prev=0
                
            if n<zeros_down[-1]:
                #Before the penultimate zero-cross
                #select the next zero-cross
                zero_post=zeros_down[np.argwhere(zeros_down>n)[0][0]]
            else:
                #finish at the last point (this not guarantee a zero cross) 
                zero_post=lenTime
                
            period=(zero_post-zero_prev)/fs #Calculate the period
            peak2peak=np.max(filteredSignal[zero_prev:zero_post])-np.min(filteredSignal[zero_prev:zero_post])
            # peak2peak=np.max(data[zero_prev:zero_post])-np.min(data[zero_prev:zero_post])
            if period>1.0/highFreq and period<1.0/lowFreq and peak2peak>p2p_threshold:
                count+=1
                zeros_prev[count]=zero_prev
                zeros_post[count]=zero_post
                SOsegmentation[zero_prev:zero_post]=count
                SOp2p.append(peak2peak)
            nzero+=1
            if nzero<len_zeros:
                n=zeros_down[nzero]
        else:
            if n+1<lenTime:
                if SOsegmentation[n+1]>0:
                    nzero+=1
                    n=zeros_down[nzero]
            n+=1
    #Refining
    nn=0
    meanSOp2p=0
    if count>0:
        SOp2p=np.array(SOp2p)
        meanSOp2p=np.mean(SOp2p)
        if SOp2p_threshold==0:
            #Set the threshold if it was not set
            SOp2p_threshold=1.25*meanSOp2p
        for n in range(1,len(SOp2p)+1):
            if SOp2p[n-1]<SOp2p_threshold:
                SOsegmentation[zeros_prev[n]:zeros_post[n]]=0
            else:
                good_so.append(n)
        for n in good_so:
            nn+=1    
            SOsegmentation[zeros_prev[n]:zeros_post[n]]=nn
            SOdetection[zeros_prev[n]:zeros_post[n]]=1
            SOs.append(np.argwhere(data[np.argwhere(SOsegmentation==nn)[:,0]]==np.min(data[np.argwhere(SOsegmentation==nn)[:,0]]))[0][0]+np.argwhere(SOsegmentation==nn)[0,0])
                    
    count=len(SOs)
    return SOdetection, SOsegmentation, SOs,count,meanSOp2p

def SOERP(data,SOs,prevTime=100,postTime=100):
    data_m=np.zeros((prevTime+postTime,))
    for i in range(len(SOs)):
        if (SOs[i]>prevTime and SOs[i]<6900):
            data_m+=data[SOs[i]-prevTime:SOs[i]+postTime]
    data_m=data_m/len(SOs)
    SOerp=data_m
    return SOerp        

# #Test 
# x=0.002*np.cos(1.9*np.pi*np.linspace(0,4,401))
# SOdetected=SODetection(x)