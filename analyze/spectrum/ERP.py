#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 20:37:59 2019

@author: felipe
"""
import sys, os
sys.path.append(os.path.abspath('../spectrum'))
import numpy as np
import scipy.signal as signal
import scipy.stats as stats
from scipy import ndimage
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import Data
import Epochs
import Wavelets
from matplotlib import rc
rc('text', usetex=True)



def getEpochs(data,markers,fs,shamData=0,startTime=-1,endTime=2,baseline=[-1,-0.5]):
    #inialization of variables
    nonZeroIndexes=np.nonzero(markers)[0]
    Epochs1=Epochs.Epochs(1,fs,startTime,endTime)
    Epochs2=Epochs.Epochs(2,fs,startTime,endTime)
    Epochs3=Epochs.Epochs(3,fs,startTime,endTime)
    Epochs4=Epochs.Epochs(4,fs,startTime,endTime)
    if shamData is not 0:
        EpochsS1=Epochs.Epochs(1,fs,startTime,endTime)
        EpochsS2=Epochs.Epochs(2,fs,startTime,endTime)
        EpochsS3=Epochs.Epochs(3,fs,startTime,endTime)
        EpochsS4=Epochs.Epochs(4,fs,startTime,endTime)
    #Time offsets
    pre =int( np.ceil(startTime*fs))
    post = int(np.ceil(endTime*fs)+1)
    #Procedure
    if len(nonZeroIndexes)>0:
		#if there are markers
        #correct the marker with valid offset time
        if pre+nonZeroIndexes[0]<0:
            initialMarker=1
        else:
            initialMarker=0
        if post+nonZeroIndexes[-1]>np.shape(data)[0]:
            finalMarker=len(nonZeroIndexes)-1
        else:
            finalMarker=len(nonZeroIndexes)
        
        for mark in nonZeroIndexes[initialMarker:finalMarker]:
            #For each marker
            if markers[mark]==1:
                Epochs1.addData(signal.detrend(data[mark+pre:mark+post]))
                if shamData is not 0:
                    EpochsS1.addData(signal.detrend(shamData[mark+pre:mark+post]))
            elif markers[mark]==2:
                Epochs2.addData(signal.detrend(data[mark+pre:mark+post]))
                if shamData is not 0:
                    EpochsS2.addData(signal.detrend(shamData[mark+pre:mark+post]))
            elif markers[mark]==3:
                Epochs3.addData(signal.detrend(data[mark+pre:mark+post]))
                if shamData is not 0:
                    EpochsS3.addData(signal.detrend(shamData[mark+pre:mark+post]))
            elif markers[mark]==4:
                Epochs4.addData(signal.detrend(data[mark+pre:mark+post]))
                if shamData is not 0:
                    EpochsS4.addData(signal.detrend(shamData[mark+pre:mark+post]))
        #Correct the baseline
        Epochs1.baselineCorrect(baseline)
        Epochs2.baselineCorrect(baseline)
        Epochs3.baselineCorrect(baseline)
        Epochs4.baselineCorrect(baseline)
        if shamData is not 0:
            EpochsS1.baselineCorrect(baseline)
            EpochsS2.baselineCorrect(baseline)
            EpochsS3.baselineCorrect(baseline)
            EpochsS4.baselineCorrect(baseline)
            return  Epochs1, Epochs2, Epochs3,Epochs4,EpochsS1, EpochsS2, EpochsS3,EpochsS4
        else:
            return Epochs1, Epochs2, Epochs3,Epochs4

def ERP(Epoch,fs=100,NxNy=256,lowFreq=0.1,highFreq=16):
    #ERP from all channels
    b,a=signal.cheby1(4,0.01,[lowFreq/(fs/2), highFreq/(fs/2)],'bandpass')
    
    #ERP in each channel
    splitPoints=np.arange(NxNy,np.shape(Epoch.data)[1],NxNy)  
    meansplitEpochs=np.mean(np.split(Epoch.data,splitPoints,axis=1),axis=0)
    erpChannels=signal.filtfilt(b,a,meansplitEpochs,axis=0)
    #Mean
    erp=np.mean(meansplitEpochs,axis=1)
    erp=signal.filtfilt(b,a,erp)
    return erp,erpChannels

def ERPSingle(Epoch,fs=100,lenTime=801,lowFreq=0.1,highFreq=16):
    #ERP from all channels
    splitPoints=np.arange(lenTime,np.shape(Epoch.data)[0],lenTime)  
    b,a=signal.cheby1(4,0.01,[lowFreq/(fs/2), highFreq/(fs/2)],'bandpass')
    erp=np.mean(np.split(Epoch.data,splitPoints),axis=0)
    stderp=np.std(np.split(Epoch.data,splitPoints),axis=0)
    erp=signal.filtfilt(b,a,erp)
    return erp

def tTestERP(Epoch1,Epoch2,NxNy=256):
	#Split data
    splitPoints=np.arange(NxNy,np.shape(Epoch1.data)[1],NxNy)  
    erpsEpoch1=np.mean(np.split(Epoch1.data,splitPoints,axis=1),axis=2)
    erpsEpoch2=np.mean(np.split(Epoch2.data,splitPoints,axis=1),axis=2)
    
    #Calcualte statisitical significant
    tvalue,pvalue=stats.ttest_ind(erpsEpoch1,erpsEpoch2,axis=0)
    return tvalue,pvalue

def tTestERPSingle(Epoch1,Epoch2,lenTime=801):
	#Split data
    splitPoints=np.arange(lenTime,np.shape(Epoch1.data)[0],lenTime)  
    erpsEpoch1=np.split(Epoch1.data,splitPoints,axis=0)
    erpsEpoch2=np.split(Epoch2.data,splitPoints,axis=0)
    #Calcualte statisitical significant
    tvalue,pvalue=stats.ttest_ind(erpsEpoch1,erpsEpoch2,axis=0)
    return tvalue,pvalue
	
def ERSP(Epoch,fs=100,NxNy=256):
    #Split each 256 nodes and average over them
    splitPoints=np.arange(NxNy,np.shape(Epoch.data)[1],NxNy)
    meansplitEpochs=np.mean(np.split(Epoch.data,splitPoints,axis=1),axis=2)
    #spectrogram window properties
    window=signal.get_window('hamming', int(2*fs))
    noverlap=int(1.9*fs)

    for n in range(len(splitPoints)+1):
        #Get the spectrogram of each trial
        f,t,Sxx_n=signal.spectrogram(meansplitEpochs[n,:],fs=fs, window=window, noverlap=noverlap, mode='complex')
        if n==0:    
            Sxx=np.zeros((len(splitPoints)+1,len(f),len(t)),dtype=complex)
        #Indexes: trial:freqs:time
        Sxx[n,:,:]=Sxx_n[0:len(f),:]
    #Calculate the average over trials
    ersp=np.mean(np.abs(Sxx),axis=0)
    return f,t,Sxx,ersp


def waveletERSP(Epoch,fs=100,NxNy=256):
    splitPoints=np.arange(NxNy,np.shape(Epoch.data)[1],NxNy)
    meansplitEpochs=np.mean(np.split(Epoch.data,splitPoints,axis=1),axis=2)
    # mother = pycwt.Morlet(6)
    # s0 = 2/fs  # Starting scale
    # octaves_suboctaves = 1 / 8 # Twelve sub-octaves per octaves
    # N = 7 / octaves_suboctaves  # Eigth powers of two per total sub-octave
    freqs = np.logspace(-1,1.478,300)
    dt=1/fs
		
    for n in range(len(splitPoints)+1):
        Sxx_n=Wavelets.Morlet(meansplitEpochs[n,:],freqs=freqs,dt=dt,omega0=15)
        scales=Sxx_n.getscales()	
        coefs=Sxx_n.getdata()
        coefs/=scales[:,None]
        if n==0:    
            Sxx=np.zeros((len(splitPoints)+1,len(freqs),np.shape(meansplitEpochs)[1]),dtype=complex)
        #trial:freqs:time
        Sxx[n,:,:]=coefs[0:len(freqs),:]
    ersp=np.mean(np.abs(Sxx),axis=0)
    t=np.arange(0,np.shape(meansplitEpochs)[1]/fs,1/fs)
    return np.flip(freqs),t,np.flip(Sxx),np.flipud(ersp)

def waveletERSPSingle(Epoch,fs=100,lenTime=801,correctF=True):
    splitPoints=np.arange(lenTime,np.shape(Epoch.data)[0],lenTime)
    meansplitEpochs=np.split(Epoch.data,splitPoints,axis=0)
    # mother = pycwt.Morlet(6)
    # s0 = 2/fs  # Starting scale
    # octaves_suboctaves = 1 / 8 # Twelve sub-octaves per octaves
    # N = 7 / octaves_suboctaves  # Eigth powers of two per total sub-octave
    freqs = np.logspace(-1,1.478,300)
    dt=1/fs
		
    for n in range(len(splitPoints)+1):
        Sxx_n=Wavelets.Morlet(meansplitEpochs[n],freqs=freqs,dt=dt,omega0=15)
        scales=Sxx_n.getscales()	
        coefs=Sxx_n.getdata()
        if correctF==True:
            coefs/=freqs[:,None]
        if n==0:    
            Sxx=np.zeros((len(splitPoints)+1,len(freqs),np.shape(meansplitEpochs)[1]),dtype=complex)
        #trial:freqs:time
        Sxx[n,:,:]=coefs[0:len(freqs),:]
    ersp=np.mean(np.abs(Sxx),axis=0)
    t=np.arange(0,np.shape(meansplitEpochs)[1]/fs,1/fs)
    return np.flip(freqs),t,np.flip(Sxx),np.flipud(ersp)

def waveletITC(Epoch,fs=100,NxNy=256,lowFreq=0.5,highFreq=16):
    
    f,t,fourierEpochs,ersp=waveletERSP(Epoch,fs=fs,NxNy=NxNy)
    frequencies=np.arange(np.argwhere(f>=lowFreq)[-1][0],np.argwhere(f>=highFreq)[-1][0])
    itcp=np.zeros((len(frequencies),np.shape(fourierEpochs)[2]))
    itcl=np.zeros((len(frequencies),np.shape(fourierEpochs)[2]))
    itcp_f=fourierEpochs[:,frequencies,:]/np.abs(fourierEpochs[:,frequencies,:]) 
    itcp[:,:]=np.abs(np.sum(np.angle(itcp_f),axis=0))/(np.shape(fourierEpochs)[0]*np.pi)
    itcl_f=np.sum(fourierEpochs[:,frequencies,:],axis=0)/np.sqrt(np.shape(fourierEpochs)[0]*np.sum(np.abs(fourierEpochs[:,frequencies,:])**2,axis=0))
    itcl[:,:]=np.abs(itcl_f)
    return f[np.argwhere(f>=lowFreq)[-1][0]:np.argwhere(f>=highFreq)[-1][0]],t,itcp, itcl

def waveletITCSingle(Epoch,fs=100,lenTime=801,lowFreq=0.5,highFreq=16):
    
    f,t,fourierEpochs,ersp=waveletERSPSingle(Epoch,fs=fs,lenTime=lenTime)
    frequencies=np.arange(np.argwhere(f>=highFreq)[-1][0],np.argwhere(f>=lowFreq)[-1][0])
    itcp=np.zeros((len(frequencies),np.shape(fourierEpochs)[2]))
    itcl=np.zeros((len(frequencies),np.shape(fourierEpochs)[2]))
    itcp_f=fourierEpochs[:,frequencies,:]/np.abs(fourierEpochs[:,frequencies,:]) 
    itcp[:,:]=np.abs(np.sum(np.angle(itcp_f),axis=0))/(np.shape(fourierEpochs)[0]*np.pi)
    itcl_f=np.sum(fourierEpochs[:,frequencies,:],axis=0)/np.sqrt(np.shape(fourierEpochs)[0]*np.sum(np.abs(fourierEpochs[:,frequencies,:])**2,axis=0))
    itcl[:,:]=np.abs(itcl_f)
    return f[np.argwhere(f>=highFreq)[-1][0]:np.argwhere(f>=lowFreq)[-1][0]],t,itcp, itcl

def ITC(Epoch,fs=100,NxNy=256,lowFreq=0.5,highFreq=16):
    
    f,t,fourierEpochs,ersp=ERSP(Epoch,fs=fs,NxNy=NxNy)
    frequencies=np.arange(np.argwhere(f>=lowFreq)[0][0],np.argwhere(f>=highFreq)[0][0])
    itcp=np.zeros((len(frequencies),np.shape(fourierEpochs)[2]))
    itcl=np.zeros((len(frequencies),np.shape(fourierEpochs)[2]))
    itcp_f=fourierEpochs[:,frequencies,:]/np.abs(fourierEpochs[:,frequencies,:]) 
    itcp[:,:]=np.abs(np.sum(np.angle(itcp_f),axis=0))/np.shape(fourierEpochs)[0]
    itcl_f=np.sum(fourierEpochs[:,frequencies,:],axis=0)/np.sqrt(np.shape(fourierEpochs)[0]*np.sum(np.abs(fourierEpochs[:,frequencies,:])**2,axis=0))
    itcl[:,:]=np.abs(itcl_f)
    return f[np.argwhere(f>=lowFreq)[0][0]:np.argwhere(f>=highFreq)[0][0]],t,itcp, itcl
    
    #return itc



def plotERPS(Epoch1,Epoch2,fs=100,label1='STIM', label2='SHAM', startTime=-0.5, endTime=2,savefile='ERP.pdf'):    
    #Calculate ERP and t-test
    erp1,erp1Ch=ERP(Epoch1)
    erp2,erp2Ch=ERP(Epoch2)
    tvalue,pvalue=tTestERP(Epoch1,Epoch2)
    pmasked=np.ma.masked_where(pvalue > 0.05, pvalue<0.05)*(np.max(tvalue)+0.1)
    pmasked2=np.ma.masked_where(pvalue > 0.01, pvalue<0.01)*(np.max(tvalue)+0.2)
    minimo=np.min((np.min(erp1),np.min(erp2)))
    maximo=np.max((np.max(erp1),np.max(erp2)))
    time=np.linspace(startTime,endTime,int(fs*(endTime-startTime)//1)+1)
    time1=np.linspace(startTime,endTime+0.25,int(fs*(endTime+0.25-startTime)//1)+1)
    xticks=time1[0:-1:int((0.25*fs)//1)]
    if startTime<0:
        zeroPoint=int(-(startTime*fs)//1)
    else:
        zeroPoint=0
    fig,ax = plt.subplots(2,1,figsize=(8,8),gridspec_kw={'height_ratios':[0.8,0.2]})
    ax[0].plot(time,erp2,linewidth=2,label=label2)
    ax[0].plot(time,erp1,linewidth=2,label=label1)
    ax[0].set_xticks(xticks)
    ax[0].set_xticklabels(time1[0:-1:int((0.25*fs)//1)])
    ax[0].plot([time[zeroPoint],time1[zeroPoint]],[minimo,maximo],'k')
    ax[0].set_ylabel('z-score')
    ax[0].legend()
    ax[1].plot(time,tvalue,label='t-value')
    ax[1].plot(time,pmasked,'*',label=r'$p<0.05$')
    ax[1].plot(time,pmasked2,'d',label=r'$p<0.01$')
    ax[1].set_xticks(xticks)
    ax[1].set_xticklabels(time1[0:-1:int((0.25*fs)//1)])
    ax[1].set_ylabel('t-values')
    ax[1].set_xlabel('Time (s)')
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(savefile,dpi=300)
    
def plotITC(Epoch1,Epoch2,fs=100,startTime=-0.5, endTime=2,lowFreq=0.5,highFreq=16,savefile='ITC.pdf'):
    
    fig,ax=plt.subplots(2,2,figsize=(8,8))
    #Epoch1 Stim
    frequencies,times,itcp1,itcl1=ITC(Epoch1,lowFreq=lowFreq,highFreq=highFreq)
    #Epoch2 Sham
    frequencies,times,itcp2,itcl2=ITC(Epoch2,lowFreq=lowFreq,highFreq=highFreq)
    
    ##Time and frequency labels
    labels_time=[]
    nticks=5
    for n in range(0,len(times),nticks):
        labels_time.append('%.2f' % (times[n]+startTime))
    time_ticks=np.arange(0,len(times)+1,nticks)
    
    frequencies_ticks=[]
    for fr in [0.5,1,5,10,15]:
        frequencies_ticks.append(np.argwhere(frequencies>=fr)[0][0])
    labels_frequency=[]
    for m in range(0,len(frequencies_ticks)):
        labels_frequency.append('%.1f' % frequencies[frequencies_ticks[m]])
    norml=colors.Normalize(vmin=0.0,vmax=1.0)
    normp=colors.Normalize(vmin=0.0,vmax=np.pi)
    #ITCP
    im0=ax[0,0].imshow(itcp1,aspect='auto',norm=normp)
    ax[0,0].set_ylabel('Frequency (Hz)')
    ax[0,0].set_xticks(time_ticks)
    ax[0,0].set_xticklabels(labels_time)
    ax[0,0].set_yticks(frequencies_ticks)
    ax[0,0].set_yticklabels(labels_frequency)
    ax[0,0].set_title('Inter-trial phase coherence STIM')
    for ticklabel in ax[0,0].get_xticklabels():
        ticklabel.set_rotation(60)
    cbar0=fig.colorbar(im0, ax=ax[0,0])
    cbar0.set_ticks(np.arange(0,np.pi+0.1,0.1))
    cbar0.set_label('ITC')
    
    im2=ax[0,1].imshow(itcp2,aspect='auto',norm=normp)
    ax[0,1].set_ylabel('Frequency (Hz)')
    ax[0,1].set_xticks(time_ticks)
    ax[0,1].set_xticklabels(labels_time)
    ax[0,1].set_yticks(frequencies_ticks)
    ax[0,1].set_yticklabels(labels_frequency)
    ax[0,1].set_title('Inter-trial phase coherence SHAM')
    for ticklabel in ax[0,1].get_xticklabels():
        ticklabel.set_rotation(60)
    cbar1=fig.colorbar(im2, ax=ax[0,1])
    cbar1.set_label('ITC')
    cbar1.set_ticks(np.arange(0,np.pi+0.1,0.1))
    #ITCL
    im1=ax[1,0].imshow(itcl1,aspect='auto',norm=norml)
    ax[1,0].set_xticks(time_ticks)
    ax[1,0].set_xticklabels(labels_time)
    ax[1,0].set_yticks(frequencies_ticks)
    ax[1,0].set_yticklabels(labels_frequency)
    ax[1,0].set_ylabel('Frequency (Hz)')
    ax[1,0].set_xlabel('Time (s)')
    ax[1,0].set_title('Inter-trial lineal coherence STIM')
    for ticklabel in ax[1,0].get_xticklabels():
        ticklabel.set_rotation(60)
    cbar2=fig.colorbar(im1, ax=ax[1,0])
    cbar2.set_label('ITC')
    cbar2.set_ticks(np.arange(0,1.1,0.1))
    
    im3=ax[1,1].imshow(itcl2,aspect='auto',norm=norml)
    ax[1,1].set_ylabel('Frequency (Hz)')
    ax[1,1].set_xticks(time_ticks)
    ax[1,1].set_xticklabels(labels_time)
    ax[1,1].set_yticks(frequencies_ticks)
    ax[1,1].set_yticklabels(labels_frequency)
    ax[1,1].set_title('Inter-trial lineal coherence SHAM')
    for ticklabel in ax[1,1].get_xticklabels():
        ticklabel.set_rotation(60)
    cbar3=fig.colorbar(im3, ax=ax[1,1])
    cbar3.set_label('ITC')
    cbar3.set_ticks(np.arange(0,1.1,0.1))
    fig.tight_layout()
    fig.savefig(savefile,dpi=300)

def plotWaveletITC(Epoch1,Epoch2,fs=100,startTime=-0.5, endTime=2,lowFreq=0.5,highFreq=16,savefile='waveletITC.pdf'):
    
    fig,ax=plt.subplots(2,2,figsize=(8,8))
    #Epoch1 Stim
    frequencies,times,itcp1,itcl1=waveletITC(Epoch1,lowFreq=lowFreq,highFreq=highFreq)
    #Epoch2 Sham
    frequencies,times,itcp2,itcl2=waveletITC(Epoch2,lowFreq=lowFreq,highFreq=highFreq)
    
    ##Time and frequency labels
    labels_time=[]
    nticks=50
    for n in range(0,len(times),nticks):
        labels_time.append('%.2f' % (times[n]+startTime))
    time_ticks=np.arange(0,len(times)+1,nticks)
    
    frequencies_ticks=[]
    for fr in [0.5,1,5,10,15]:
        frequencies_ticks.append(np.argwhere(frequencies>=fr)[0][0])
    labels_frequency=[]
    for m in range(0,len(frequencies_ticks)):
        labels_frequency.append('%.1f' % frequencies[frequencies_ticks[m]])
    norml=colors.Normalize(vmin=0.0,vmax=1.0)
    normp=colors.Normalize(vmin=0.0,vmax=np.pi)
    #ITCP
    im0=ax[0,0].imshow(itcp1,aspect='auto',norm=normp)
    ax[0,0].set_ylabel('Frequency (Hz)')
    ax[0,0].set_xticks(time_ticks)
    ax[0,0].set_xticklabels(labels_time)
    ax[0,0].set_yticks(frequencies_ticks)
    ax[0,0].set_yticklabels(labels_frequency)
    ax[0,0].set_title('Inter-trial phase coherence STIM')
    for ticklabel in ax[0,0].get_xticklabels():
        ticklabel.set_rotation(60)
    cbar0=fig.colorbar(im0, ax=ax[0,0])
    cbar0.set_ticks(np.arange(0,np.pi+0.1,0.1))
    cbar0.set_label('ITC')
    
    im2=ax[0,1].imshow(itcp2,aspect='auto',norm=normp)
    ax[0,1].set_ylabel('Frequency (Hz)')
    ax[0,1].set_xticks(time_ticks)
    ax[0,1].set_xticklabels(labels_time)
    ax[0,1].set_yticks(frequencies_ticks)
    ax[0,1].set_yticklabels(labels_frequency)
    ax[0,1].set_title('Inter-trial phase coherence SHAM')
    for ticklabel in ax[0,1].get_xticklabels():
        ticklabel.set_rotation(60)
    cbar1=fig.colorbar(im2, ax=ax[0,1])
    cbar1.set_label('ITC')
    cbar1.set_ticks(np.arange(0,np.pi+0.1,0.1))
    #ITCL
    im1=ax[1,0].imshow(itcl1,aspect='auto',norm=norml)
    ax[1,0].set_xticks(time_ticks)
    ax[1,0].set_xticklabels(labels_time)
    ax[1,0].set_yticks(frequencies_ticks)
    ax[1,0].set_yticklabels(labels_frequency)
    ax[1,0].set_ylabel('Frequency (Hz)')
    ax[1,0].set_xlabel('Time (s)')
    ax[1,0].set_title('Inter-trial lineal coherence STIM')
    for ticklabel in ax[1,0].get_xticklabels():
        ticklabel.set_rotation(60)
    cbar2=fig.colorbar(im1, ax=ax[1,0])
    cbar2.set_label('ITC')
    cbar2.set_ticks(np.arange(0,1.1,0.1))
    
    im3=ax[1,1].imshow(itcl2,aspect='auto',norm=norml)
    ax[1,1].set_ylabel('Frequency (Hz)')
    ax[1,1].set_xticks(time_ticks)
    ax[1,1].set_xticklabels(labels_time)
    ax[1,1].set_yticks(frequencies_ticks)
    ax[1,1].set_yticklabels(labels_frequency)
    ax[1,1].set_title('Inter-trial lineal coherence SHAM')
    for ticklabel in ax[1,1].get_xticklabels():
        ticklabel.set_rotation(60)
    cbar3=fig.colorbar(im3, ax=ax[1,1])
    cbar3.set_label('ITC')
    cbar3.set_ticks(np.arange(0,1.1,0.1))
    fig.tight_layout()
    fig.savefig(savefile,dpi=300)

def plotERSP(Epoch1,Epoch2,startTime=-0.5,endTime=2,savefile='ERSP.pdf'):
    frequencies,times,SxxStim,erspStim=ERSP(Epoch1)
    frequencies,times,SxxSham,erspSham=ERSP(Epoch2)
    print('Trials: ')
    print(np.shape(SxxStim)[0])
    print(np.shape(SxxSham)[0])
    #Time and frequency labels
    labels_frequency=[]
    frequencies_ticks=[]
    for fr in [0.5,1,5,10,15,20,30]:
        frequencies_ticks.append(np.argwhere(frequencies>=fr)[0][0])
    for m in range(0,len(frequencies_ticks)):
        labels_frequency.append('%.1f' % frequencies[frequencies_ticks[m]])
    
    labels_time=[]
    nticks=5
    for n in range(0,len(times),nticks):
        labels_time.append('%.2f' % (times[n]+startTime))
    time_ticks=np.arange(0,len(times)+1,nticks)
    
    fig, ax=plt.subplots(3,1,figsize=(8,8))
    im0=ax[0].imshow(10*np.log10(erspStim[0:frequencies_ticks[-1],:]),aspect='auto')
    ax[0].set_xticks(time_ticks)
    ax[0].set_xticklabels(labels_time)
    ax[0].set_yticks(frequencies_ticks)
    ax[0].set_yticklabels(labels_frequency)
    ax[0].set_title('STIM ERSP')
    cbar0=fig.colorbar(im0, ax=ax[0])
    cbar0.set_label('dB')
    
    im1=ax[1].imshow(10*np.log10(erspSham[0:frequencies_ticks[-1],:]),aspect='auto')
    ax[1].set_xticks(time_ticks)
    ax[1].set_xticklabels(labels_time)
    ax[1].set_yticks(frequencies_ticks)
    ax[1].set_yticklabels(labels_frequency)	
    ax[1].set_title('SHAM ERSP')
    cbar1=fig.colorbar(im1, ax=ax[1])		
    cbar1.set_label('dB')

    #Calculate t-Test
    tvalue,pvalue=stats.ttest_ind(SxxStim[:,0:frequencies_ticks[-1],:],SxxSham[:,0:frequencies_ticks[-1],:],axis=0)
    pmasked=np.ma.masked_where(pvalue > 0.05, pvalue<0.05)
    open_pvalue = ndimage.binary_opening(pmasked)
    closed_pvalue = ndimage.binary_closing(open_pvalue)
    im2=ax[2].imshow(np.abs(tvalue),aspect='auto')
    if len(np.nonzero(closed_pvalue))>2:
        ax[2].contour(closed_pvalue,[0.5],linewidths=2, colors='r')
    ax[2].set_xticks(time_ticks)
    ax[2].set_xticklabels(labels_time)
    ax[2].set_yticks(frequencies_ticks)
    ax[2].set_yticklabels(labels_frequency)	
    ax[2].set_title('t-test')
    cbar2=fig.colorbar(im2, ax=ax[2])	
    cbar2.set_label('t-value')
    fig.tight_layout()
    fig.savefig(savefile,dpi=300)

def plotERSPROI(Epoch1,startTime=-0.5,endTime=2):
    frequencies,times,SxxStim,erspStim=waveletERSP(Epoch1)
    timeZero=np.argwhere(np.abs(times+startTime)<1e-6)[0][0]
    timeTwo=np.argwhere((times+startTime)>2)[0][0]
    freqdot5=np.argwhere(frequencies>0.5)[0][0]
    freqone=np.argwhere(frequencies>1.25)[0][0]
    freqeleven=np.argwhere(frequencies>11)[0][0]
    freqsixteen=np.argwhere(frequencies>16)[0][0]
    
    erspRoi1=erspStim[freqdot5:freqone,timeZero:timeTwo]
    erspRoi2=erspStim[freqeleven:freqsixteen,timeZero:timeTwo]
    normRoi=colors.Normalize(vmin=np.min(erspStim),vmax=np.max(erspStim))
    fig,ax=plt.subplots(2,1)
    im1=ax[0].imshow(erspRoi1,aspect='auto',norm=normRoi)
    im2=ax[1].imshow(erspRoi2,aspect='auto',norm=normRoi)
    cbar1=fig.colorbar(im1, ax=ax[0])		
    cbar1.set_label('a. u.')
    cbar2=fig.colorbar(im2, ax=ax[1])		
    cbar2.set_label('a. u.')
    plt.figure()
    plt.plot(np.sum(erspRoi1,axis=0))
    plt.plot(np.sum(erspRoi2,axis=0))
    plt.show()
    print(np.argwhere(erspRoi1==np.max(erspRoi1)))
    print(np.argwhere(erspRoi2==np.max(erspRoi2)))
    print(np.sum(erspRoi1))
    print(np.sum(erspRoi2))

def ERSPROI(Epoch1,startTime=-0.5,endTime=2):
    frequencies,times,SxxStim,erspStim=waveletERSP(Epoch1)
    timeZero=np.argwhere(np.abs(times+startTime)<1e-6)[0][0]
    timeTwo=np.argwhere((times+startTime)>2)[0][0]
    freqdot5=np.argwhere(frequencies>0.5)[0][0]
    freqone=np.argwhere(frequencies>1.25)[0][0]
    freqeleven=np.argwhere(frequencies>11)[0][0]
    freqsixteen=np.argwhere(frequencies>16)[0][0]
    erspRoi1=erspStim[freqdot5:freqone,timeZero:timeTwo]
    erspRoi2=erspStim[freqeleven:freqsixteen,timeZero:timeTwo]
    sumRoi1=np.sum(erspRoi1)
    sumRoi2=np.sum(erspRoi2)
    return sumRoi1,sumRoi2

    
   
def plotWaveletERSP(Epoch1,Epoch2,startTime=-0.5,endTime=2,savefile='waveletERSP.pdf'):
    frequencies,times,SxxStim,erspStim=waveletERSP(Epoch1)
    frequencies,times,SxxSham,erspSham=waveletERSP(Epoch2)
    #Time and frequency labels
    labels_frequency=[]
    frequencies_ticks=[]
    for fr in [0.5,1,5,10,15,20,30]:
        frequencies_ticks.append(np.argwhere(frequencies>=fr)[0][0])
    for m in range(0,len(frequencies_ticks)):
        labels_frequency.append('%.1f' % frequencies[frequencies_ticks[m]])
    
    labels_time=[]
    nticks=50
    
    for n in range(0,len(times),nticks):
        labels_time.append('%.2f' % (times[n]+startTime))
    time_ticks=np.arange(0,len(times)+1,nticks)
    
    
    fig, ax=plt.subplots(3,1,figsize=(8,8))
    im0=ax[0].imshow(erspStim[0:frequencies_ticks[-1],:],aspect='auto')
    ax[0].set_xticks(time_ticks)
    ax[0].set_xticklabels(labels_time)
    ax[0].set_yticks(frequencies_ticks)
    ax[0].set_yticklabels(labels_frequency)
    ax[0].set_title('STIM ERSP')
    cbar0=fig.colorbar(im0, ax=ax[0])
    cbar0.set_label('a. u.')
    
    im1=ax[1].imshow(erspSham[0:frequencies_ticks[-1],:],aspect='auto')
    ax[1].set_xticks(time_ticks)
    ax[1].set_xticklabels(labels_time)
    ax[1].set_yticks(frequencies_ticks)
    ax[1].set_yticklabels(labels_frequency)	
    ax[1].set_title('SHAM ERSP')
    cbar1=fig.colorbar(im1, ax=ax[1])		
    cbar1.set_label('a. u.')

    #Calculate t-Test
    tvalue,pvalue=stats.ttest_ind(SxxStim[:,0:frequencies_ticks[-1],:],SxxSham[:,0:frequencies_ticks[-1],:],axis=0)
    pmasked=np.ma.masked_where(pvalue > 0.05, pvalue<0.05)
    open_pvalue = ndimage.binary_opening(pmasked)
    closed_pvalue = ndimage.binary_closing(open_pvalue)
    im2=ax[2].imshow(np.abs(tvalue),aspect='auto')
    if len(np.nonzero(closed_pvalue))>2:
        ax[2].contour(closed_pvalue,[0.5],linewidths=2, colors='r')
    ax[2].set_xticks(time_ticks)
    ax[2].set_xticklabels(labels_time)
    ax[2].set_yticks(frequencies_ticks)
    ax[2].set_yticklabels(labels_frequency)	
    ax[2].set_title('t-test')
    cbar2=fig.colorbar(im2, ax=ax[2])	
    cbar2.set_label('t-value')
    fig.tight_layout()
    fig.savefig(savefile,dpi=300)
    
def main():
    if len(sys.argv) > 0:
        fs=100
        #Load txt file
        filename = sys.argv[1]
        filenameStim=filename[0:-7]+'Stim.txt'
        indexPath=0		
        for m in range(len(filename)-1,0,-1):
            if filename[m]=='/':
                indexPath=m
                break
        filenameBaseline=filename[0:indexPath+1]+'SWS-baseline-Phi.txt'
        
        data, datazScore =Data.loadStoredData(filename)
        print('Loaded data')
        
        dataBaseline, datazScoreBaseline=Data.loadStoredData(filenameBaseline)
        print('Loaded baseline')
		
        dataStim, markers, time=Data.loadStoredStim(filenameStim)
        print('Loaded stimuli')

        Epochs1, Epochs2, Epochs3,Epochs4,EpochsS1, EpochsS2, EpochsS3,EpochsS4=getEpochs(datazScore, markers,fs=100, shamData=datazScoreBaseline, startTime=-3, endTime=3, baseline=[-2,-1])
        plotERSPROI(Epochs1,startTime=-3, endTime=3)
        plotWaveletERSP(Epochs1,EpochsS1,startTime=-3,endTime=3,savefile=filename[:-4]+'-waveletERSP.pdf')
        plotWaveletITC(Epochs1,EpochsS1,startTime=-3,endTime=3,savefile=filename[:-4]+'-waveletITC.pdf')
        plotERPS(Epochs1,EpochsS1,fs=100,label1='STIM', label2='SHAM', startTime=-3, endTime=3,savefile=filename[:-4]+'-ERP.pdf')   
        #plotITC(Epochs1,EpochsS1,startTime=-2,endTime=3,savefile=filename[:-4]+'-ITC.pdf')
        #plotERSP(Epochs1,EpochsS1,startTime=-2,endTime=3,savefile=filename[:-4]+'-ERSP.pdf')

if __name__=='__main__':
    main()  
