#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 00:11:08 2019
@author: felipe
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.signal as signal
import scipy.stats as stats
import Data
import sys
from matplotlib import rc
rc('text', usetex=True)

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise(ValueError, "smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise(ValueError, "Input vector needs to be bigger than window size.")

    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise(ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len//2-1):-(window_len//2)]

def detectBand(meanData,fs=100,lowFrequency=9,highFrequency=16,threshold=0.26,tmin=0.45,gap_max=0.1,windowLength=0.1):
    b,a=signal.cheby1(4,1e-6,[lowFrequency/(fs/2), highFrequency/(fs/2)],'bandpass')
    filtered=signal.filtfilt(b,a,meanData)
    rms=np.zeros((len(filtered),1))
    smoothed_rms=np.zeros((len(filtered),1))
    band_detection=np.zeros((len(filtered),1))
    timeWindow=int(np.ceil(windowLength*fs))
    #Spindle Characterisitics
    tmin=tmin*fs
    gap_max=gap_max*fs
    hold=0
    start=0
    accum=0
    gap=0
    count=0
    flag_filt=0
    for n in range(timeWindow,len(filtered)):
        rms[n]=np.sqrt(np.sum((filtered[n-timeWindow:n])**2)/timeWindow)
        if n>timeWindow:
            smooth_rms=smooth(rms[n-timeWindow:n,0],window_len=timeWindow,window='hamming')
            if flag_filt==0:
                smoothed_rms[0:timeWindow,0]=smooth_rms
                flag_filt=1
            else:
                smoothed_rms[n]=smooth_rms[timeWindow-1]
            if smoothed_rms[n]>=threshold:
                if hold==0:
                    start=n
                accum+=1
                hold=1
            elif hold==1:
                accum+=1
                gap+=1
                if gap>gap_max:
                    gap=0
                    hold=0
                    if accum>=tmin:
                        band_detection[start:n]=1
                        count+=1
                    accum=0
    return filtered,smoothed_rms,band_detection,count

def detectSpindles(meanData,fs=100,threshold=0.26,gap_max=0.1,time_Window=0.2):
    b,a=signal.cheby1(4,1e-6,[9/(fs/2), 16/(fs/2)],'bandpass')
    filtered=signal.filtfilt(b,a,meanData)
    rms=np.zeros((len(filtered),1))
    smoothed_rms=np.zeros((len(filtered),1))
    spindle_detection=np.zeros((len(filtered),1))
    timeWindow=int(np.ceil(time_Window*fs));
    #Spindle Characterisitics
    tmin=0.49*fs
    tmax=2.01*fs
    gap_max=gap_max*fs #Inter-spindle refractoriness
    hold=0
    start=0
    accum=0
    gap=0
    count=0
    flag_filt=0
    for n in range(timeWindow,len(filtered)):
        rms[n]=np.sqrt(np.sum((filtered[n-timeWindow:n])**2)/timeWindow)
        if n>timeWindow:
            smooth_rms=smooth(rms[n-timeWindow:n,0],window_len=timeWindow,window='hamming')
            if flag_filt==0:
                smoothed_rms[0:timeWindow,0]=smooth_rms
                flag_filt=1
            else:
                smoothed_rms[n,0]=smooth_rms[timeWindow-1]
            if smoothed_rms[n]>=threshold:
                if hold==0:
                    start=n
                accum+=1
                hold=1
            elif hold==1:
                accum+=1
                gap+=1
                if gap>gap_max:
                    gap=0
                    hold=0
                    if accum>=tmin and accum<tmax:
                        spindle_detection[start:n]=1
                        count+=1
                    accum=0
    return filtered,smoothed_rms,spindle_detection,count

def detectSpindlesSTALTA(meanData,fs=100,threshold=0.26,gap_max=0.1,short_time=0.05,long_time=0.2,time_smooth=0.2):
    b,a=signal.cheby1(4,1e-6,[9/(fs/2), 16/(fs/2)],'bandpass')
    filtered=signal.filtfilt(b,a,meanData)
    stalta=np.zeros((len(filtered),1))
    smoothed_rms=np.zeros((len(filtered),1))
    spindle_detection=np.zeros((len(filtered),1))
    timeWindow=int(np.ceil(time_smooth*fs))
    stime=int(np.ceil(short_time*fs))
    ltime=int(np.ceil(long_time*fs))
    #Spindle Characterisitics
    tmin=0.49*fs
    tmax=2.01*fs
    gap_max=gap_max*fs #Inter-spindle refractoriness
    hold=0
    start=0
    accum=0
    gap=0
    count=0
    flag_filt=0
    for n in range(ltime,len(filtered)):
        rms=np.sqrt(np.sum(filtered[n-stime:n]**2)*ltime/(np.sum(filtered[n-ltime:n]**2)*stime))
        stalta[n]=np.log(1+rms**2)
        if n>stime:
            smooth_rms=smooth(stalta[n-timeWindow:n,0],window_len=timeWindow,window='hamming')
            if flag_filt==0:
                smoothed_rms[0:timeWindow,0]=smooth_rms
                flag_filt=1
            else:
                smoothed_rms[n,0]=smooth_rms[timeWindow-1]
            if smoothed_rms[n]>=threshold:
                if hold==0:
                    start=n
                accum+=1
                hold=1
            elif hold==1:
                accum+=1
                gap+=1
                if gap>gap_max:
                    gap=0
                    hold=0
                    if accum>=tmin and accum<tmax:
                        spindle_detection[start:n]=1
                        count+=1
                    accum=0
    return filtered,smoothed_rms,spindle_detection,count

def detectSpindlesSTALTA2(meanData,fs=100,threshold=0.26,gap_max=0.1,short_time=0.05,long_time=0.2,time_smooth=0.2):
    b,a=signal.cheby1(4,1e-6,[9/(fs/2), 16/(fs/2)],'bandpass')
    filtered=signal.filtfilt(b,a,meanData)
    logrms=np.zeros((len(filtered),1))
    smoothed_rms=np.zeros((len(filtered),1))
    spindle_detection=np.zeros((len(filtered),1))
    timeWindow=int(np.ceil(time_smooth*fs))
    stime=int(np.ceil(short_time*fs))
    ltime=int(np.ceil(long_time*fs))
    #Spindle Characterisitics
    tmin=0.49*fs
    tmax=2.01*fs
    gap_max=gap_max*fs #Inter-spindle refractoriness
    hold=0
    start=0
    accum=0
    gap=0
    count=0
    flag_filt=0
    for n in range(ltime,len(filtered)):
        # rms=np.sqrt(np.sum((filtered[n-timeWindow:n])**2)/timeWindow)
        rms=np.sqrt(np.sum(filtered[n-stime:n]**2)*ltime/(np.sum(filtered[n-ltime:n]**2)*stime))
        logrms[n]=1/(1+np.exp(-rms))
        if n>stime:
            smooth_rms=smooth(logrms[n-timeWindow:n,0],window_len=timeWindow,window='hamming')
            if flag_filt==0:
                smoothed_rms[0:timeWindow,0]=smooth_rms
                flag_filt=1
            else:
                smoothed_rms[n,0]=smooth_rms[timeWindow-1]
            if smoothed_rms[n]>=threshold:
                if hold==0:
                    start=n
                accum+=1
                hold=1
            elif hold==1:
                accum+=1
                gap+=1
                if gap>gap_max:
                    gap=0
                    hold=0
                    if accum>=tmin and accum<tmax:
                        spindle_detection[start:n]=1
                        count+=1
                    accum=0
    return filtered,smoothed_rms,spindle_detection,count

def detectsegmentSpindles(meanData,fs=100,threshold=0.26,gap_max=0.1):
    b,a=signal.cheby1(4,1e-6,[9/(fs/2), 16/(fs/2)],'bandpass')
    filtered=signal.filtfilt(b,a,meanData)
    rms=np.zeros((len(filtered),1))
    smoothed_rms=np.zeros((len(filtered),1))
    spindle_detection=np.zeros((len(filtered),1))
    spindle_segmentation=np.zeros((len(filtered),1))
    timeWindow=int(np.ceil(0.2*fs));
    #Spindle Characterisitics
    tmin=0.49*fs
    tmax=2.01*fs
    gap_max=gap_max*fs; #Inter-spindle refractoriness
    hold=0
    start=0
    accum=0
    gap=0
    count=0
    flag_filt=0
    for n in range(timeWindow,len(filtered)):
        rms[n]=np.sqrt(np.sum((filtered[n-timeWindow:n])**2)/timeWindow)
        if n>timeWindow:
            smooth_rms=smooth(rms[n-timeWindow:n,0],window_len=timeWindow,window='hamming')
            if flag_filt==0:
                smoothed_rms[0:timeWindow,0]=smooth_rms
                flag_filt=1
            else:
                smoothed_rms[n]=smooth_rms[timeWindow-1]
            if smoothed_rms[n]>=threshold:
                if hold==0:
                    start=n
                accum+=1
                hold=1
            elif hold==1:
                accum+=1
                gap+=1
                if gap>gap_max:
                    gap=0
                    hold=0
                    if accum>=tmin and accum<tmax:
                        spindle_detection[start:n]=1
                        count+=1
                        spindle_segmentation[start:n]=count
                    accum=0
    return filtered,smoothed_rms,spindle_detection,spindle_segmentation,count

def meanRMS(meanData,fs=100,lowFrequency=9,highFrequency=16):
    b,a=signal.cheby1(4,1e-6,[lowFrequency/(fs/2), highFrequency/(fs/2)],'bandpass')
    filtered=signal.filtfilt(b,a,meanData)
    rms=np.zeros((len(filtered),1))
    smoothed_rms=np.zeros((len(filtered),1))
    spindle_detection=np.zeros((len(filtered),1))
    timeWindow=int(np.ceil(0.2*fs));
    #Spindle Characterisitics
    flag_filt=0
    for n in range(timeWindow,len(filtered)):
        rms[n]=np.sqrt(np.sum((filtered[n-timeWindow:n])**2)/timeWindow)
        if n>timeWindow:
            smooth_rms=smooth(rms[n-timeWindow:n,0],window_len=timeWindow,window='hamming')
            if flag_filt==0:
                smoothed_rms[0:timeWindow,0]=smooth_rms
                flag_filt=1
            else:
                smoothed_rms[n]=smooth_rms[timeWindow-1]
    return np.mean(smoothed_rms)

def plotDetection(meanData,detection,fs=100,label='Spindles',filename='spindles'):
    plt.figure()
    plt.plot(np.arange(0,np.shape(meanData)[0])*1/fs, meanData, label=r'$\phi_e$')
    plt.plot(np.arange(0,np.shape(meanData)[0])*1/fs,detection, label=label)
    plt.ylabel('z-score')
    plt.xlabel('Time (s)')
    plt.legend()
    plt.savefig(filename+'-detection.pdf',dpi=300)

def eventsxMinute(segmentedSignal,fs=100):
    events_samples=np.argwhere(segmentedSignal!=0)
    minutes=1
    previous_events=0
    current_events=0
    number_events=np.zeros((int(len(segmentedSignal)//(60*fs)),))
    percentage_events=np.zeros((int(len(segmentedSignal)//(60*fs)),))
    for n in range(len(segmentedSignal)):
        if n in events_samples:
            current_events=segmentedSignal[n]
        if n>60*minutes*fs:
            number_events[minutes-1]=current_events-previous_events
            percentage_events[minutes-1]=np.count_nonzero(segmentedSignal[60*(minutes-1)*fs:60*minutes*fs])/(60.0*fs)
            previous_events=current_events
            minutes+=1
    return number_events,percentage_events
    
def main():
    if len(sys.argv) > 0:
        fs=100
        #Load txt file
        filename = sys.argv[1]
        print('Loading and plotting '+filename+'...')
        data, datazScore=Data.loadStoredData(filename)
        filtered, rms,detectionSpindles=detectSpindles(datazScore,fs=fs)
        plotDetection(datazScore,detectionSpindles,filename=filename[:-4])
        filtered, rms,detectionDelta=detectBand(datazScore,fs=fs,lowFrequency=4.5,highFrequency=8.0,threshold=.26,tmin=0.45,gapMax=0.05,windowLength=0.1)
        plotDetection(datazScore,detectionDelta, label='Theta',filename=filename[:-4]+'-theta')
        

if __name__=='__main__':
    main()  