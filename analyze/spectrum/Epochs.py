#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 18:36:57 2019

@author: felipe
"""
import numpy as np

class Epochs:
    def __init__(self,markerType,fs,startTime,endTime):
        self.stim=markerType
        self.fs=fs
        self.data=np.array([(),()])
        self.time=np.linspace(startTime,endTime,(endTime-startTime)*self.fs+1)
        
    def addData(self,newData):
        if len(self.data)==2:
            self.data=newData
        else:
            self.data=np.hstack((self.data,newData))
        
    def baselineCorrect(self,baseline):
        zero=np.argwhere(np.abs(self.time)<1/self.fs)[0][0]
        init=int(np.ceil(baseline[0]*self.fs)+zero)
        ending=int(np.ceil(baseline[1]*self.fs)+zero)
        if len(np.shape(self.data))>1: 
            meanCorrection=np.nanmean(self.data[init:ending,:])
        else:
            meanCorrection=np.nanmean(self.data[init:ending])
        self.data=self.data-meanCorrection
        
    def setTime(self,startTime,endTime):
        self.time=np.linspace(startTime,endTime,(endTime-startTime)*self.fs+1)