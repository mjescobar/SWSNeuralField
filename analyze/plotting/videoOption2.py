#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 02:46:50 2019

@author: felipe
"""
import sys, os
sys.path.append(os.path.abspath('../spectrum'))
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import pandas as pd
import Data
import scipy.signal as signal
from matplotlib import rc
rc('text', usetex=True)

def handle_close(evt):
        sys.exit()
class VData(object):
    def __init__(self,data,norm,lims,fs):
        
        self.meanData=np.mean(data,axis=1)
        self.data=np.reshape(data,(np.shape(data)[0],16,16))
        self.norm=norm
        self.minimo=lims[0]
        self.maximo=lims[1]
        self.fs=fs
        self.fig = plt.figure(figsize=(6,6),frameon=False)
        gs1=gridspec.GridSpec(2,1,figure=self.fig,width_ratios=[1],height_ratios=[1,0.2])
        self.ax=self.fig.add_subplot(gs1[0],frame_on=False)
        self.ax1=self.fig.add_subplot(gs1[1],frame_on=True)
        clb=plt.colorbar(cm.ScalarMappable(norm=self.norm, cmap=plt.cm.RdYlBu),ax=self.ax)
        clb.set_label('z-score')
    def make_frame(self,t):
        n=int(t*self.fs)
        t_start=np.floor(t/6)
        n_second=int(t_start*6*self.fs)
        n_now=n-n_second
        tarray=np.linspace(0,6,6*self.fs+1)
        self.ax.clear()
        self.ax1.clear()
        self.ax.tick_params(axis='both',which='both',length=0)
        self.ax.imshow(self.data[n,:,:].T,cmap=plt.cm.RdYlBu,interpolation='bilinear',norm=self.norm)
        self.ax.set_xticklabels([])
        self.ax.set_yticklabels([])
        self.ax1.plot(tarray[0:n_now],self.meanData[n_second:n_second+n_now])
        self.ax1.set_xlabel('time (s)')
        self.ax1.set_ylabel('z-score')
        self.ax1.set_ylim([self.minimo,self.maximo])
        self.ax1.set_xlim([0,6])
        self.ax1.set_xticklabels(np.arange(t_start*6,t_start*6+6))
        if n==0:
            plt.show()
        return mplfig_to_npimage(self.fig)
        
    
def main():
    if len(sys.argv) > 0:
        fs=100
        fps=60
        #Load txt file
        filename = sys.argv[1]
        print('Loading and plotting '+filename+' ...')
        data, datazScore,std=Data.loadStoredData(filename)
        data=datazScore.to_numpy(dtype=float)
        lowFrequency=0.1
        highFrequency=40
        b,a=signal.cheby1(4,1e-6,[lowFrequency/(fs/2), highFrequency/(fs/2)],'bandpass')
        data=signal.filtfilt(b,a,data,axis=0)
        minimo=np.min(data)
        maximo=np.max(data)
        if minimo<-maximo:
            minimo=minimo//2*2
            maximo=-minimo
        else:
            maximo=maximo//2*2
            minimo=-maximo
        lims=(minimo,maximo)
        norm= colors.Normalize(vmin=minimo, vmax=maximo)
        a=VData(data,norm,lims,fs)        
        time=np.shape(data)[0]
        print(time)
        anim = VideoClip(a.make_frame, duration=int((time//fs)))
        anim.write_videofile('SWSOption2_bilinear.mp4', fps=fps)
        a.fig.canvas.mpl_connect('close_event',handle_close)
    else:
        print("Error: filename not valid")
if __name__=='__main__':
    main()        

