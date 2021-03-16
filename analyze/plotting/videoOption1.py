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
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import pandas as pd
import Data
import sys
from matplotlib import rc
rc('text', usetex=True)

def handle_close(evt):
        sys.exit()
class VData(object):
    def __init__(self,data,norm,lims,fs):
        self.data=data
        self.norm=norm
        self.minimo=lims[0]
        self.maximo=lims[1]
        self.fs=fs
        self.fig, self.ax = plt.subplots(figsize=(6,6))
        gs1=gridspec.GridSpec(1,1,figure=self.fig,width_ratios=[1],height_ratios=[1])
        self.ax=self.fig.add_subplot(gs1[0])
    def make_frame(self,t):
        n=int(t*self.fs)
        self.ax.clear()
        self.ax.imshow(np.reshape(self.data[n,:],(16,16)).T,cmap=plt.cm.RdYlBu,interpolation='gaussian',norm=self.norm)
        self.ax.text(0.5,0.5,'z-score interval = $\pm$%.2f'%self.maximo)
        self.ax.set_axis_off()
        plt.draw()
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
        anim.write_videofile('SWS.mp4', fps=fps)
        a.fig.canvas.mpl_connect('close_event',handle_close)
    else:
        print("Error: filename not valid")
if __name__=='__main__':
    main()        

