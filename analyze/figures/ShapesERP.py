#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 14:40:13 2021

@author: felipe
"""
import sys, os
sys.path.append(os.path.abspath('../spectrum'))
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import Data
import ERP
import Epochs
import scipy.signal as signal
import matplotlib.gridspec as gridspec
from matplotlib import rc
rc('text', usetex=True)

erps=np.zeros((301,6))
channelerps=np.zeros((301,6))
channelerps1=np.zeros((301,6))
# irs=np.zeros((999,6))


#Without noise
#[0,1,2,3,4,13];
#['rectangular','rising ramp','decreasing ramp','triangular','Gaussian','rectangular trapezoid']
# for n,nshape in enumerate([3,6,5,1,4,2]):
#     filename='/home/felipe/Dropbox/UTFSM/NeuralField/timeseries/shapes/ShapePulseLowNoise-%d-Phi.txt'%nshape
#     data, datazScore,std=Data.loadStoredData(filename)
#     data=data.to_numpy(dtype=float)
#     data=np.mean(data,axis=1)-np.mean(data)
#     irs[:,n]=data
#     del data


#With noise, ERP
fs=100
names=['decreaseRamp','halfTrapece','Gaussian','rectangular','triangle','riseRamp']
amplitudes=['3.46e+01','2.45e+01','3.18e+01','2.00e+01','3.46e+01','3.46e+01']
for n in range(6):
    filename='/media/felipe/TOSHIBA2T/Shapes-Long/SWSlong-seed-1-%s-F0.20-A%s-D1.00e-01-sdA-OFF-sdF-OFF-pF0.02-pR2.00-Phi.txt'%(names[n],amplitudes[n])
    filenameStim='/media/felipe/TOSHIBA2T/Shapes-Long/SWSlong-seed-1-%s-F0.20-A%s-D1.00e-01-sdA-OFF-sdF-OFF-pF0.02-pR2.00-Stim.txt'%(names[n],amplitudes[n])
    data, datazScore,std=Data.loadStoredData(filename,L=np.arange(0,256))
    data=data.to_numpy(dtype=float)
    datasingle=np.mean(data,axis=1)-np.mean(data)
    dataStim, markers, time=Data.loadStoredStim(filenameStim,L=256)
    Epochs1all, Epochs2, Epochs3,Epochs4=ERP.getEpochs(data,markers,fs,shamData=0,startTime=-1,endTime=2,baseline=[-0.45,0])
    Epochs1, Epochs2, Epochs3,Epochs4=ERP.getEpochs(datasingle,markers,fs,shamData=0,startTime=-1,endTime=2,baseline=[-0.45,0])
    allerp,allerpChannels=ERP.ERP(Epochs1all,fs=100,NxNy=256,lowFreq=0.1,highFreq=40)
    erp=ERP.ERPSingle(Epochs1,fs=100,lenTime=len(Epochs1.time),lowFreq=0.1,highFreq=40)
    channelerps[:,n]=np.mean(allerpChannels[:,218:222],axis=1)
    channelerps1[:,n]=np.mean(allerpChannels[:,219],axis=1)
    erps[:,n]=erp
    del data, datazScore, std, Epochs1all, Epochs1
#%%
fig=plt.figure(figsize=(5.2,3.7))
gs1 = gridspec.GridSpec(1, 2, figure=fig, 
                        height_ratios=[1],
                        width_ratios=[1,1],wspace=0.5,hspace=0.4)
axA=fig.add_subplot(gs1[0])
axB=fig.add_subplot(gs1[1])

time=np.arange(-0.1,0.5,0.01)

axA.plot(time,channelerps[90:150,0],color=plt.cm.tab10(0),label='1')
# axC.plot(time,channelerps[80:180,0],color=plt.cm.tab10(0),label='1' )

# axC.plot(time,channelerps[80:180,1],color=plt.cm.tab10(1),label='2' )
axA.plot(time,channelerps[90:150,1],color=plt.cm.tab10(1),label='2')

# axC.plot(time,channelerps[80:180,2],color=plt.cm.tab10(2),label='3' )
axA.plot(time,channelerps[90:150,2],color=plt.cm.tab10(2),label='3')

# axD.plot(time,channelerps[80:180,3],color=plt.cm.tab10(3),label='4' )
axB.plot(time,channelerps[90:150,3],color=plt.cm.tab10(3),label='4')

# axD.plot(time,channelerps[80:180,4],color=plt.cm.tab10(4),label='5' )
axB.plot(time,channelerps[90:150,4],color=plt.cm.tab10(4),label='5')

# axD.plot(time,channelerps[80:180,5],color=plt.cm.tab10(5),label='6' )
axB.plot(time,channelerps[90:150,5],color=plt.cm.tab10(5),label='6')


# axA.plot(time,channelerps1[90:150,0],color=plt.cm.tab10(0),label='1')

# axA.plot(time,channelerps1[90:150,1],color=plt.cm.tab10(1),label='2')

# axA.plot(time,channelerps1[90:150,2],color=plt.cm.tab10(2),label='3')

# axB.plot(time,channelerps1[90:150,3],color=plt.cm.tab10(3),label='4')

# axB.plot(time,channelerps1[90:150,4],color=plt.cm.tab10(4),label='5')

# axB.plot(time,channelerps1[90:150,5],color=plt.cm.tab10(5),label='6')


axA.set_xlabel('Time (s)',fontsize=8)
axB.set_xlabel('Time (s)',fontsize=8)
# axC.set_xlabel('Time (s)',fontsize=8)
# axD.set_xlabel('Time (s)',fontsize=8)

axA.set_ylabel('ERP ($s^{-1}$)',fontsize=8)
axB.set_ylabel('ERP ($s^{-1}$)',fontsize=8)
# axC.set_ylabel('ERP $s^{-1}$',fontsize=8)
# axD.set_ylabel('ERP $s^{-1}$',fontsize=8)


axA.set_xticks([-0.1,0,0.1,0.2,0.3,0.4,0.5])
axB.set_xticks([-0.1,0,0.1,0.2,0.3,0.4,0.5])
# axC.set_xticks([-0.2,0,0.2,0.4,0.6,0.8])
# axD.set_xticks([-0.2,0,0.2,0.4,0.6,0.8])
axA.set_xticklabels([-0.1,0,0.1,0.2,0.3,0.4,0.5])
axB.set_xticklabels([-0.1,0,0.1,0.2,0.3,0.4,0.5])
# axC.set_xticklabels([-0.2,0,0.2,0.4,0.6,0.8])
# axD.set_xticklabels([-0.2,0,0.2,0.4,0.6,0.8])

axA.set_ylim([-0.02,0.01])
axB.set_ylim([-0.02,0.01])
# axC.set_ylim([-0.011,0.061])
# axD.set_ylim([-0.011,0.031])


hA,lA=axA.get_legend_handles_labels()
hB,lB=axB.get_legend_handles_labels()
for h,l in zip(hB,lB):
    hA.append(h)
    lA.append(l)

axA.vlines([0,0.1],ymin=-0.02,ymax=0.01,linestyle=[(0,(1,0)),(1,(2,3))],color='k')
axB.vlines([0,0.1],ymin=-0.02,ymax=0.01,linestyle=[(0,(1,0)),(1,(2,3))],color='k')
axA.legend(hA,lA,ncol=6, loc='upper left',fontsize=8,bbox_to_anchor=(-0.05,1.01,1,0.1),columnspacing=2.9,handletextpad=0.9)

axA.text(-0.3,1.05,'A',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axA.transAxes)
axB.text(-0.3,1.05,'B',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axB.transAxes)
# axC.text(-0.2,1,'C',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axC.transAxes)
# axD.text(-0.2,1,'D',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axD.transAxes)


fig.savefig('FigS1.eps',dpi=300,bbox_inches='tight')

