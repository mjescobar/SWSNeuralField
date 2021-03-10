#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 03:15:06 2021

@author: felipe
"""
import sys, os
sys.path.append(os.path.abspath(os.path.join('..', 'plotting')))
from cycler import cycler
import numpy as np
import scipy.stats as stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import statistic_plots as sp
import regression as regression
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
rc('text', usetex=True)


#Load data
file_shapes=np.load('../collectedData/CollectResults_shapesLong.npz')
collect_Iso_shapes=file_shapes['collect_Iso']
collect_Ispindles_shapes=file_shapes['collect_Ispindles']
collectr_Iso_shapes=file_shapes['collectr_Iso']
collectr_Ispindles_shapes=file_shapes['collectr_Ispindles']
collect_pc_shapes=file_shapes['collect_pc']
collect_pcp_shapes=file_shapes['collect_pcp']
collect_pso_shapes=file_shapes['collect_pso']
collect_psp_shapes=file_shapes['collect_psp']
collect_coincidences_shapes=file_shapes['collect_coincidences']
collect_so_shapes=file_shapes['collect_so']
collect_spindles_shapes=file_shapes['collect_spindles']
collectr_pc_shapes=file_shapes['collectr_pc']
collectr_pcp_shapes=file_shapes['collectr_pcp']
collectr_pso_shapes=file_shapes['collectr_pso']
collectr_psp_shapes=file_shapes['collectr_psp']
collectr_coincidences_shapes=file_shapes['collectr_coincidences']
collectr_so_shapes=file_shapes['collectr_so']
collectr_spindles_shapes=file_shapes['collectr_spindles']
y_axis_shapes=[0,1,2,3,4,5]
labels_shapes=['DR','HT','Gauss','Rect','Tri','RR']


file_energy=np.load('../collectedData/CollectResults_energyLong.npz')
collect_Iso_energy=file_energy['collect_Iso']
collect_Ispindles_energy=file_energy['collect_Ispindles']
collectr_Iso_energy=file_energy['collectr_Iso']
collectr_Ispindles_energy=file_energy['collectr_Ispindles']
collect_pc_energy=file_energy['collect_pc']
collect_pcp_energy=file_energy['collect_pcp']
collect_pso_energy=file_energy['collect_pso']
collect_psp_energy=file_energy['collect_psp']
collect_coincidences_energy=file_energy['collect_coincidences']
collect_so_energy=file_energy['collect_so']
collect_spindles_energy=file_energy['collect_spindles']
collectr_pc_energy=file_energy['collectr_pc']
collectr_pcp_energy=file_energy['collectr_pcp']
collectr_pso_energy=file_energy['collectr_pso']
collectr_psp_energy=file_energy['collectr_psp']
collectr_coincidences_energy=file_energy['collectr_coincidences']
collectr_so_energy=file_energy['collectr_so']
collectr_spindles_energy=file_energy['collectr_spindles']
y_axis_energy=np.array([10,20,40,60,80,100])  
labels_energy=['10','20','40','60','80','100']  


file_duration=np.load('../collectedData/CollectResults_durationLong.npz')
collect_Iso_duration=file_duration['collect_Iso']
collect_Ispindles_duration=file_duration['collect_Ispindles']
collectr_Iso_duration=file_duration['collectr_Iso']
collectr_Ispindles_duration=file_duration['collectr_Ispindles']
collect_pc_duration=file_duration['collect_pc']
collect_pcp_duration=file_duration['collect_pcp']
collect_pso_duration=file_duration['collect_pso']
collect_psp_duration=file_duration['collect_psp']
collect_coincidences_duration=file_duration['collect_coincidences']
collect_so_duration=file_duration['collect_so']
collect_spindles_duration=file_duration['collect_spindles']
collectr_pc_duration=file_duration['collectr_pc']
collectr_pcp_duration=file_duration['collectr_pcp']
collectr_pso_duration=file_duration['collectr_pso']
collectr_psp_duration=file_duration['collectr_psp']
collectr_coincidences_duration=file_duration['collectr_coincidences']
collectr_so_duration=file_duration['collectr_so']
collectr_spindles_duration=file_duration['collectr_spindles']
y_axis_duration=np.array([0.05,0.075,0.15,0.30])
labels_duration=['0.05 s','0.07 s','0.15 s','0.30 s']

#Preprare figure layout
fig3=plt.figure(figsize=(5.2,3.75))
gs3=gridspec.GridSpec(2,3,figure=fig3,width_ratios=[0.5,0.5,0.5],height_ratios=[1,1])
gs3.update(wspace=0.38, hspace=0.25) # set the spacing between axes.
axA=fig3.add_subplot(gs3[0,0])
axB=fig3.add_subplot(gs3[0,1])
axC=fig3.add_subplot(gs3[1,0])
axD=fig3.add_subplot(gs3[1,1])
axE=fig3.add_subplot(gs3[0,2])
axF=fig3.add_subplot(gs3[1,2])

#Calculate shapes
A=np.array([10,17.32,17.32,17.32,15.897,15.885,14.883,14.142,16.330,14.142,20.186,30.001,12.248,12.248])
A=A*2
duration=0.1
fs=1000

sigma=duration**2/10 #5 sigma width
#time array
t=np.linspace(0,duration,int(duration*fs+1))
t1=np.linspace(0,duration+2/fs,int(duration*fs+3))
zeroshape=np.zeros((np.shape(t)[0]+2,))
first_half=np.zeros_like(t)
second_half=np.zeros_like(t)
first_quarter=np.zeros_like(t)
last_quarter=np.zeros_like(t)
middle_half=np.zeros_like(t)
first_half[0:int(duration*fs//2)+1]=1
second_half[int(duration*fs//2+1):int(duration*fs+1)]=t[int(duration*fs//2)+1:int(duration*fs+1)]-duration/2
first_quarter[0:int(duration*fs//4+1)]=t[0:int(duration*fs//4+1)]
middle_half[int(duration*fs//4+1):int(duration*fs+1)]=1
last_quarter[int(duration*fs*3//4+1):int(duration*fs+1)]=t[int(duration*fs*3//4+1):int(duration*fs+1)]-3*duration/4
#shapes
rectangular=np.ones_like(t)*A[0]
triangular=A[1]-2/duration*A[1]*np.abs(t-duration/2)
increase_ramp=A[2]*t/duration
decrease_ramp=A[3]*(1-t/duration)
half_trapece=np.ones_like(t)*A[12]-A[12]*(2*second_half/duration)
gaussian=A[4]*np.exp(-(t-duration/2)**2/sigma)


#%%
#collect onset amplitudes
xa=[]
xb=[]
xc=[]
xd=[]
xe=[]
xf=[]
ya=[]
yb=[]
yc=[]
yd=[]
ye=[]
yf=[]


#Plot Shape
axA.plot(decrease_ramp[0],np.mean(collect_Iso_shapes[:,0]),'s',label='1')
axB.plot(decrease_ramp[0],np.mean(collect_so_shapes[:,0]),'s')
axC.plot(decrease_ramp[0],np.mean(collect_Ispindles_shapes[:,0]),'s')
axD.plot(decrease_ramp[0],np.mean(collect_spindles_shapes[:,0]),'s')
axE.plot(decrease_ramp[0],np.mean(collect_pso_shapes[:,0]),'s')
axF.plot(decrease_ramp[0],np.mean(collect_pcp_shapes[:,0]),'s')
xa.append(decrease_ramp[0])
xb.append(decrease_ramp[0])
xc.append(decrease_ramp[0])
xd.append(decrease_ramp[0])
xe.append(decrease_ramp[0])
xf.append(decrease_ramp[0])

ya.append(np.mean(collect_Iso_shapes[:,0]))
yb.append(np.mean(collect_so_shapes[:,0]))
yc.append(np.mean(collect_Ispindles_shapes[:,0]))
yd.append(np.mean(collect_spindles_shapes[:,0]))
ye.append(np.mean(collect_pso_shapes[:,0]))
yf.append(np.mean(collect_pcp_shapes[:,0]))

axA.plot(half_trapece[0],np.mean(collect_Iso_shapes[:,1]),'s',label='2')
axB.plot(half_trapece[0],np.mean(collect_so_shapes[:,1]),'s')
axC.plot(half_trapece[0],np.mean(collect_Ispindles_shapes[:,1]),'s')
axD.plot(half_trapece[0],np.mean(collect_spindles_shapes[:,1]),'s')
axE.plot(half_trapece[0],np.mean(collect_pso_shapes[:,1]),'s')
axF.plot(half_trapece[0],np.mean(collect_pcp_shapes[:,1]),'s')
xa.append(half_trapece[0])
xb.append(half_trapece[0])
xc.append(half_trapece[0])
xd.append(half_trapece[0])
xe.append(half_trapece[0])
xf.append(half_trapece[0])

ya.append(np.mean(collect_Iso_shapes[:,1]))
yb.append(np.mean(collect_so_shapes[:,1]))
yc.append(np.mean(collect_Ispindles_shapes[:,1]))
yd.append(np.mean(collect_spindles_shapes[:,1]))
ye.append(np.mean(collect_pso_shapes[:,1]))
yf.append(np.mean(collect_pcp_shapes[:,1]))


axA.plot(gaussian[0],np.mean(collect_Iso_shapes[:,2]),'s',label='3')
axB.plot(gaussian[0],np.mean(collect_so_shapes[:,2]),'s')
axC.plot(gaussian[0],np.mean(collect_Ispindles_shapes[:,2]),'s')
axD.plot(gaussian[0],np.mean(collect_spindles_shapes[:,2]),'s')
axE.plot(gaussian[0],np.mean(collect_pso_shapes[:,2]),'s')
axF.plot(gaussian[0],np.mean(collect_pcp_shapes[:,2]),'s')
xa.append(gaussian[0])
xb.append(gaussian[0])
xc.append(gaussian[0])
xd.append(gaussian[0])
xe.append(gaussian[0])
xf.append(gaussian[0])

ya.append(np.mean(collect_Iso_shapes[:,2]))
yb.append(np.mean(collect_so_shapes[:,2]))
yc.append(np.mean(collect_Ispindles_shapes[:,2]))
yd.append(np.mean(collect_spindles_shapes[:,2]))
ye.append(np.mean(collect_pso_shapes[:,2]))
yf.append(np.mean(collect_pcp_shapes[:,2]))


axA.plot(rectangular[0],np.mean(collect_Iso_shapes[:,3]),'s',label='4')
axB.plot(rectangular[0],np.mean(collect_so_shapes[:,3]),'s')
axC.plot(rectangular[0],np.mean(collect_Ispindles_shapes[:,3]),'s')
axD.plot(rectangular[0],np.mean(collect_spindles_shapes[:,3]),'s')
axE.plot(rectangular[0],np.mean(collect_pso_shapes[:,3]),'s')
axF.plot(rectangular[0],np.mean(collect_pcp_shapes[:,3]),'s')
xa.append(rectangular[0])
xb.append(rectangular[0])
xc.append(rectangular[0])
xd.append(rectangular[0])
xe.append(rectangular[0])
xf.append(rectangular[0])

ya.append(np.mean(collect_Iso_shapes[:,3]))
yb.append(np.mean(collect_so_shapes[:,3]))
yc.append(np.mean(collect_Ispindles_shapes[:,3]))
yd.append(np.mean(collect_spindles_shapes[:,3]))
ye.append(np.mean(collect_pso_shapes[:,3]))
yf.append(np.mean(collect_pcp_shapes[:,3]))


axA.plot(triangular[0],np.mean(collect_Iso_shapes[:,4]),'s',label='5')
axB.plot(triangular[0],np.mean(collect_so_shapes[:,4]),'s')
axC.plot(triangular[0],np.mean(collect_Ispindles_shapes[:,4]),'s')
axD.plot(triangular[0],np.mean(collect_spindles_shapes[:,4]),'s')
axE.plot(triangular[0],np.mean(collect_pso_shapes[:,4]),'s')
axF.plot(triangular[0],np.mean(collect_pcp_shapes[:,4]),'s')
xa.append(triangular[0])
xb.append(triangular[0])
xc.append(triangular[0])
xd.append(triangular[0])
xe.append(triangular[0])
xf.append(triangular[0])

ya.append(np.mean(collect_Iso_shapes[:,4]))
yb.append(np.mean(collect_so_shapes[:,4]))
yc.append(np.mean(collect_Ispindles_shapes[:,4]))
yd.append(np.mean(collect_spindles_shapes[:,4]))
ye.append(np.mean(collect_pso_shapes[:,4]))
yf.append(np.mean(collect_pcp_shapes[:,4]))


axA.plot(increase_ramp[0],np.mean(collect_Iso_shapes[:,5]),'s',label='6')
axB.plot(increase_ramp[0],np.mean(collect_so_shapes[:,5]),'s')
axC.plot(increase_ramp[0],np.mean(collect_Ispindles_shapes[:,5]),'s')
axD.plot(increase_ramp[0],np.mean(collect_spindles_shapes[:,5]),'s')
axE.plot(increase_ramp[0],np.mean(collect_pso_shapes[:,5]),'s')
axF.plot(increase_ramp[0],np.mean(collect_pcp_shapes[:,5]),'s')
xa.append(increase_ramp[0])
xb.append(increase_ramp[0])
xc.append(increase_ramp[0])
xd.append(increase_ramp[0])
xe.append(increase_ramp[0])
xf.append(increase_ramp[0])

ya.append(np.mean(collect_Iso_shapes[:,5]))
yb.append(np.mean(collect_so_shapes[:,5]))
yc.append(np.mean(collect_Ispindles_shapes[:,5]))
yd.append(np.mean(collect_spindles_shapes[:,5]))
ye.append(np.mean(collect_pso_shapes[:,5]))
yf.append(np.mean(collect_pcp_shapes[:,5]))

#plot Energy
duration=0.1
second_half=np.zeros_like(t)
second_half[int(duration*fs//2+1):int(duration*fs+1)]=t[int(duration*fs//2)+1:int(duration*fs+1)]-duration/2
for n,energy in enumerate(y_axis_energy):
    An=np.sqrt(energy/duration)*A[2]/20.0
    decrease_ramp=An*(1-t/duration)
    zeroshape[1:-1]=decrease_ramp
    # axA2.plot(t1,zeroshape,label=labels_energy[n],color=plt.cm.Reds((1.0-n*.12)))
    axA.plot(decrease_ramp[0],np.mean(collect_Iso_energy[:,n]),'o',color=plt.cm.Blues((1.0-n*.12)),label=labels_energy[n])
    axB.plot(decrease_ramp[0],np.mean(collect_so_energy[:,n]),'o',color=plt.cm.Blues((1.0-n*.12)))
    axC.plot(decrease_ramp[0],np.mean(collect_Ispindles_energy[:,n]),'o',color=plt.cm.Blues((1.0-n*.12)))
    axD.plot(decrease_ramp[0],np.mean(collect_spindles_energy[:,n]),'o',color=plt.cm.Blues((1.0-n*.12)))
    axE.plot(decrease_ramp[0],np.mean(collect_pso_energy[:,n]),'o',color=plt.cm.Blues((1.0-n*.12)))
    axF.plot(decrease_ramp[0],np.mean(collect_pcp_energy[:,n]),'o',color=plt.cm.Blues((1.0-n*.12)))
    xa.append(decrease_ramp[0])
    xb.append(decrease_ramp[0])
    xc.append(decrease_ramp[0])
    xd.append(decrease_ramp[0])
    xe.append(decrease_ramp[0])
    xf.append(decrease_ramp[0])
    
    ya.append(np.mean(collect_Iso_energy[:,n]))
    yb.append(np.mean(collect_so_energy[:,n]))
    yc.append(np.mean(collect_Ispindles_energy[:,n]))
    yd.append(np.mean(collect_spindles_energy[:,n]))
    ye.append(np.mean(collect_pso_energy[:,n]))
    yf.append(np.mean(collect_pcp_energy[:,n]))

#plot Duration
for n,duration in enumerate(y_axis_duration):
    An=np.sqrt(40/duration)
    t=np.linspace(0,duration,int(duration*fs+1))
    t1=np.linspace(0,duration+2/fs,int(duration*fs+3))
    zeroshape=np.zeros((np.shape(t)[0]+2,))
    rectangular=np.ones_like(t)*An
    zeroshape[1:-1]=rectangular
    # axA3.plot(t1,zeroshape,label=labels_duration[n],color=plt.cm.Greens((1.0-n*.12)))
    axA.plot(rectangular[0],np.mean(collect_Iso_duration[:,n]),'o',color=plt.cm.Reds((1.0-n*.12)),label=labels_duration[n])
    axB.plot(rectangular[0],np.mean(collect_so_duration[:,n]),'o',color=plt.cm.Reds((1.0-n*.12)))
    axC.plot(rectangular[0],np.mean(collect_Ispindles_duration[:,n]),'o',color=plt.cm.Reds((1.0-n*.12)))
    axD.plot(rectangular[0],np.mean(collect_spindles_duration[:,n]),'o',color=plt.cm.Reds((1.0-n*.12)))
    axE.plot(rectangular[0],np.mean(collect_pso_duration[:,n]),'o',color=plt.cm.Reds((1.0-n*.12)))
    axF.plot(rectangular[0],np.mean(collect_pcp_duration[:,n]),'o',color=plt.cm.Reds((1.0-n*.12)))
    xa.append(rectangular[0])
    xb.append(rectangular[0])
    xc.append(rectangular[0])
    xd.append(rectangular[0])
    xe.append(rectangular[0])
    xf.append(rectangular[0])
    
    ya.append(np.mean(collect_Iso_shapes[:,n]))
    yb.append(np.mean(collect_so_shapes[:,n]))
    yc.append(np.mean(collect_Ispindles_shapes[:,n]))
    yd.append(np.mean(collect_spindles_shapes[:,n]))
    ye.append(np.mean(collect_pso_shapes[:,n]))
    yf.append(np.mean(collect_pcp_shapes[:,n]))    


xa=np.array(xa)
xb=np.array(xb)
xc=np.array(xc)
xd=np.array(xd)
xe=np.array(xe)
xf=np.array(xf)
ya=np.array(ya)
yb=np.array(yb)
yc=np.array(yc)
yd=np.array(yd)
ye=np.array(ye)
yf=np.array(yf)
x1=np.arange(0,61)

wa,erra,pra=regression.regression(xa, ya, 1)
ya1=regression.prediction(x1, wa)
axA.plot(x1,ya1,'k')
axA.set_ylabel('$I^{(SO)}$',fontsize=8,labelpad=0)
axA.text(0,0.4,'r:%.3f'%pra)
axA.text(-0.35,1,'A',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axA.transAxes)
axA.set_xlim([-1,61])
axA.set_ylim([-0.05,0.6])
axA.tick_params(axis='both', which='major', labelsize=8)

wb,errb,prb=regression.regression(xb, yb, 1)
yb1=regression.prediction(x1, wb)
axB.plot(x1,yb1,'k')
axB.set_ylabel('$N_{SO}$/min.',fontsize=8,labelpad=0)
axB.text(0,20,'r:%.3f'%prb)
axB.text(-0.35,1,'B',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axB.transAxes)
axB.set_xlim([-1,61])
axB.set_ylim([10,24])
axB.tick_params(axis='both', which='major', labelsize=8)

wc,errc,prc=regression.regression(xc, yc, 1)
yc1=regression.prediction(x1, wc)
axC.plot(x1,yc1,'k')
axC.set_ylabel('$I^{(SP)}$',fontsize=8,labelpad=0)
axC.text(0,0.2,'r:%.3f'%prc)
axC.text(-0.35,1,'D',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axC.transAxes)
axC.set_xlabel('Pulse onset amplitude',fontsize=8,labelpad=0)
axC.set_xlim([-1,61])
axC.set_ylim([-0.05,0.3])
axC.tick_params(axis='both', which='major', labelsize=8)

wd,errd,prd=regression.regression(xd, yd, 1)
yd1=regression.prediction(x1, wd)
axD.plot(x1,yd1,'k')
axD.set_ylabel('$N_{SP}$/min.',fontsize=8,labelpad=0)
axD.text(0,12,'r:%.3f'%prd)
axD.text(-0.35,1,'E',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axD.transAxes)
axD.set_xlabel('Pulse onset amplitude',fontsize=8,labelpad=0)
axD.set_xlim([-1,61])
axD.set_ylim([5,15])
axD.tick_params(axis='both', which='major', labelsize=8)

we,erre,pre=regression.regression(xe, ye, 1)
ye1=regression.prediction(x1, we)
axE.plot(x1,ye1,'k')
axE.set_ylabel('$P(C)$',fontsize=8,labelpad=0)
axE.text(0,0.3,'r:%.3f'%pre)
axE.text(-0.35,1,'C',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axE.transAxes)
axE.set_xlim([-1,61])
axE.set_ylim([0.1,0.4])
axE.tick_params(axis='both', which='major', labelsize=8)

wf,errf,prf=regression.regression(xf, yf, 1)
yf1=regression.prediction(x1, wf)
axF.plot(x1,yf1,'k')
axF.set_ylabel('$P(C|SP)$',fontsize=8,labelpad=0)
axF.set_xlabel('Pulse onset amplitude',fontsize=8,labelpad=0)
axF.text(0,0.4,'r: %.3f'%prf)
axF.text(-0.35,1,'F',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axF.transAxes)
axF.set_xlim([-1,61])
axF.set_ylim([0.1,0.6])
axF.tick_params(axis='both', which='major', labelsize=8)

hA,lA=axA.get_legend_handles_labels()
legfig3=fig3.legend(hA,lA,ncol=8,fontsize=8,loc='upper left',bbox_to_anchor=(0.07,0.95,0.1,0.1),columnspacing=1.1,handletextpad=0.3)

fig3.savefig('./output/S6_Fig.eps',dpi=300,bbox_inches='tight',bbox_extra_artists=[legfig3])