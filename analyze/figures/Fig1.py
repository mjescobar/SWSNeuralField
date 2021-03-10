#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 4 19:10:51 2019

@author: felipe
"""
import sys, os
sys.path.append(os.path.abspath(os.path.join('..', 'plotting')))
import numpy as np
import scipy.io as sio
import h5py
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.patches as patches
from matplotlib import rc
rc('text', usetex=True)

import steadyFunctions

#%%%
#interpolation point
phi_ni=np.zeros((101,))
phi_ei=np.zeros((101,))
phi_ri=np.zeros((101,))
phi_si=np.zeros((101,))
rho_ei=np.zeros((101,))
rho_ri=np.zeros((101,))
rho_si=np.zeros((101,))
Geei=np.zeros((101,))
Geii=np.zeros((101,))
Gesi=np.zeros((101,))
Grei=np.zeros((101,))
Grsi=np.zeros((101,))
Gsei=np.zeros((101,))
Gsri=np.zeros((101,))
Gsni=np.zeros((101,))
Gsrsi=np.zeros((101,))
Gesei=np.zeros((101,))
Gesrei=np.zeros((101,))
XYZi=np.zeros((3,101))
#Selected point
phin=np.arange(1e-15,2.1,0.1)

phi_n=np.zeros_like(phin)
phi_e=np.zeros_like(phin)
phi_r=np.zeros_like(phin)
phi_s=np.zeros_like(phin)
rho_e=np.zeros_like(phin)
rho_r=np.zeros_like(phin)
rho_s=np.zeros_like(phin)
Gee=np.zeros_like(phin)
Gei=np.zeros_like(phin)
Ges=np.zeros_like(phin)
Gre=np.zeros_like(phin)
Grs=np.zeros_like(phin)
Gse=np.zeros_like(phin)
Gsr=np.zeros_like(phin)
Gsn=np.zeros_like(phin)
Gsrs=np.zeros_like(phin)
Gese=np.zeros_like(phin)
Gesre=np.zeros_like(phin)
XYZ=np.zeros((3,len(phin)))

Qmax=340
sigma_p=0.0038
nus=np.array([5.541, -5.652, 1.53, 0.286, 1.12, 2.67, -1.73, 9.22])*1e-3
nu1=np.array([3.06e-3, -3.24e-3, 0.92e-3, 0.26e-3, 2.88e-3, 4.73e-3, -1.95e-3, 2.70e-3])
nu2=np.array([6.81e-3, -6.89e-3, 1.85e-3, 0.3e-3, 0.21e-3, 1.61e-3, -1.62e-3, 12.6e-3])

alpha=45
beta=186
gamma=116
r=0.086
k0=10
t0=0.085
Lx=0.5
Ly=0.5

omega=np.arange(0,2*np.pi*50+0.001,np.pi/100)
Power=np.zeros((len(omega),len(phin)))
Poweri=np.zeros((len(omega),101))

#%% 
print("Load fsolve data")
for step in range(1,102):
    #h5py requires the absolute path
    matfile=h5py.File('/home/felipe/Dropbox/UTFSM/NeuralField/analyze/collectedData/NF-interpolation/interpolation_step_'+str(step)+'-phin1.mat','r')
    solutions_varios=np.array(matfile['soulutions_e'])
    nusi=np.array(matfile['nus'])
    if step==66:
        solutions_e=solutions_varios

    #Phin=0
    index_phin=np.argwhere(solutions_varios[:,0]>0)[0][0]
    phi_ni[step-1]=solutions_varios[index_phin,0]
    phi_ei[step-1]=solutions_varios[index_phin,1]
    phi_ri[step-1]=solutions_varios[index_phin,2]
    phi_si[step-1]=solutions_varios[index_phin,3]
    rho_ei[step-1]=phi_ei[step-1]/sigma_p*(1-(phi_ei[step-1]/Qmax));
    rho_ri[step-1]=phi_ri[step-1]/sigma_p*(1-(phi_ri[step-1]/Qmax));
    rho_si[step-1]=phi_si[step-1]/sigma_p*(1-(phi_si[step-1]/Qmax));
    Geei[step-1]=nusi[0,step-1]*rho_ei[step-1]
    Geii[step-1]=nusi[1,step-1]*rho_ei[step-1]
    Gesi[step-1]=nusi[2,step-1]*rho_ei[step-1]
    Grei[step-1]=nusi[3,step-1]*rho_ri[step-1]
    Grsi[step-1]=nusi[4,step-1]*rho_ri[step-1]
    Gsei[step-1]=nusi[5,step-1]*rho_si[step-1]
    Gsri[step-1]=nusi[6,step-1]*rho_si[step-1]
    Gsni[step-1]=nusi[7,step-1]*rho_si[step-1]
    Gesrei[step-1]=Gesi[step-1]*Gsri[step-1]*Grei[step-1]
    Gesei[step-1]=Gesi[step-1]*Gsei[step-1]
    Gsrsi[step-1]=Gsri[step-1]*Grsi[step-1]
    rho=np.array([[rho_ei[step-1]],[rho_ri[step-1]],[rho_si[step-1]]])
    XYZi[:,step-1]=np.array([Geei[step-1]/(1-Geii[step-1]), (Gesi[step-1]*Gsei[step-1]+Gesi[step-1]*Gsri[step-1]*Grei[step-1])/((1-Gsri[step-1]*Grsi[step-1])*(1-Geii[step-1])),-Gsri[step-1]*Grsi[step-1]*(alpha*beta)/(alpha+beta)**2])
    Poweri[:,step-1]=steadyFunctions.SpatialSpectrum(omega,Lx,Ly,alpha,beta,gamma,r,rho,nusi[:,step-1],t0,k0)   

#%%
print("Calculate spectrum with phi_n variation")
for value,n in zip(phin,range(len(phin))):
    index_phin=np.argwhere(solutions_e[:,0]>=value)[0][0]
    phi_n[n]=solutions_e[index_phin,0]
    phi_e[n]=solutions_e[index_phin,1]
    phi_r[n]=solutions_e[index_phin,2]
    phi_s[n]=solutions_e[index_phin,3]
    rho_e[n]=phi_e[n]/sigma_p*(1-(phi_e[n]/Qmax));
    rho_r[n]=phi_r[n]/sigma_p*(1-(phi_r[n]/Qmax));
    rho_s[n]=phi_s[n]/sigma_p*(1-(phi_s[n]/Qmax));
    Gee[n]=nus[0]*rho_e[n]
    Gei[n]=nus[1]*rho_e[n]
    Ges[n]=nus[2]*rho_e[n]
    Gre[n]=nus[3]*rho_r[n]
    Grs[n]=nus[4]*rho_r[n]
    Gse[n]=nus[5]*rho_s[n]
    Gsr[n]=nus[6]*rho_s[n]
    Gsn[n]=nus[7]*rho_s[n]
    Gesre[n]=Ges[n]*Gsr[n]*Gre[n]
    Gese[n]=Ges[n]*Gse[n]
    Gsrs[n]=Gsr[n]*Grs[n]
    rho=np.array([[rho_e[n]],[rho_r[n]],[rho_s[n]]])
    XYZ[:,n]=np.array([Gee[n]/(1-Gei[n]), (Ges[n]*Gse[n]+Ges[n]*Gsr[n]*Gre[n])/((1-Gsr[n]*Grs[n])*(1-Gei[n])),-Gsr[n]*Grs[n]*(alpha*beta)/(alpha+beta)**2])
    Power[:,n]=steadyFunctions.SpatialSpectrum(omega,Lx,Ly,alpha,beta,gamma,r,rho,nus,t0,k0)

#%%
#Plots 
fig=plt.figure(figsize=(5.2,5.2))
gs=gridspec.GridSpec(2,2,figure=fig,width_ratios=[1,1],height_ratios=[1,1],wspace=0.35, hspace=0.3)

ax1=fig.add_subplot(gs[1,0])
ax2=fig.add_subplot(gs[0,1])
ax4=fig.add_subplot(gs[1,1])
correctedPower=np.zeros((np.shape(Power)[0]-1,np.shape(Power)[1]))
#Spectrums
for n in range(0,len(phin),2):
    correctedPower[:,n]=Power[1::,n]*omega[1::]/(2*np.pi)
    if n==0:
        ax1.loglog(omega[20:6001]/(2*np.pi),correctedPower[20:6001,n],color='blue',label=r'$\phi_n^{(0)}=0$')
    elif n==10:
        ax1.loglog(omega[20:6001]/(2*np.pi),correctedPower[20:6001,n],'k',linewidth=3,label='$\phi_n^{(0)}=1$')
    elif n==len(phin)-1:
        ax1.loglog(omega[20:6001]/(2*np.pi),correctedPower[20:6001,n],color='red',label='$\phi_n^{(0)}=2$')
    elif n==11:
        ax1.loglog(omega[20:6001]/(2*np.pi),correctedPower[20:6001,n],color='#cccccc',label='$Other\ \phi_n^{(0)}$')
    else:
        ax1.loglog(omega[20:6001]/(2*np.pi),correctedPower[20:6001,n],color='#cccccc')

ax1.loglog([0.85,0.85],[7e-4,7e-4],'^',color='black',linewidth=2,label='f=0.85 Hz')
ax1.loglog([0.5,1.25],[7e-4,7e-4],linestyle=(0, (3, 2)),color='black',linewidth=2)
ax1.loglog([9,16],[7e-4,7e-4],linestyle=(0, (2, 1)),color='black',linewidth=2)
ax1.text(0.7,1e-3,'SO',fontsize=8)
ax1.text(10,1e-3,'SP',fontsize=8)
ax1.set_xlim([0.1,30])
ax1.set_xlabel('Frequency (Hz)',fontsize=8)
ax1.set_ylim([3e-8,2e-3])
ax1.set_ylabel(r'Power ($s^{-2}$/Hz $\times$ Hz)',fontsize=8)
ax1.legend(ncol=2,fontsize=8,loc='lower left',handletextpad=0.5,handlelength=0.5,columnspacing=0.3,borderaxespad=0.1)
ax1.tick_params(labelsize=8)
ax1.text(-0.1,1,'C',fontsize=10,fontweight=1000,verticalalignment='bottom',transform=ax1.transAxes)

#

ax2.loglog(omega[20:6001]/(2*np.pi),Poweri[20:6001,0]*omega[20:6001]/(2*np.pi),color='C1',linewidth=2)
ax2.loglog(omega[20:6001]/(2*np.pi),Poweri[20:6001,100]*omega[20:6001]/(2*np.pi),color='C2',linewidth=2)
ax2.loglog(omega[20:6001]/(2*np.pi),Poweri[20:6001,50]*omega[20:6001]/(2*np.pi),color='magenta',linewidth=2)
ax2.loglog(omega[20:6001]/(2*np.pi),Poweri[20:6001,66]*omega[20:6001]/(2*np.pi),color='blue',linewidth=3)
ax2.loglog(omega[20:6001]/(2*np.pi),Poweri[20:6001,80]*omega[20:6001]/(2*np.pi),color='cyan',linewidth=2)
ax2.loglog([0.85,0.85],[7e-4,7e-4],'^',color='black',linewidth=2)
ax2.loglog([0.5,1.25],[7e-4,7e-4],linestyle=(0, (3, 2)),color='black',linewidth=2)
ax2.loglog([11,16],[7e-4,7e-4],linestyle=(0, (2, 1)),color='black',linewidth=2)
ax2.text(0.7,1e-3,'SO',fontsize=8)
ax2.text(10,1e-3,'SP',fontsize=8)
ax2.set_xlim([0.1,30])
ax2.set_ylim([5e-8,2e-3])
#ax2.set_xticklabels([1,10],fontsize=8)
#ax2.set_yticklabels([1e-5,1e-4,1e-3],fontsize=8)
ax2.legend(['N2','N3','50\%','66\%','80\%','f=0.85 Hz'], fontsize=8,ncol=2, loc='lower left',handletextpad=0.5,handlelength=0.5,borderaxespad=0.1)
ax2.set_xlabel('Frequency (Hz)',fontsize=8)
ax2.set_ylabel(r'Power ($s^{-2}/Hz \times Hz$)',fontsize=8)
ax2.tick_params(labelsize=8)
ax2.text(-0.2,1,'B',fontsize=10,fontweight=1000,verticalalignment='bottom',transform=ax2.transAxes)



# 3d plot 

ax3=fig.add_subplot(gs[0,0], projection='3d')
ax3.set_frame_on(False)
# ax3.set_axis_off()
#Phin
ax3.scatter(XYZ[0,:], XYZ[1,:], XYZ[2,:], s=6,marker='o',c='blue')
ax3.scatter(XYZi[0,:], XYZi[1,:], XYZi[2,:], s=6,marker='o')

#Phin
ax3.scatter(XYZ[0,10], XYZ[1,10], XYZ[2,10], marker='d',s=80, c='black',edgecolors='black')
ax3.scatter(XYZ[0,20], XYZ[1,20], XYZ[2,20], marker='d',s=40, c='red',edgecolors='red')

#Interpolation
ax3.scatter(XYZi[0,0], XYZi[1,0], XYZi[2,0], marker='d',s=40, c='C1',edgecolors='C1')
ax3.scatter(XYZi[0,100], XYZi[1,100], XYZi[2,100], marker='d',s=40, c='C2',edgecolors='C2')
ax3.scatter(XYZi[0,50], XYZi[1,50], XYZi[2,50], marker='d',s=40, c='magenta',edgecolors='magenta')

ax3.scatter(XYZi[0,66], XYZi[1,66], XYZi[2,66], marker='d',s=40, c='blue',edgecolors='blue')
ax3.scatter(XYZi[0,80], XYZi[1,80], XYZi[2,80], marker='d',s=40, c='cyan',edgecolors='cyan')
ax3.text(XYZi[0,0], XYZi[1,0]+0.007, XYZi[2,0], 'N2', None,fontsize=8)
ax3.text(XYZi[0,100], XYZi[1,100]+0.007, XYZi[2,100]-0.05, 'N3', None,fontsize=8)
ax3.text(XYZ[0,0], XYZ[1,0]+0.007, XYZ[2,0]+0.007, r'$\phi_n^{(0)}=0$', None,fontsize=8)
ax3.text(XYZ[0,10], XYZ[1,10]+0.007, XYZ[2,10]+0.007, r'$\phi_n^{(0)}=1$', None,fontsize=8)
ax3.text(XYZ[0,20], XYZ[1,20]+0.007, XYZ[2,20]+0.007, r'$\phi_n^{(0)}=2$', None,fontsize=8)
ax3.set_xlabel('X',fontsize=8)
ax3.set_ylabel('Y',fontsize=8)
ax3.set_zlabel('Z',fontsize=8,labelpad=1)
ax3.set_xlim([0.84, 0.94])
ax3.set_ylim([-0.04, 0.04])
ax3.set_zlim([0, 1.0])

ax3.view_init(30, 45)
ax3.tick_params(labelsize=8)
for tick in ax3.get_xticklabels():
    tick.set_rotation(-35)
for tick in ax3.get_yticklabels():
    tick.set_rotation(30)
ax1.text(-0.2,2.31,'A',fontsize=10,fontweight=1000,verticalalignment='bottom',transform=ax1.transAxes)



#Gains
color=next(ax2._get_lines.prop_cycler)['color']
ax4.plot(phin,(Gee-Gee[10])/Gee[10]*100,label=r'$G_{ee}$',color=color)
ax4.plot(1,0,'d',color='k')
color=next(ax4._get_lines.prop_cycler)['color']
ax4.plot(phin,(Gei-Gei[10])/Gei[10]*100,label=r'$-G_{ei}$',color=color)
# ax4.plot(1,(Gei[10]-Gei[10])/Gei[10]*100,'d',color=color)
color=next(ax4._get_lines.prop_cycler)['color']
ax4.plot(phin,(Gese-Gese[10])/Gese[10]*100,label=r'$G_{es}G_{se}$',color=color)
# ax4.plot(1,(Gese[10]-Gese[10])/Gese[10]*100,'d',color=color)
color=next(ax4._get_lines.prop_cycler)['color']
ax4.plot(phin,(Gesre-Gesre[10])/Gesre[10]*100,label=r'$-G_{es}G_{re}G_{sr}$',color=color)
# ax4.plot(1,(Gesre[10]-Gesre[10])/Gesre[10]*100,'d',color=color)
color=next(ax4._get_lines.prop_cycler)['color']
ax4.plot(phin,(Gsrs-Gsrs[10])/Gsrs[10]*100,label=r'$-G_{sr}G_{rs}$',color=color)
# ax4.plot(1,(Gsrs[10]-Gsrs[10])/Gsrs[10]*100,'d',color=color)

color=next(ax4._get_lines.prop_cycler)['color']
ax4.plot(phin,(phi_e-phi_e[10])/phi_e[10]*100,'--',label=r'$\phi_e^{(0)}$',color=color)
# ax4.plot(1,(phi_e[10]-phi_e[10])/phi_e[10]*100,'^',color=color)
color=next(ax4._get_lines.prop_cycler)['color']
ax4.plot(phin,(phi_r-phi_r[10])/phi_r[10]*100,'--',label=r'$\phi_r^{(0)}$',color=color)
# ax4.plot(1,(phi_r[10]-phi_r[10])/phi_r[10]*100,'^',color=color)
color=next(ax4._get_lines.prop_cycler)['color']
ax4.plot(phin,(phi_s-phi_s[10])/phi_s[10]*100,'--',label=r'$\phi_s^{(0)}$',color=color)
# ax4.plot(1,(phi_s[10]-phi_s[10])/phi_s[10]*100,'^',color=color)

ax4.set_xlabel(r'$\phi_n^{(0)}$ $(s^{-1})$',labelpad=0.1,fontsize=8)
ax4.set_ylabel('Loop Gains (\%)',labelpad=0.1,fontsize=8)
ax4.set_ylim([-101,151])
ax4.set_xlim([-0.01,2.01])
ax4.set_xticks([0,0.5,1,1.5,2.0])
ax4.set_xticklabels([r'0.0',r'0.5',r'1.0',r'1.5',r'2.0'],fontsize=8)
ax4.set_yticks([-100,-50,0,50,100,150])
ax4.set_yticklabels(['-100','-50','0','50','100','150'],fontsize=8)

ax4.tick_params(labelsize=8)
hB4,lB4=ax4.get_legend_handles_labels()
hB4A=[]
hB4B=[]
lB4A=[]
lB4B=[]
for i,h in enumerate(hB4):
    if i<5:
        hB4A.append(h)
    else:
        hB4B.append(h)
for i,l in enumerate(lB4):
    if i<5: 
        lB4A.append(l)
    else:
        lB4B.append(l)
legA=ax4.legend(hB4A,lB4A,fontsize=8,ncol=1,handletextpad=0.5,handlelength=0.5,loc='upper left',borderaxespad=0.1)
legB=ax4.legend(hB4B,lB4B,fontsize=8,ncol=1,handletextpad=0.5,handlelength=0.5,loc='lower right',borderaxespad=0.1)
ax4.add_artist(legA)
ax4.add_artist(legB)
ax4.text(-0.2,1,'D',fontsize=10,fontweight=1000,verticalalignment='bottom',transform=ax4.transAxes)
fig.savefig('./output/Fig1.eps',dpi=300)
