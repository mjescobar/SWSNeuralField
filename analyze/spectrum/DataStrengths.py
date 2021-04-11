#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 16:43:38 2020

@author: felipe
"""
import pandas as pd
import numpy as np
import scipy.stats as stats
import scipy.signal as signal
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
import Wavelets
import sys

def loadStoredStrengths(filename,sep='\t'):
    strengths_index=np.arange(0,8)
    phi_index=np.arange(8,13)
    ratios_index=np.arange(13,21)
    textfile = open(filename,'r')
    data=pd.read_csv(textfile, sep=sep)
    data=data.to_numpy()
    strengths=np.array(data[:,strengths_index],dtype='float')
    phi=np.array(data[:,phi_index],dtype='float')
    ratios=np.array(data[:,ratios_index],dtype='float')
    return strengths, phi, ratios

def calculateRho(phi,Qmax=340,sigma_rho=0.0038):
    rho=phi[:,:]/sigma_rho*(np.ones_like(phi)-(phi[:,:]/Qmax))
    return rho

def calculateGains(strengths,phi):
    #strengths: nu_ee,nu_ei,nu_es,nu_re,nu_rs,nu_se,nu_sr,nu_sn
    #phi:  phi_e, phi_i, phi_r, phi_s, phi_n
    rho=calculateRho(phi)
    G=np.zeros((np.shape(strengths)[0],8))
    G[:,0]=strengths[:,0]*rho[:,0]#Gee
    G[:,1]=strengths[:,1]*rho[:,0]#Gei
    G[:,2]=strengths[:,2]*rho[:,0]#Ges
    G[:,3]=strengths[:,3]*rho[:,2]#Gre
    G[:,4]=strengths[:,4]*rho[:,2]#Grs
    G[:,5]=strengths[:,5]*rho[:,3]#Gse
    G[:,6]=strengths[:,6]*rho[:,3]#Gsr
    G[:,7]=strengths[:,7]*rho[:,3]#Gsn
    return G

def calculateXYZ(G,alpha=45,beta=186):
    #X=Gee/(1-Gei)
    X=G[:,0]/(np.ones_like(G[:,0])-G[:,1])
    #Y=(Ges*Gse+Ges*Gsr*Gre)/((1-Gsr*Grs)*(1-Gei))
    Y=(G[:,2]*G[:,5]+G[:,2]*G[:,6]*G[:,3])/((np.ones_like(G[:,6])-G[:,6]*G[:,4])*(np.ones_like(G[:,1])-G[:,1]))
    #Z=-Gsr*Grs*(alpha*beta)/(alpha+beta)^2
    Z=-G[:,6]*G[:,4]*(alpha*beta)/(alpha+beta)**2
    return X,Y,Z
    

    
# strengths, phi, ratios =loadStoredStrengths('/home/felipe/PlasticityandAAS/SWSAASPhi01-seed-1-baseline-Strengths.txt')
# plt.plot(strengths)
# plt.figure()
# plt.plot(phi)
# plt.show()
# G=calculateGains(strengths,phi)
# X,Y,Z=calculateXYZ(G)

# plt.figure()
# plt.plot(G)
# #%%
# plt.figure()
# for m in range(0,80,2):
#     plt.plot(X[m*10:(m+1)*10],Y[m*10:(m+1)*10],'b')
#     plt.plot(X[(m+1)*10:(m+2)*10],Y[(m+1)*10:(m+2)*10],'k')
# plt.figure()    
# for m in range(0,80,2):
#     plt.plot(1-X[m*10:(m+1)*10],Z[m*10:(m+1)*10],'r')
#     plt.plot(1-X[(m+1)*10:(m+2)*10],Z[(m+1)*10:(m+2)*10],'g')
#     plt.plot(Y[m*10:(m+1)*10],Z[m*10:(m+1)*10],'m')
#     plt.plot(Y[(m+1)*10:(m+2)*10],Z[(m+1)*10:(m+2)*10],'c')
# plt.plot(Y,Z)
