#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 17:07:52 2019

@author: felipe
"""
import numpy as np


def Dtemporal(omega,alpha,beta):
    D=(1-omega*1j/alpha)*(1-omega*1j/beta);
    return D

def CharacteristicEquation(omega,k,alpha,beta,gamma,r,rho,strengths,t0):
    Gsrs=rho[2]*rho[1]*strengths[6]*strengths[4];
    Gee=rho[0]*strengths[0];
    Gei=rho[0]*strengths[1];
    Gese=rho[0]*rho[2]*strengths[2]*strengths[5];
    Gesre=rho[0]*rho[1]*rho[2]*strengths[2]*strengths[6]*strengths[3];
    L=1/Dtemporal(omega,alpha,beta);
    q2r2=(1-omega*1j/gamma)**2-(1/(1-Gei*L))*(Gee*L+((Gese*L**2+Gesre*L**3)*np.exp(1j*omega*t0))/(1-Gsrs*L**2))
    s=(1-Gei*L)*(1-Gsrs*L**2)*(q2r2+k**2*r**2);
    return s
    
def TransferFunction(omega,k,alpha,beta,gamma,r,rho,strengths,t0):
    Ges=rho[0]*strengths[2];
    Gsn=rho[2]*strengths[7];
    L=1/Dtemporal(omega,alpha,beta);
    denominator=CharacteristicEquation(omega,k,alpha,beta,gamma,r,rho,strengths,t0);
    numerator=Ges*Gsn*(L**2)*np.exp(1j*omega*t0/2);
    h=numerator/denominator
    h=h/np.sum(numerator)
    return h

    
def SpatialSpectrum(omega,Lx,Ly,alpha,beta,gamma,r,rho,strengths,t0,k0):
    p_a=0
    for m in range(25):
        for n in range(25):
             kx=2*np.pi*m/Lx
             ky=2*np.pi*n/Ly
             k=np.sqrt(kx**2+ky**2);
             p_a+=np.abs(TransferFunction(omega,k,alpha,beta,gamma,r,rho,strengths,t0))**2*np.exp(-k**2/k0)*2*np.pi/Lx*2*np.pi/Ly;
    P=p_a
    return P