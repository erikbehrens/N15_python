#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 21:24:35 2019

@author: behrense
"""

from netCDF4 import Dataset 
import matplotlib.pyplot as plt
import numpy as np


def initial_conditions(din):
    rn15std=0.0036765
    tmp=1.005*rn15std*din/(1+1.005*rn15std)
    return tmp

def delta(C,C15):
    rn15std=0.0036765
    tmp=(C15/(C-C15)/rn15std-1)*1000
    
    return tmp

d1=Dataset('DIN.nc')
lon=d1.variables['nav_lon'][:]
lat=d1.variables['nav_lat'][:]
DIN=d1.variables['DIN'][:]
d3=Dataset('tmask.nc')
tmask=d3.variables['tmaskutil'][:]


#%% create initial conditions for tracer concentrations as per NEMO code

PHN=np.zeros(DIN.shape)
PHN[:14,:,:]=0.1 # over the frist 200 m else 0
PHD=PHN.copy()
ZMI=PHN.copy()
ZME=PHN.copy()
DET=PHN.copy()

DIN_N15=initial_conditions(DIN)
PHN_N15=initial_conditions(PHN)
PHD_N15=initial_conditions(PHD)
ZMI_N15=initial_conditions(ZMI)
ZME_N15=initial_conditions(ZME)
DET_N15=initial_conditions(DET)

#%%  initial tracer concentrations next time step
PHN_next=np.zeros(PHN.shape)
PHD_next=np.zeros(PHN.shape)
ZMI_next=np.zeros(PHN.shape)
ZME_next=np.zeros(PHN.shape)
DIN_next=np.zeros(PHN.shape)
DET_next=np.zeros(PHN.shape)

PHN_N15_next=np.zeros(PHN.shape)
PHD_N15_next=np.zeros(PHN.shape)
ZMI_N15_next=np.zeros(PHN.shape)
ZME_N15_next=np.zeros(PHN.shape)
DIN_N15_next=np.zeros(PHN.shape)
DET_N15_next=np.zeros(PHN.shape)


#%% check that things are ok
plt.close('all')
fig = plt.figure(figsize=(15,10))
ax1 = fig.add_subplot(231)
plt.pcolormesh(np.where(tmask[0,:,:]==0,1,delta(PHN[0,0,:,:],PHN_N15[0,0,:,:])),vmin=4.9,vmax=5.1)
plt.colorbar()
plt.text(10,10,'PHN',color='white')
ax1 = fig.add_subplot(232)
plt.pcolormesh(np.where(tmask[0,:,:]==0,1,delta(PHD[0,0,:,:],PHD_N15[0,0,:,:])),vmin=4.9,vmax=5.1)
plt.colorbar()
plt.text(10,10,'PHD',color='white')
ax1 = fig.add_subplot(233)
plt.pcolormesh(np.where(tmask[0,:,:]==0,1,delta(ZMI[0,0,:,:],ZMI_N15[0,0,:,:])),vmin=4.9,vmax=5.1)
plt.colorbar()
plt.text(10,10,'ZMI',color='white')
ax1 = fig.add_subplot(234)
plt.pcolormesh(np.where(tmask[0,:,:]==0,1,delta(ZME[0,0,:,:],ZME_N15[0,0,:,:])),vmin=4.9,vmax=5.1)
plt.colorbar()
plt.text(10,10,'ZME',color='white')
ax1 = fig.add_subplot(235)
plt.pcolormesh(np.where(tmask[0,:,:]==0,1,delta(DIN[0,0,:,:],DIN_N15[0,0,:,:])),vmin=4.9,vmax=5.1)
plt.colorbar()
plt.text(10,10,'DIN',color='white')
ax1 = fig.add_subplot(236)
plt.pcolormesh(np.where(tmask[0,:,:]==0,1,delta(DET[0,0,:,:],DET_N15[0,0,:,:])),vmin=4.9,vmax=5.1)
plt.colorbar()
plt.text(10,10,'DET',color='white')


plt.savefig('Initial_surface_concentration.png',dpi=200,pading=.1,bbox_inches='tight')

#%% code as similar to Fortran

def calculate_ratios():  
    global fprn,zphn,fprd,zphd
    global DIN,DIN_N15,PHN,PHN_N15,PHD,PHD_N15,ZMI,ZMI_N15,ZME,ZME_N15
    global rn15std
    global un,rdin,bassimphn,fcassimphn,bassimphd,fcassimphn,rtdin15,rtphnn15,rtphdn15
    global rtzmin15,rtzmen15,rtdetrn15,rmi,bmiexcr,fcmiexcr,rme,bmeexcr,fcmeexcr

    eps_assim = 0.
    eps_excr= 0.

    rmin=0.0000001
    rmax=0.9999999
    rdt=40*60
   
    #keep code as similar to Fortran so we can copy across
    
    zdim=0
    tdim=0
    
    for jj in range(fprn.shape[0]):
        for ji in range(fprn.shape[1]):
              #uno3 = npp*dtbio/biono3
              un[jj,ji]=(fprn[jj,ji] * zphn[jj,ji] + fprd[jj,ji] * zphd[jj,ji]) * rdt /max(DIN[tdim,zdim,jj,ji],rmin)
              un[jj,ji] = min(un[jj,ji], rmax)
              un[jj,ji] = max(un[jj,ji], rmin)
              #rno3 = biodin15/(biono3-biodin15)
              rdin[jj,ji] = DIN_N15[tdim,zdim,jj,ji]/(max(DIN[tdim,zdim,jj,ji]-DIN_N15[tdim,zdim,jj,ji],rmin))
              rdin[jj,ji] = min(rdin[jj,ji],2*rn15std)
              rdin[jj,ji] = max(rdin[jj,ji],rn15std/2)
              #bassim = rno3 + eps_assim*(1-uno3)/uno3*log(1-uno3)*rno3/1000.
              bassimphn[jj,ji]=rdin[jj,ji]+rdin[jj,ji]*eps_assim*(1-un[jj,ji])/un[jj,ji]*np.log(1-un[jj,ji])/1000.
              fcassimphn[jj,ji]=bassimphn[jj,ji]/(1+bassimphn[jj,ji])
              #bassim = rno3 + eps_assim*(1-uno3)/uno3*log(1-uno3)*rno3/1000.
              bassimphd[jj,ji]=rdin[jj,ji]+rdin[jj,ji]*eps_assim*(1-un[jj,ji])/un[jj,ji]*np.log(1-un[jj,ji])/1000.
              fcassimphd[jj,ji]=bassimphd[jj,ji]/(1+bassimphd[jj,ji])
              #rtdin15 = biodin15/biono3
              rtdin15[jj,ji] = DIN_N15[tdim,zdim,jj,ji] /max(DIN[tdim,zdim,jj,ji],rmin)
              rtdin15[jj,ji] = min(rtdin15[jj,ji],2*rn15std/(1+rn15std))
              rtdin15[jj,ji] = max(rtdin15[jj,ji],rn15std/(1+rn15std)/2.)
              #rtphytn15 = biophytn15/biophyt
              rtphnn15[jj,ji] = PHN_N15[tdim,zdim,jj,ji] /max(PHN[tdim,zdim,jj,ji],rmin)
              rtphdn15[jj,ji] = PHD_N15[tdim,zdim,jj,ji] /max(PHD[tdim,zdim,jj,ji],rmin)
              rtphnn15[jj,ji] = min(rtphnn15[jj,ji],2.*rn15std/(1+rn15std))
              rtphnn15[jj,ji] = max(rtphnn15[jj,ji],rn15std/(1+rn15std)/2.)
              rtphdn15[jj,ji] = min(rtphdn15[jj,ji],2.*rn15std/(1+rn15std))
              rtphdn15[jj,ji] = max(rtphdn15[jj,ji],rn15std/(1+rn15std)/2.)
              #rtzoopn15 = biozoopn15/biozoop
              rtzmin15[jj,ji] = ZMI_N15[tdim,zdim,jj,ji] /max(ZMI[tdim,zdim,jj,ji],rmin)
              rtzmen15[jj,ji] = ZME_N15[tdim,zdim,jj,ji] /max(ZME[tdim,zdim,jj,ji],rmin)
              rtzmin15[jj,ji] = min(rtzmin15[jj,ji],2.*rn15std/(1+rn15std))
              rtzmin15[jj,ji] = max(rtzmin15[jj,ji],rn15std/(1+rn15std)/2.)
              rtzmen15[jj,ji] = min(rtzmen15[jj,ji],2.*rn15std/(1+rn15std))
              rtzmen15[jj,ji] = max(rtzmen15[jj,ji],rn15std/(1+rn15std)/2.)
              #rtdetrn15 = biodetrn15/biodetr
              rtdetrn15[jj,ji] = DET_N15[tdim,zdim,jj,ji] /max(DET[tdim,zdim,jj,ji],rmin)
              rtdetrn15[jj,ji] = min(rtdetrn15[jj,ji],2.*rn15std/(1+rn15std))
              rtdetrn15[jj,ji] = max(rtdetrn15[jj,ji],rn15std/(1+rn15std)/2.)
              #rzoop = biozoopn15/(biozoop-biozoopn15)
              rmi[jj,ji] = ZMI_N15[tdim,zdim,jj,ji] /max(ZMI[tdim,zdim,jj,ji]-ZMI_N15[tdim,zdim,jj,ji],rmin)
              rmi[jj,ji] = min(rmi[jj,ji], 2.*rn15std)
              rmi[jj,ji] = max(rmi[jj,ji], rn15std/2.)
              #bexcr = rzoop - eps_excr*rzoop/1000.
              bmiexcr[jj,ji] = rmi[jj,ji] - eps_excr*rmi[jj,ji]/1000.
              fcmiexcr[jj,ji] = bmiexcr[jj,ji]/(1+bmiexcr[jj,ji])
              #rzoop = biozoopn15/(biozoop-biozoopn15)
              rme[jj,ji] = ZME_N15[tdim,zdim,jj,ji] /max(ZME[tdim,zdim,jj,ji]-ZME_N15[tdim,zdim,jj,ji],rmin)
              rme[jj,ji] = min(rme[jj,ji], 2.*rn15std)
              rme[jj,ji] = max(rme[jj,ji], rn15std/2.)
              #bexcr = rzoop - eps_excr*rzoop/1000.
              bmeexcr[jj,ji] = rme[jj,ji] - eps_excr*rme[jj,ji]/1000.
              fcmeexcr[jj,ji] = bmeexcr[jj,ji]/(1+bmeexcr[jj,ji])
              
              return 
    

#%%

def compute_fluxes():
    global fprn,zphn,fprd,zphd
    global DIN,DIN_N15,PHN,PHN_N15,PHD,PHD_N15,ZMI,ZMI_N15,ZME,ZME_N15
    global PHN_N15_flux,PHD_N15_flux,ZMI_N15_flux,ZME_N15_flux,DET_N15_flux,DIN_N15_flux
    global rn15std
    global un,rdin,bassimphn,fcassimphn,bassimphd,fcassimphn,rtdin15,rtphnn15,rtphdn15
    global rtzmin15,rtzmen15,rtdetrn15,rmi,bmiexcr,fcmiexcr,rme,bmeexcr,fcmeexcr,fmiexcr,fmeexcr,freminn
    global PHN_flux,PHD_flux,ZMI_flux,ZME_flux,DET_flux,DIN_flux
    

    fn_cons_n15=np.zeros(fprn.shape)
    fn_prod_n15=np.zeros(fprn.shape)
    fn_cons=np.zeros(fprn.shape)
    fn_prod=np.zeros(fprn.shape)
    b0=1
    xfdfrac1=0.333
    xfdfrac2=1
    xbetan=0.77
    xphi=.2

    
        
    for jj in range(fprn.shape[0]):
        for ji in range(fprn.shape[1]):
               PHN_flux[jj,ji] = b0 * ( fprn[jj,ji] * zphn[jj,ji] \
                                     - fdpn[jj,ji] \
                                     - fdpn[jj,ji] \
                                     - fdpn2[jj,ji] \
                                     - fgmipn[jj,ji] \
                                     - fgmepn[jj,ji] ) 
             
                
               PHD_flux[jj,ji] = b0 * ( fprd[jj,ji] * zphd[jj,ji] \
                                      - fdpd[jj,ji] \
                                      - fdpd2[jj,ji] \
                                      - fgmepd[jj,ji])
                   
               ZMI_flux[jj,ji] = b0 * (fmigrown15[jj,ji] \
                                  - fgmezmi[jj,ji] \
                                  - fdzmi[jj,ji]\
                                  - fdzmi2[jj,ji])
                   
               ZME_flux[jj,ji] = b0 * (fmegrown15[jj,ji]            \
                                    - fdzme[jj,ji] \
                                    - fdzme2[jj,ji])

               DET_flux[jj,ji] = b0 * (fdpn[jj,ji] \
                                    +  (1.0 - xfdfrac1) * fdpd[jj,ji] \
                                    + fdzmi[jj,ji]                    \
                                    + (1.0 - xfdfrac2) *fdzme[jj,ji] \
                                    + ((1.0 - xbetan) *  (finmi[jj,ji] + finme[jj,ji])) \
                                    - fgmid[jj,ji]                   \
                                    - fgmed[jj,ji]                   \
                                    - fdd[jj,ji]                     \
                                    + fslowgain[jj,ji] - fslowloss[jj,ji]             \
                                    - (f_sbenin_n[jj,ji] / fse3t[jj,ji])           \
                                    + ffast2slown[jj,ji] )            
  
               fn_cons[jj,ji] = -(fprn[jj,ji] * zphn[jj,ji])\
                                      -(fprd[jj,ji] * zphd[jj,ji])        
                 
               fn_prod[jj,ji] =(xphi * (fgmipn[jj,ji] + fgmid[jj,ji])) \
                          + (xphi * (fgmepn[jj,ji]   \
                                  + fgmepd[jj,ji]    \
                                  + fgmezmi[jj,ji]   \
                                  + fgmed[jj,ji]))  \
                                  + fmiexcr[jj,ji]                           \
                                  + fmeexcr[jj,ji]                           \
                                  + fdd[jj,ji]              \
                                  + freminn[jj,ji]                           \
                                  + fdpn2[jj,ji]             \
                                  + fdpd2[jj,ji]             \
                                  + fdzmi2[jj,ji]            \
                                  + rtzmen15[jj,ji]*fdzme2[jj,ji]             
               DIN_flux[jj,ji] = b0 * ( fn_prod[jj,ji] + fn_cons[jj,ji] )
               
               PHN_N15_flux[jj,ji] = b0 * ( fcassimphn[jj,ji]*fprn[jj,ji] * zphn[jj,ji] \
                                     - rtphnn15[jj,ji]*fdpn[jj,ji] \
                                     - rtphnn15[jj,ji]*fdpn[jj,ji] \
                                     - rtphnn15[jj,ji]*fdpn2[jj,ji] \
                                     - rtphnn15[jj,ji]*fgmipn[jj,ji] \
                                     - rtphnn15[jj,ji]*fgmepn[jj,ji] ) 
             
                
               PHD_N15_flux[jj,ji] = b0 * ( fcassimphd[jj,ji]*fprd[jj,ji] * zphd[jj,ji] \
                                      - rtphdn15[jj,ji]*fdpd[jj,ji] \
                                      - rtphdn15[jj,ji]* fdpd2[jj,ji] \
                                      - rtphdn15[jj,ji]*fgmepd[jj,ji])
                   
               ZMI_N15_flux[jj,ji] = b0 * (fmigrown15[jj,ji] \
                                  - rtzmin15[jj,ji]*fgmezmi[jj,ji] \
                                  - rtzmin15[jj,ji]*fdzmi[jj,ji]\
                                  - rtzmin15[jj,ji]*fdzmi2[jj,ji])
                   
               ZME_N15_flux[jj,ji] = b0 * (fmegrown15[jj,ji]            \
                                    - rtzmen15[jj,ji]*fdzme[jj,ji] \
                                    - rtzmen15[jj,ji]*fdzme2[jj,ji])

               DET_N15_flux[jj,ji] = b0 * (rtphnn15[jj,ji]*fdpn[jj,ji] \
                                    +  (1.0 - xfdfrac1) *rtphdn15[jj,ji]* fdpd[jj,ji] \
                                    + rtzmin15[jj,ji]*fdzmi[jj,ji]                    \
                                    + (1.0 - xfdfrac2) * rtzmen15[jj,ji]*fdzme[jj,ji] \
                                    + ((1.0 - xbetan) *  (rtzmin15[jj,ji]*finmi[jj,ji] + rtzmen15[jj,ji]*finme[jj,ji])) \
                                    - rtdetrn15[jj,ji]*fgmid[jj,ji]                   \
                                    - rtdetrn15[jj,ji]*fgmed[jj,ji]                   \
                                    - rtdetrn15[jj,ji]*fdd[jj,ji]                     \
                                    + fslowgain[jj,ji] - fslowloss[jj,ji]             \
                                    - (f_sbenin_n[jj,ji] / fse3t[jj,ji])           \
                                    + ffast2slown[jj,ji] )            
  
               fn_cons_n15[jj,ji] = -(fcassimphn[jj,ji]*fprn[jj,ji] * zphn[jj,ji])\
                                      -(fcassimphd[jj,ji]*fprd[jj,ji] * zphd[jj,ji])        
                 
               fn_prod_n15[jj,ji] =(xphi * (rtphnn15[jj,ji]*fgmipn[jj,ji] + rtdetrn15[jj,ji]*fgmid[jj,ji])) \
                          + (xphi * (rtphnn15[jj,ji]*fgmepn[jj,ji]   \
                                  + rtphdn15[jj,ji]*fgmepd[jj,ji]    \
                                  + rtzmin15[jj,ji]*fgmezmi[jj,ji]   \
                                  + rtdetrn15[jj,ji]*fgmed[jj,ji]))  \
                                  + fmiexcr[jj,ji]                           \
                                  + fmeexcr[jj,ji]                           \
                                  + rtdetrn15[jj,ji]*fdd[jj,ji]              \
                                  + freminn[jj,ji]                           \
                                  + rtphnn15[jj,ji]*fdpn2[jj,ji]             \
                                  + rtphdn15[jj,ji]*fdpd2[jj,ji]             \
                                  + rtzmin15[jj,ji]*fdzmi2[jj,ji]            \
                                  + rtzmen15[jj,ji]*fdzme2[jj,ji]             
               DIN_N15_flux[jj,ji] = b0 * ( fn_prod_n15[jj,ji] + fn_cons_n15[jj,ji] )
    return


def compute_new_concentrations():
    dayinsec=86400.
    global DIN_next,DET_next,PHN_next,PHD_next,ZMI_next,ZME_next
    global DIN_N15_next,DET_N15_next,PHN_N15_next,PHD_N15_next,ZMI_N15_next,ZME_N15_next
    DIN_next[0,0,:,:]=DIN[0,0,:,:]+DIN_flux/dayinsec
    DET_next[0,0,:,:]=DET[0,0,:,:]+DET_flux/dayinsec
    PHN_next[0,0,:,:]=PHN[0,0,:,:]+PHN_flux/dayinsec
    PHD_next[0,0,:,:]=PHD[0,0,:,:]+PHD_flux/dayinsec
    ZMI_next[0,0,:,:]=ZMI[0,0,:,:]+ZMI_flux/dayinsec
    ZME_next[0,0,:,:]=ZME[0,0,:,:]+ZME_flux/dayinsec    
    DIN_N15_next[0,0,:,:]=DIN_N15[0,0,:,:]+DIN_N15_flux/dayinsec
    DET_N15_next[0,0,:,:]=DET_N15[0,0,:,:]+DET_N15_flux/dayinsec
    PHN_N15_next[0,0,:,:]=PHN_N15[0,0,:,:]+PHN_N15_flux/dayinsec
    PHD_N15_next[0,0,:,:]=PHD_N15[0,0,:,:]+PHD_N15_flux/dayinsec
    ZMI_N15_next[0,0,:,:]=ZMI_N15[0,0,:,:]+ZMI_N15_flux/dayinsec
    ZME_N15_next[0,0,:,:]=ZME_N15[0,0,:,:]+ZME_N15_flux/dayinsec    
    return
         


#%% define fluxes and set some random values
rn15std=0.0036765
fprn=np.zeros((DIN.shape[2],DIN.shape[3]))+.1
zphn=np.zeros(fprn.shape)+.1
zphd=np.zeros(fprn.shape)+.2
fprd=np.zeros(fprn.shape)+.3

fdpn=np.zeros(fprn.shape)+.1
fdpn2=np.zeros(fprn.shape)+.2
fgmipn=np.zeros(fprn.shape)+.3
fgmepn=np.zeros(fprn.shape)+.4


fdpd=np.zeros(fprn.shape)+.2
fdpd2=np.zeros(fprn.shape)+.3
fgmepd=np.zeros(fprn.shape)+.5

fmigrown15=np.zeros(fprn.shape)+.05
fgmezmi=np.zeros(fprn.shape)+.15
fdzmi=np.zeros(fprn.shape)+.25
fdzmi2=np.zeros(fprn.shape)+.35

fmegrown15=np.zeros(fprn.shape)+.15
fdzme=np.zeros(fprn.shape)+.25
fdzme2=np.zeros(fprn.shape)+.35

finmi=np.zeros(fprn.shape)+.1
finme=np.zeros(fprn.shape)+.2
fgmid=np.zeros(fprn.shape)+.3
fgmed=np.zeros(fprn.shape)+.2
fdd=np.zeros(fprn.shape)+.22
fslowgain=np.zeros(fprn.shape)
fslowloss=np.zeros(fprn.shape)
f_sbenin_n=np.zeros(fprn.shape)
fse3t=np.zeros(fprn.shape)+.1
ffast2slown=np.zeros(fprn.shape)

fmiexcr=np.zeros(fprn.shape)+.1
fmeexcr=np.zeros(fprn.shape)+.2
freminn=np.zeros(fprn.shape)+.2

bmiexcr=np.zeros(fprn.shape)
bmeexcr=np.zeros(fprn.shape)
#%% define ratio fields 
un=np.zeros(fprn.shape)
bassimphn=np.zeros(fprn.shape)
bassimphd=np.zeros(fprn.shape)
rdin=np.zeros(fprn.shape)
fcassimphd=np.zeros(fprn.shape)
fcassimphn=np.zeros(fprn.shape)
rtdin15=np.zeros(fprn.shape)
rtphnn15=np.zeros(fprn.shape)
rtphdn15=np.zeros(fprn.shape)
rtzmin15=np.zeros(fprn.shape)
rtzmen15=np.zeros(fprn.shape)
rtdetrn15=np.zeros(fprn.shape)
fcmiexcr=np.zeros(fprn.shape)
fcmeexcr=np.zeros(fprn.shape)
rmi=np.zeros(fprn.shape)
rme=np.zeros(fprn.shape)     
#%% define tracer fluxes
PHN_flux=np.zeros(fprn.shape)
PHD_flux=np.zeros(fprn.shape)
ZMI_flux=np.zeros(fprn.shape)
ZME_flux=np.zeros(fprn.shape)
DET_flux=np.zeros(fprn.shape)
DIN_flux=np.zeros(fprn.shape)

PHN_N15_flux=np.zeros(fprn.shape)
PHD_N15_flux=np.zeros(fprn.shape)
ZMI_N15_flux=np.zeros(fprn.shape)
ZME_N15_flux=np.zeros(fprn.shape)
DET_N15_flux=np.zeros(fprn.shape)
DIN_N15_flux=np.zeros(fprn.shape)
#%% model run
calculate_ratios()
compute_fluxes()
compute_new_concentrations()


#%%


plt.close('all')
fig = plt.figure(figsize=(15,10))
ax1 = fig.add_subplot(231)
plt.pcolormesh(np.where(tmask[0,:,:]==0,1,delta(PHN_next[0,0,:,:],PHN_N15_next[0,0,:,:])),vmin=4.9,vmax=5.1)
plt.colorbar()
plt.text(10,10,'PHN',color='white')
ax1 = fig.add_subplot(232)
plt.pcolormesh(np.where(tmask[0,:,:]==0,1,delta(PHD_next[0,0,:,:],PHD_N15_next[0,0,:,:])),vmin=4.9,vmax=5.1)
plt.colorbar()
plt.text(10,10,'PHD',color='white')
ax1 = fig.add_subplot(233)
plt.pcolormesh(np.where(tmask[0,:,:]==0,1,delta(ZMI_next[0,0,:,:],ZMI_N15_next[0,0,:,:])),vmin=4.9,vmax=5.1)
plt.colorbar()
plt.text(10,10,'ZMI',color='white')
ax1 = fig.add_subplot(234)
plt.pcolormesh(np.where(tmask[0,:,:]==0,1,delta(ZME_next[0,0,:,:],ZME_N15_next[0,0,:,:])),vmin=4.9,vmax=5.1)
plt.colorbar()
plt.text(10,10,'ZME',color='white')
ax1 = fig.add_subplot(235)
plt.pcolormesh(np.where(tmask[0,:,:]==0,1,delta(DIN_next[0,0,:,:],DIN_N15_next[0,0,:,:])),vmin=4.9,vmax=5.1)
plt.colorbar()
plt.text(10,10,'DIN',color='white')
ax1 = fig.add_subplot(236)
plt.pcolormesh(np.where(tmask[0,:,:]==0,1,delta(DET_next[0,0,:,:],DET_N15_next[0,0,:,:])),vmin=4.9,vmax=5.1)
plt.colorbar()
plt.text(10,10,'DET',color='white')

plt.savefig('Updated_surface_concentration.png',dpi=200,pading=.1,bbox_inches='tight')
