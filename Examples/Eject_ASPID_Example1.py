"""
Created on Thu Agu 17 11:00:00 2023
@author: Matias Chacón de la Cruz
"""
import os
dir1= r'D:\OneDrive\Papers\Timber\_Hysteresis_ASPID\Supp_Material\Python_codes'
os.chdir(dir1)

import numpy as NP
import matplotlib.pyplot as PL
from   matplotlib import gridspec
import matplotlib as mpl
import ASPID as MAT


#===================================
#A. Model parameters
#===================================
#-------------------------
#Positive envelope
#-------------------------
modP= 4
Ko= 4e0
Fo= 3e1
r1= 1/10
xu= 10
xf= 25
ParamP= [modP,Ko,Fo,r1,xu,xf]

#-------------------------
#Negative envelope
#-------------------------
modN= 4
Ko= 4e0
Fo= 2e1
r1= 1/10
xu= 5
xf= 35
ParamN= [modN,Ko,Fo,r1,xu,xf]


#-------------------------
#Hysteretics
#-------------------------
Po=  [2e0,-5e0]
Kpo= [1e0,4e0]
Kro= [6e0,7e0]
Ku=  [1e1,1e1]
Lr=  [0.5,0.4]
Lu=  [0.2,0.4]
Lt=  [0.2,0.3]
Ap=  [0.2,0.2]
Ar=  [0.05,0.05]
Ao=  [0.15,0.2]
Bx=  3e-4
Be=  7e-3

prop=  [Po,Kpo,Kro,Ku,
        Lr,Lu,Lt,
        Ap,Ar,Ao,Bx,Be,
        ParamP,ParamN]


#=========================
#B. Algorithm
#=========================
#-----------------------
#Displacement input
#-----------------------
D=     NP.zeros([0])
Ncyc=  50
Dpeak= NP.array([0,5,-5,10,-10,15,-15,20,-20])

for i in range(len(Dpeak)-1):
  d1= Dpeak[i]
  d2= Dpeak[i+1]
  Di= NP.linspace(d1,d2,Ncyc)
  D=  NP.hstack([D,Di[:-1]])
#
Nstep= len(D)
Disp_= D

#-------------------------
#State variables
#-------------------------
Nstate= 13
Ustatev= NP.zeros([Nstate])
Ustatev[5]= Po[0]
Ustatev[6]= Po[1]
Force= NP.zeros([1])

#-------------------------
#Cycle solution
#-------------------------
X=  NP.zeros([Nstep])
Y=  NP.zeros([Nstep])
Ed= NP.zeros([Nstep])

for i in range(Nstep-1):
  Disp=   NP.array([Disp_[i]])       #previous strain step n
  Disp_n= NP.array([Disp_[i+1]])     #actual strain step n+1
  dDisp=  Disp_n - Disp
  
  Force,Ustatev= MAT.Main(Disp,dDisp,Ustatev,prop,Force)
  X[i]=  Disp_n[0]
  Y[i]=  Force[0]
  Ed[i]= Ustatev[2]
#


#=========================
#C. Plots outputs
#=========================
#-------------------------
#Settings
#-------------------------
PL.close('all')
fig= PL.figure(1,figsize=(5*1.6,5), dpi=80)
fig.patch.set_facecolor('w')    #axes background color
gs=  gridspec.GridSpec(2,2)
mpl.rcParams['text.usetex']=True

ax1= fig.add_subplot(gs[0,0])
ax2= fig.add_subplot(gs[1,0],sharex=ax1)
ax3= fig.add_subplot(gs[0,1])
ax4= fig.add_subplot(gs[1,1])
gs.update(wspace=0.4, hspace=0.0)
for ax in [ax1]:
    PL.setp(ax.get_xticklabels(),visible=False)
#

#-------------------------
#Data
#-------------------------
U= X[:-1]  #displ.
F= Y[:-1]  #force
T= NP.arange(1,Nstep) #time

ax1.plot(U,F,'orangered',linewidth=0.7)
ax3.plot(T,Ed[:-1]/1e3,'orange')
ax4.plot(T,F,'teal')

ax1.set_xlabel('Displacement, '+r'$u$'+' [mm]')
ax1.set_ylabel('Force, '+r'$F$'+' [kN]')
ax3.set_ylabel('Dissipated energy, '+r'$\Upsilon_d$'+' [kJ]')
ax4.set_xlabel('Pseudo-time step, '+r'$t$'+' [-]')
ax4.set_ylabel('Force, '+r'$F$'+' [kN]')

PL.savefig('Fig_Example1.pdf')
