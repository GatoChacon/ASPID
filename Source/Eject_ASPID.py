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
modP= 3
Ko= 2.6e0
Fo= 7e0
r1= 1/25
r2= -1/70
xu= 25
xf= 45
ll= -0.02
ParamP= [modP,Ko,Fo,r1,r2,xu,xf,ll]

modP= 1
Xp= NP.array([[0,0],[3,4],[10,5],[20,4.5],[40,1]])
ParamP= [modP,Xp]

#-------------------------
#Negative envelope
#-------------------------
modN= 3
Ko= 3.6e0
Fo= 3e0
r1= 1/25
r2= -1/20
xu= 27
xf= 50
ll= -0.02
ParamN= [modN,Ko,Fo,r1,r2,xu,xf,ll]


#-------------------------
#Hysteretics
#-------------------------
Po=  [1.2,-2]
Kpo= [1e-1,5e-1]
Kro= [6e0,6e0]
Ku=  [7e0,7e0]
Lr=  [0.45,0.65]
Lu=  [0.075,0.11]
Lt=  [0.2,0.2]
Ap=  [0.065,0.065]
Ar=  [0.075,0.075]
Ao=  [-0.01,0.03]
Bx=  1e-4
Be=  4e-4

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
Ncyc=  20
D=     NP.zeros([0])
Dpeak= NP.array([0,0,60])
Dpeak= NP.array([0,0,5,-5,10,-10,20,-20,30,-30,40,-40,40,-40,40,-40,40])
# Dpeak= NP.array([0,0,10,-3,20,-3,30,-3,25,-3,40,-3,40,-3,40,-3,40,-3])

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