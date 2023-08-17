"""
Created on Thu Agu 17 11:00:00 2023
@author: Matias Chacón de la Cruz
"""

import numpy as NP
from   math import exp
from   numpy import sign
from   FunctionsR1 import Heaviside
from   FunctionsPaths import Fmono, Fpinch, Ftransit


def Main(Disp,dDisp,Ustatev,Prop,Force):
    #================================
    #Initial settings
    #================================
    #----------------------
    # Auxiliar variables
    #----------------------
    Zero= 0.0
    One=  1.0
    Half= 0.5
    Tol1= 1.0e-10
    Tol2= 1.0e10
    
    #----------------------
    # State variables
    #----------------------
    sens_t= int(Ustatev[0])  #flag of last non-zero direction
    flag_t= int(Ustatev[1])  #flag of state curve
    Ed_t=   Ustatev[2]       #dissipated energy
    rp_t=   Ustatev[3]       #maximum displacement +
    rn_t=   Ustatev[4]       #minimum displacement -
    fp_t=   Ustatev[5]       #final force +
    fn_t=   Ustatev[6]       #final force -
    xr_t=   Ustatev[7]       #reversal displacement
    fr_t=   Ustatev[8]       #reversal force
    ri_t=   Ustatev[9]       #initial of transit curve disp
    fi_t=   Ustatev[10]      #initial of transit curve force
    wp_t=   Ustatev[11]      #damage +
    wn_t=   Ustatev[12]      #damage -
    
    #----------------------
    # Parameters
    #----------------------
    Po=  Prop[0]
    Kpo= Prop[1]
    Kro= Prop[2]
    Ku=  Prop[3]
    Lr=  Prop[4]
    Lu=  Prop[5]
    Lt=  Prop[6]
    Ap=  Prop[7]
    Ar=  Prop[8]
    Ao=  Prop[9]
    Bx=  Prop[10]
    Be=  Prop[11]
    ParamP= Prop[12]
    ParamN= Prop[13]
    
    
    #================================
    #Fields
    #================================
    #----------------------
    # Disp, increment, force
    #----------------------
    x_t= Disp[0]
    f_t= Force[0]
    dx=  dDisp[0]
    x=   x_t + dx
    
    #---------------------
    #Non-zero disp. increment
    #---------------------
    hdx=  Heaviside(abs(dx)-Tol1,Zero,0)
    flow= int(hdx)*sign(dx)
    sens= sens_t
    if hdx==1:
        sens= flow
    #
    
    #---------------------
    #r+-, rmax
    #---------------------
    hp= Heaviside(x,Zero,0)
    hn= Heaviside(-x,Zero,0)
    rp= max(rp_t,x*hp)
    rn= min(rn_t,x*hn)
    rm= max(rp,abs(rn),Tol1)
    
    #---------------------
    #Initial values
    #---------------------
    fp= fp_t
    fn= fn_t
    wp= wp_t
    wn= wn_t
    Xr= [xr_t,fr_t]
    Xi= [ri_t,fi_t]
    
    #================================
    #Damaged zone and load case
    #================================
    #---------------------
    #Zone: 0:no damaged, 1: damaged
    #---------------------
    zone= 0
    if rp_t==0 and rn_t==0:
        zone= 0
    elif x>=rn_t and x<=rp_t:
        zone= 1
    #
    tagp= int(Heaviside(rp-rp_t,Zero,0))
    tagn= int(Heaviside(-rn+rn_t,Zero,0))

    #---------------------
    #Load case= +-1: monotonic, +-2:pinching, +-3:transition
    #---------------------
    flag= flow*(zone + 1)
    if zone==0 and abs(flag_t)>=2:
        flag= 3*flow
    #

    #================================
    #Updated force
    #================================
    if flow==0:
        f= f_t
    else:
        #---------------------
        #Damaged strength at point F
        #---------------------
        fp= (One - wp_t)*fp_t
        fn= (One - wn_t)*fn_t
        
        #---------------------
        #Point P+-
        #---------------------
        Pop= Po[0]
        Pon= Po[1]
        dP=  Pop - Pon
        Plp= Pop + sign(Ao[0])*dP*Half
        Pln= Pon + sign(Ao[1])*dP*Half
        Pp=  Plp + (Pop-Plp)*exp(-abs(Ao[0])*rm)
        Pn=  Pln + (Pon-Pln)*exp(-abs(Ao[1])*rm)
        P=   [Pp,Pn]
        
        #---------------------
        #Damaged pinching stiffness
        #---------------------
        Kpp= Kpo[0]*exp(-Ap[0]*rm)
        Kpn= Kpo[1]*exp(-Ap[1]*rm)
        Kp=  [Kpp,Kpn]

        #---------------------
        #Damaged re-loading stiffness
        #---------------------
        Krp= Kro[0]*exp(-Ar[0]*rm)
        Krn= Kro[1]*exp(-Ar[1]*rm)
        Kr=  [Krp,Krn]
        
        Cons= [P,Kp,Kr,Ku,Lr,Lu,Lt]
        
        #---------------------
        #Reversal point
        #---------------------
        if sens_t*flow<0:
            Xr= [x_t,f_t]
        #
        
        #---------------------
        #Final point
        #---------------------
        Xf= [rp,fp]
        if flow==-1:
            Xf= [rn,fn]
        #
        
        #---------------------
        #Initial point of transition
        #---------------------
        if abs(flag_t)==2:
            Xi= [rp_t,(One-wp_t)*fp_t]
            if flow==-1:
                Xi= [rn_t,(One-wn_t)*fn_t]
            #
        #
        
        #---------------------
        #U/R curve
        #---------------------
        if zone==1:
            f= Fpinch(x,Xr,Xf,Cons,flow)
        
        #---------------------
        #Envelopes/transition
        #---------------------
        else:
            #---------------------
            #Transition
            #---------------------
            if abs(flag)==3:
                f=  Ftransit(x,Xr,Xi,ParamP,ParamN,Cons,flow)
                if f==Tol2:
                    flag= flow
                #
            #
            #---------------------
            #Envelopes
            #---------------------
            if abs(flag)==1:
                fmp= Fmono(x,ParamP,0)
                fmn= Fmono(-x,ParamN,0)
                f=   hp*fmp - hn*fmn
            #
        #
        #---------------------
        #Effective Final force
        #---------------------
        if tagp==1:
            fp= f
        elif tagn==1:
            fn= f
        #
    #
    #---------------------
    #Updated disp, force
    #---------------------
    Disp[0]= x
    Force[0]= f
    
    
    #================================
    # Dissipated energy, damage variables
    #================================
    Ed= Ed_t + Half*(f_t + f)*dx
    w=  (One - exp(-Bx*rm))*(One - exp(-Be*Ed))*abs(dx)
    wp= max(w,wp_t)
    wn= max(w,wn_t)
    
    
    
    #================================
    # update state variables
    #================================
    Ustatev[0]=  sens  #flag of non-zero direction
    Ustatev[1]=  flag  #flag of curve
    Ustatev[2]=  Ed    #dissipated energy
    Ustatev[3]=  rp    #maximum displacement +
    Ustatev[4]=  rn    #minimum displacement -
    Ustatev[5]=  fp    #maximum effective force +
    Ustatev[6]=  fn    #minimum effective force -
    Ustatev[7]=  Xr[0] #reversal displacement
    Ustatev[8]=  Xr[1] #reversal force
    Ustatev[9]=  Xi[0] #initial of transit curve disp
    Ustatev[10]= Xi[1] #initial of transit curve force
    Ustatev[11]= wp    #damage +
    Ustatev[12]= wn    #damage -
    
    return Force, Ustatev
#
