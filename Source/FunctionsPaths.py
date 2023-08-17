"""
Created on Thu Agu 17 11:00:00 2023
@author: Matias Chacón de la Cruz
"""
import numpy as NP
from   math import exp
from   FunctionsR1 import Poly3, Bezier3, Inrange
from scipy.interpolate import interp1d
# from scipy.interpolate import CubicSpline
from scipy.interpolate import PchipInterpolator



#===================================
#A. Envelope curves
#===================================
def Fmono(x,Param,opt):
    #-------------------
    #Constants
    #-------------------
    Zero= 0.0
    One=  1.0
    
    #-------------------
    #Pameters
    #-------------------
    mod=  int(Param[0])
    
    #-------------------
    #0. Linear
    #-------------------
    if mod==0:
        Ko= Param[1]
        
        z=  abs(x)
        f=  Ko*x
        df= Ko
    
    #-------------------
    #I. Linear-interpolation, n-points
    #-------------------
    elif mod==1:
        XY= Param[1]
        n=  len(XY)
        X=  NP.zeros(n)
        Y=  NP.zeros(n)
        for i in range(n):
            X[i]= XY[i][0]
            Y[i]= XY[i][1]
        #
        zend= X[-1]
        fend= Y[-1]
        Inter1= interp1d(X,Y)
        
        z= abs(x)
        if z<=zend:
            f=  Inter1(z)
            df= Zero           #in progress
        else:
            f=  fend
            df= Zero
        #
    
    #-------------------
    #II. Cubic-spline, n-points
    #-------------------
    elif mod==2:
        XY= Param[1]
        n=  len(XY)
        X=  NP.zeros(n)
        Y=  NP.zeros(n)
        for i in range(n):
            X[i]= XY[i][0]
            Y[i]= XY[i][1]
        #
        zend= X[-1]
        fend= Y[-1]
        #Inter3=  CubicSpline(X,Y)       #Suffer Rungue's problem
        Inter3m= PchipInterpolator(X,Y)  #ok
        
        z= abs(x)
        if z<=zend:
            f=  Inter3m(z,0)
            df= Inter3m(z,1)
        else:
            f=  fend
            df= Zero
        #
    
    #-------------------
    #III. EPHM model
    #-------------------
    elif mod==3:
        Ko= Param[1]
        Fo= Param[2]
        r1= Param[3]
        r2= Param[4]
        xu= Param[5]
        xf= Param[6]
        ll= Param[7]
        xo= Fo/Ko
        y=  x/xo
        Fu= (Fo + r1*Ko*xu)*(One - exp(-xu/xo))
        Ff= Fu + r2*Ko*(xf-xu)
        Fr= Zero
        
        if x>=Zero and x<=xu:
            f=  (Fo + r1*Ko*x)*(One - exp(-y))
            df= r1*Ko*(1-exp(-y)) + (Fo + r1*Ko*x)*exp(-y)/xo
        elif x>xu and x<=xf:
            f=  Fu+ r2*Ko*(x - xu)
            df= r2*Ko
        else:
            f=  (Ff-Fr)*exp(ll*(x-xf))
            df= (Ff-Fr)*ll*exp(ll*(x-xf))
        #
    #-------------------
    #IV. Proposed by authors
    #-------------------
    elif mod==4:
        Ko= Param[1]
        Fo= Param[2]
        r1= Param[3]
        xu= Param[4]
        xf= Param[5]
        xo= Fo/Ko
        y=  x/xo
        Fu= (Fo + r1*Ko*xu)*(One - exp(-xu/xo))
        Fr= Zero
        
        if x>=Zero and x<=xu:
            f=  (Fo + r1*Ko*x)*(One - exp(-y))
            df= r1*Ko*(1-exp(-y)) + (Fo + r1*Ko*x)*exp(-y)/xo
        elif x>xu and x<xf:
            yu= xu/xo
            dfu=r1*Ko*(1-exp(-yu)) + (Fo + r1*Ko*xu)*exp(-yu)/xo
            x1= xu;  f1= Fu;  df1= dfu
            x2= xf;  f2= Fr;  df2= 0
            f=  Poly3(x,x1,f1,df1,x2,f2,df2,0)
            df= Poly3(x,x1,f1,df1,x2,f2,df2,1)
        else:
            f= 0
            df=0
        #
    #

    if opt==0:
        return f
    else:
        return df
    #
#



#===================================
#B. Loading/unloading curves
#===================================
def Fpinch(x,Xr,Xf,Cons,flow):
    Tol1= 1.0e-15
    
    xr= Xr[0];  yr= Xr[1]
    xf= Xf[0];  yf= Xf[1]
    P=  Cons[0]
    Kp= Cons[1]
    Kr= Cons[2]
    Ku= Cons[3]
    Lr= Cons[4]
    Lu= Cons[5]
    Pp=  P[0];     Pn=  P[1]
    Kpp= Kp[0];    Kpn= Kp[1]
    Krp= Kr[0];    Krn= Kr[1]
    Kup= Ku[0];    Kun= Ku[1]
    Lrp= Lr[0];    Lrn= Lr[1]
    Lup= Lu[0];    Lun= Lu[1]
    drp= max(Tol1,Lrp*abs(xf-xr))
    drn= max(Tol1,Lrn*abs(xf-xr))
    dup= max(Tol1,Lup*abs(xf-xr))
    dun= max(Tol1,Lun*abs(xf-xr))
    
    if flow==1:
        xj= xf - drp
        xk= xr + dup
        yj= Pp + Kpp*xj
        yk= Pp + Kpp*xk
        
        #---------------
        if xj>=xk:
            if x>=xr and x<=xk:
                x1= xr;  f1= yr;  df1= Kup
                x2= xk;  f2= yk;  df2= Kpp
                f= Bezier3(x,x1,f1,df1,x2,f2,df2,0)
            elif x>xk and x<=xj:
                f= Pp + Kpp*x
            else:
                x1= xj;  f1= yj;  df1= Kpp
                x2= xf;  f2= yf;  df2= Krp
                f= Bezier3(x,x1,f1,df1,x2,f2,df2,0)
            #
        else:
            x1= xr;  f1= yr;  df1= Kup
            x2= xf;  f2= yf;  df2= Krp
            f= Bezier3(x,x1,f1,df1,x2,f2,df2,0)
        #
    else:
        xj= xf + drn
        xk= xr - dun
        yj= Pn + Kpn*xj
        yk= Pn + Kpn*xk
        
        #---------------
        if xj<=xk:
            if x>=xf and x<=xj:
                x1= xf;  f1= yf;  df1= Krn
                x2= xj;  f2= yj;  df2= Kpn
                f= Bezier3(x,x1,f1,df1,x2,f2,df2,0)
            elif x>xj and x<=xk:
                f= Pn + Kpn*x
            else:
                x1= xk;  f1= yk;  df1= Kpn
                x2= xr;  f2= yr;  df2= Kun
                f= Bezier3(x,x1,f1,df1,x2,f2,df2,0)
            #
        else:
             x1= xf;  f1= yf;  df1= Krn
             x2= xr;  f2= yr;  df2= Kun
             f= Bezier3(x,x1,f1,df1,x2,f2,df2,0)
        #
    #
    return f
#



#===================================
#C. Transition curves
#===================================
def Ftransit(x,Xr,Xf,ParamP,ParamN,Cons,flow):
    Tol1= 1e-5
    Tol2= 1e10
    
    r_t= Xf[0]
    dr=  abs(Xf[0]-Xr[0])
    Kr=  Cons[2]
    Lt=  Cons[6]
    drp= max(Tol1,Lt[0]*dr)
    drn= max(Tol1,Lt[1]*dr)
    Krp= Kr[0]
    Krn= Kr[1]
    
    if flow==1:
        x1=  r_t
        x2=  r_t + drp
        f1=  Fpinch(x1,Xr,Xf,Cons,flow)
        f2=  Fmono(x2,ParamP,0)
        df1= Krp
        df2= Fmono(x2,ParamP,1)
    else:
        x1=  r_t - drn
        x2=  r_t
        f1=  -Fmono(-x1,ParamN,0)
        f2=  Fpinch(x2,Xr,Xf,Cons,flow)
        df1= Fmono(-x1,ParamN,1)
        df2= Krn
    #
    F=   Bezier3(x,x1,f1,df1,x2,f2,df2,0)
    out= Inrange(x,x1,x2)
    
    if out==0:
        F= Tol2
    #
    
    return F
#