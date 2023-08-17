"""
Created on Thu Agu 17 11:00:00 2023
@author: Matias Chacón de la Cruz
"""
from math import exp,acos,atan,sqrt
from numpy import sign


#===============================================
# Heaviside
#===============================================
def Heaviside(x,h0,Out):
    #-------------------
    #Heaviside function
    #-------------------
    Zero= 0.0
    One=  1.0

    #-------------------
    #Exact
    #-------------------
    if x>Zero:
        H=  One
    elif x<Zero:
        H=  Zero
    else:
        H=  h0
    #
    dH= Zero

    #-------------------
    #output
    #-------------------
    if (Out==0):
        return H
    else:
        return dH
    #
#


#==================================
# Cubic Polynomial, C-1 class
#==================================
def Poly3(x,x1,f1,df1,x2,f2,df2,opt):
    #------------------------
    #polinomio 3 orden genera continuidad C1 en curvas
    #------------------------
    #x1, f1, df1: valor x1, f(x1) y df(x1) inicio
    #x1, f1, df1: valor x2, f(x2) y df(x2) final
    #------------------------
    #Constantes
    #------------------------
    Two=    2.0
    Three=  3.0
    sqTiny= 1e-15
    dx=     x2 - x1
    
    #------------------------
    #Near-zero increment dx
    #------------------------
    if abs(dx)>sqTiny:
        #------------------------
        #constants
        #------------------------
        c1= ((Two*df1 + df2)*dx - Three*(f2-f1))/dx**2
        c2= ((    df1 + df2)*dx -   Two*(f2-f1))*x1/dx**3
        
        a=  ((df1 + df2)*dx - Two*(f2 - f1))/dx**3
        b=  -(c1 + Three*c2)
        c=  df1 + x1*(Two*c1 + Three*c2)
        d=  f1  - x1*(df1 + x1*(c1 + c2))
        
        #------------------------
        #polinomio 3 orden y su derivada
        #------------------------
        f=        a*x**3 +     b*x**2 + c*x + d
        df= Three*a*x**2 + Two*b*x    + c
    #------------------------
    #Near-zero increment dx
    #------------------------
    else:
        f=  f1
        df= df1
    #
    #------------------------
    #output
    #------------------------
    if opt==0:
        return f
    else:
        return df
    #
#


#==================================
#Cuadratic Bezier curve C-1 class
#==================================
def Bezier3(x,x0,f0,df0,x2,f2,df2,opt):
    #--------------------
    #Constants
    #--------------------
    One=  1.0
    Two=  2.0
    Half= 0.5
    
    #--------------------
    #Slope
    #--------------------
    dx=  x2 - x0
    dy=  f2 - f0
    idx= Ivar(dx)
    m=   dy*idx
    sm=  sign(m)
    ss=  sign(df0*df2)
    
    #--------------------
    #Zone and shape case
    #--------------------
    zone= 0
    if x>=x0 and x<=x2:
        zone= 1
    #
    
    #--------------------
    #Adjust slope if necesary
    #--------------------
    E0= abs(df0)
    E2= abs(df2)
    s0= sign(df0 - m)
    s2= sign(df2 - m)
    flag= sm*s2
    
    if s0*s2==1:
        E0=   m
        E2=   m
        flag= 0
    #
    
    #--------------------
    #P1
    #--------------------
    if flag==1:
        dx1= (ss*E0*dx-dy)/(ss*E0-E2)
        dy1= E0*(E2*dx-dy)/(E2-ss*E0)
        x1=  x2 - dx1
        f1=  f0 + ss*sm*dy1
    elif flag==-1:
        dx1= (ss*E2*dx-dy)/(ss*E2-E0)
        dy1= E2*(E0*dx-dy)/(E0-ss*E2)
        x1=  x0 + dx1
        f1=  f2 - ss*sm*dy1
    #
    
    #--------------------
    #f
    #--------------------
    if abs(flag)*zone>0:
        d1=  x0 + x2 - Two*x1
        d2=  x1**2 - x0*x2 + d1*x
        id1= Ivar(d1)
        id2= Ivar(d2)
        t=   (x0 - x1 + sqrt(d2))*id1
        f=   f1 + (f0-f1)*(One-t)**2 + (f2-f1)*t**2
    else:
        t= (x-x0)*idx
        f= f0 + t*(f2-f0)
    #
    
    #--------------------
    #df
    #--------------------
    if abs(flag)*zone>0:
        dyt= Two*(One-t)*(f1-f0) + Two*t*(f2-f1)
        dtx= Half*sqrt(abs(id2))
        df=  dyt*dtx
    else:
        df= (f2-f0)*idx
    #
    
    if opt==0:
        return f
    else:
        return df
    #
#


#==================================
#Inverse of variable
#==================================
def Ivar(x):
    #--------------------
    #Constants
    #--------------------
    Zero= 0.0
    One=  1.0

    ix= Zero
    if x!=Zero:
        ix= One/x
    #
    return ix
#


#==================================
#Find if variable fall in a stated range
#==================================
def Inrange(x,x1,x2):
    out= 0
    if x1<x2:
        if x>=x1 and x<=x2:
            out= 1
        #
    else:
        if x>=x2 and x<=x1:
            out= 1
        #
    #
    return out
#