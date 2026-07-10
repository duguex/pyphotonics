# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 16:04:47 2023

@author: 87930
"""
import numpy as np
import matplotlib.pyplot as plt

import matplotlib

def hermite_poly(n,x):
    if n==0 :
        return 1.0+0.0*x
    elif n==1:
        return 2*x
    else:
        return 2*x* hermite_poly(n-1,x)-2*(n-1)*hermite_poly(n-2,x)
    
def vibration_wave_function(n,x):
    if n<10:
        return np.pi**-0.25/ np.sqrt(2**n*np.math.factorial(n)) * hermite_poly(n,x) * np.exp(-x*x/2.0)
    
    if n>=10:
        d=np.sqrt(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5)*(n-6)*(n-7)*(n-8)*(n-9))
        return np.pi**-0.25/ np.sqrt(2**n*np.math.factorial(n-10)) * hermite_poly(n,x) * np.exp(-x*x/2.0-np.log(d))
    
def vibration_wave_function2(n,x,w):
    return w**0.25*vibration_wave_function(n, w**0.5*x)

def dot(n1,n2,w1,w2):
    ratio=w2/w1
    rang=1000
    res=10000
    dx=2*rang/res
    x = np.linspace(rang,-rang,res,dtype=np.float64)
    intergal=0.0
    for xx in x:
        intergal+=vibration_wave_function(n1,xx)*vibration_wave_function2(n2,xx,ratio)*dx
        
    return intergal

def test():
    x = np.linspace(-15,15,1000,dtype=np.float64)
    fig, ax = plt.subplots(figsize=(8,6),dpi=300)
    #ax.plot(x,hermite_poly(0,x),label="0")
    ax.plot(x,vibration_wave_function(16,x),label=str(30))
    ax.legend()
    fig.show()
    
#for k in range(0,50):
    #print(k,round(dot(5 , k, 1.0,0.8),5))
#test()






