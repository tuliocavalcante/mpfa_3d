# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 22:01:01 2019

@author: Túlio

Aqui vão algumas funções úteis para a forma explícita do cálculo dos pesos. 

"""
import numpy as np

def binvneta(self,neta):
    #
    r = [0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5]
    b = np.zeros((6,6))
    p = neta[0,1]*neta[1,2]*neta[2,3]*neta[3,4]*neta[4,5]*neta[5,0]
    P = neta[0,0]*neta[1,1]*neta[2,2]*neta[3,3]*neta[4,4]*neta[5,5]
    p1=1
    for i in range(5):
        if neta(r(i),r(i+1))!=0:
            p1=p1*neta(r(i),r(i+1))

    c = p1/(p-P)
    
    for i in range(5):
        for j in range(5):
            if i>j:
                N = i-j-1
            elif i==j:
                N = 5
            else:
                N = 5-j+i

            num=1; den=1
            for n in range(j,j+N-1):
                num = num*(-neta(r(n),r(n)))
                
            for n in range(j-1,j+N-1):
                if neta(r(n),r(n+1))!=0:
                    den = den*neta(r(n),r(n+1))
             
            b(i,j) = c*(num/den);
    
    return b

def calcced(self,csi,sigma,b,nbf):
    #
    d=0
    for j in range(5):
        s=0
        for i in range(5):
            s = s + csi[i]*b[i,j]

        d = d + sigma[j,0]*s

    c[0] = d
    
    d = 0; e = 0
    
    t = 1
    for u in [0,2,4]:
        for j in range(u,u+1):
            s = 0
            for i in range(5):
                s = s + csi[i]*b[i,j]
    
            d = d + sigma[j,1]*s
            e = e + nbf[j]*(1+s);
    
        c(t) = d;
        f(t-1) = e;
        t = t + 1
    
    d = sum(f); 
    
    return (c, d)
