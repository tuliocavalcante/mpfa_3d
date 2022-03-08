# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 15:26:44 2019

@author: Túlio

Essa rotina é a que calcula os pesos associados ao MPFA3D, segundo o 
desenvolvimento matemático de Túlio Cavalcante.

"""
import numpy as np
import expweightsfunc as ewf
import weightsfunc as wf


def weights(M,node,t='semiexp'):
    
    w = np.empty(len(M.nodes.adjacences)).tolist()
    a = np.empty(len(M.nodes.adjacences)).tolist()
    
    velem = M.nodes[i].elemadjacences
    
    for j in range(len(velem)):
        
        # O nó Q deve ser interpolado. 
        # Olhando para elemento k, que compartilha Q, devo fazer o seguinte:
        # 1 - Determinar os nós I, J e K.
        # 2 - Determinar as faces de k e que são opostos a I (f1), J (f2) e 
        # K(f3) com seus valores de flag de Neumann, se for o caso.
        # 3 - Determinar os elementos que compartilham face com k e que são
        # opostos a I (o1), J (o2) e K(o3).
                
        T1 = 0.5*(Q+I); T2 = (1/3)*(Q+K+I); T3 = 0.5*(Q+K); T4 = (1/3)*(Q+K+J);
        T5 = 0.5*(Q+J); T6 = (1/3)*(Q+J+I);
        t1 = [T1,T6,T2,Q]; t2 = [T2,T4,T3,Q]; t3 = [T4,T6,T5,Q]; t4 = [T2,T6,T4,Q]

        csi = wf.calccsi( t1,t2,t3,t4,Kk )
        
        t1 = [T1,T2,k,Q]; t2 = [T2,T3,k,Q]; t3 = [T3,T4,k,Q]; t4 = [T4,T5,k,Q]
        t5 = [T5,T6,k,Q]; t6 = [T6,T1,k,Q]
              
        [ neta, sigma, nbf ] = wf.calcnetasigma( t1,t2,t3,t4,t5,t6,Kk,Ko1,Ko2, \
                                            Ko3,g1,g2,g3,o1,o2,o3,centelem,k )
        
        if t == 'semiexp':
            g = np.dot(csi,inv(neta))
            c = np.dot(g,sigma)
            d = sum(nbf) + np.dot(g,nbf)
        elif t == 'compexp':
            [ b ] = ewf.binvneta( neta ) # Inverte a matriz neta na mão!!!
            [ c, d ] = ewf.calcced( csi,sigma,b,nbf )
        
        a = a + d;
        w(j) = w(j) + c(1)
        if o1!=0:
            tp=find(velem==o1)
            w(tp) = w(tp) + c(2) 
        if o2!=0:
            tp=find(velem==o2)
            w(tp) = w(tp) + c(3) 
        if o3!=0:
            tp=find(velem==o3)
            w(tp) = w(tp) + c(4)


    if sum(w)!=0:
        a = a/sum(w)
        w = w/sum(w)
    
    return (w, a)