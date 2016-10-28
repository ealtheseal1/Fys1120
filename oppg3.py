# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 21:59:06 2016

@author: Eirik
"""

import numpy as np
import matplotlib.pyplot as plt
import pylab as p
import mpl_toolkits.mplot3d.axes3d as p3

e = 1.6e-19
m_p = 1.672621637e-27 #Reference Pearson
r_D = 50.0E-03# radius
D = 90.0E-06 #valleygap
c = 3.0E08 #speed of light
start = 0
stop = 300.0e-09
dt = 100.0e-15
t = np.linspace(start,stop,(stop-start)/dt)
t_test = np.linspace(start,100000,10000000)
n = len(t)

E_o = (25.0E03/90.0e-6) #V/m
B = np.array([0.0,0.0,2.0])
w = (np.abs(e)/m_p)*np.linalg.norm(B)
f = w/(2*np.pi)

def oppgave3():
    
    A = np.zeros((n,3))
    V = np.zeros((n,3))
    R = np.zeros((n,3))
    E = np.zeros((n,3))    
    tid = np.zeros(n)    
    
    V[0] = np.array([0.0 , 0.0 , 0.0])
    R[0] = np.array([0.0 , 0.0 , 0.0]) 
    tid[0] = 0
    
     #omega
    f = w/(2*np.pi) #cyclotron frekvens
    

    for i in xrange(n-1):
    
        if R[i,0] > -D/2 and R[i,0]< D/2:
            E = np.array([E_o*np.cos(w*tid[i]) , 0, 0])
                    
        else:
            E = np.zeros(3)
            
               
        if not np.linalg.norm(R[i,:]) > R_D:
            
            A = (e/m_p)*(E + np.cross(V[i],B))
        
        else:
            A = 0
        
        
        V[i+1] = V[i] + A*dt
        R[i+1] = R[i] +V[i+1]*dt
        tid[i+1] = tid[i] + dt
    
    dropout = np.linalg.norm(R[-1,:])
    print "Escape velocity is: %g m/s" %(np.linalg.norm(V[-1,:]))
    print "Percentage achieved of speed of light: %g" %((np.linalg.norm(V[-1,:])/c) *100) 
    return R,V,tid, dropout,E,A
"""
def sammenligning():
    R_D =(4.5 - 2.1) #extract minus inject
    Breal = np.array([0.0 , 0.0 , 8*1.1]) #8 magnets with 1.1 T each
    harmonic_num = 6
    freq = 50.0e06 #Hz
    omega = freq*2*np.pi
    
    rel_gam = (1.0/npsqrt(1-(vel/c)**2))
    rel_m = rel_gam*m_p
    rel_w =rel_gam*omega
    rel_r = vel/rel_w 
    
    for i in xrange(np.len(t_test)-1):
    
        if R2[i,0] > -D/2 and R[i,0]< D/2:
            E = np.array([E_o*np.cos(w*tid[i]) , 0, 0])
                    
        else:
            E = np.zeros(3)
            
               
        if not np.linalg.norm(R[i,:]) > r_D:
            
            A = (e/m_p)*(E + np.cross(V[i],B))
        
        else:
            A = 0
"""

if __name__ == "__main__":
    R,V,tid,dropout,E,A = oppgave3()
    
    fig1 = plt.figure("Oppgave 3", figsize=(9,9))
    ax1  = fig1.add_subplot(1,1,1)
    #ax1.scatter(R[:,0], R[:,1], color = 'blue')
    ax1.set_xlabel("X Pos[m]"), ax1.set_ylabel("Y Pos[m]")
    p.title("Oppg 3 - $\delta t = 100fS$")  
    p.plot(R[:,0], R[:,1])
    p.show()
    
    fig2 = plt.figure("Oppgave 3", figsize=(7,7))
    plt.plot(tid,R[:,0],label="$x(t)$")  
    plt.plot(tid,R[:,1],label="$y(t)$")
    plt.plot(tid,R[:,2],label="$z(t)$")
    plt.title("Oppg3 - $\delta t = 100fS$")
    plt.xlabel("Tid $nS$")
    plt.ylabel("Posisjon [m]")
    plt.legend()
    plt.show()    
    
    
    fig3 = plt.figure("Oppgave 3", figsize=(7,7))
    plt.plot(tid, V[:,0], "y", label="$V_{x}(t)$")
    plt.plot(tid, V[:,1], "r", label = "$V_{y}(t)$")
    plt.plot(tid, V[:,2], "b", label = "$V_{z}$(t)")
    plt.title("Oppg 3 - $\delta t = 100fS$")
    plt.xlabel("Tid $nS$")
    plt.ylabel("V, $m/S$")
    plt.legend()
    plt.show()    

    
              
            