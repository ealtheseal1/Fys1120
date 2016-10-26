    # -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 23:03:21 2016

@author: Eirik
"""
import numpy as np
import matplotlib.pyplot as plt
import pylab as p
import mpl_toolkits.mplot3d.axes3d as p3

e = 1.6e-19
m_e = 9.11e-31

E = np.array([-5.0, 0.0, 0.0])#oppgave 1
E2 = np.array([-1.0, -2.0 ,5.0])#oppgave 1 plot
B = np.array([0.0, 0.0, 2.0])

start = 0
stop = 1.0e-06
stop2 = 30e-12

dt1 = 1.0e-09
dt2 = 100.0e-09
dt3 = 1.0e-15

t = np.linspace(start,stop,(stop-start)/dt1)
t2 = np.linspace(start,stop2,(stop2-start)/dt3)

n = len(t)
n2 = len(t2)

def oppgave1():
    
    R = np.zeros((n,3))
    V = np.zeros((n,3))
    R[0] = np.array([0.0, 0.0, 0.0]) 
    V[0] = np.array([0.0, 0.0, 0.0])
    tid = np.zeros(n)
    tid[0]=0

    for i in xrange(n-1):
        a=((e*E2)/m_e)
        V[i+1] = V[i] + a*dt1
        R[i+1] = R[i] + V[i+1]*dt1
        tid[i+1] = tid[i]+dt1

    return R,V,tid


def oppgave2():
    r_m = np.zeros((n2,3))
    v_m = np.zeros((n2,3))
    t_m = np.zeros(n2)
    r_m[0] = np.array([0.0,0.0,0.0])
    v_m[0] = np.array([5.0e3,0.0,2.0e03])    
    t_m[0] = 0
    for i in xrange(n2-1):
        A = (e/m_e)*(np.cross(v_m[i],B))
        v_m[i+1] = v_m[i] + A*dt3
        r_m[i+1] = r_m[i] + v_m[i+1]*dt3
        t_m[i+1] = t_m[i]+dt3
    return r_m, v_m, t_m

        
def analytisk(tid,t_m):
    a=((e*E2)/m_e)
    A = (e/m_e)*(np.cross(v_m[0],B)) 
    
    losnenX = 0.5*a[0]*tid**2
    losnenY = 0.5*a[1]*tid**2
    losnenZ = 0.5*a[2]*tid**2
    #losnto
    return  losnenX,losnenY,losnenZ, #losnto
    


def solveequation(e,m_e,r_m,v_m,t_m,dt3):
    Q = np.linalg.norm(e)
    BB = np.linalg.norm(B)
    from sympy import Symbol, solve
    T = Symbol("T")
    Eq = (T-((2*np.pi*m_e)/(Q*BB)))
    EQ = solve(Eq, T)
    
    angleone = np.arctan2(v_m[0,1],v_m[0,0])
    angletwo = np.arctan2(-v_m[1:,1],v_m[1:,0]) 
    en = plt.plot(t_m,r_m[:,0])
    to = plt.plot(t_m,np.zeros(len(t_m)))
    plt.show()
    k = np.argmin(abs(angleone-angletwo))*dt3
    print "Numeriske perioden er: %g" %k
    print "Analytisk: %g" %(np.abs(EQ))
    print "Differansen mellom Analytisk og Numerisk er:%g" %(np.abs(EQ)-k) 
    


    
if __name__ == "__main__":
    R,V,tid = oppgave1()
    r_m, v_m, t_m = oppgave2()
    losnenX,losnenY,losnenZ = analytisk(tid,t_m)
    
    fig1 = plt.figure("Oppgave 1", figsize=(7,7)) 
    plt.plot(tid, R[:,0], "y", label="x-retning")
    plt.plot(tid, R[:,1], "r", label = "y-retning")
    plt.plot(tid, R[:,2], "g",label = "z-retning")
    plt.plot(tid,losnenX, label = "Analytisk-X")
    plt.plot(tid,losnenY, label = "Analytisk-Y")
    plt.plot(tid,losnenZ, label = "Analytisk-Z")    
    plt.title("Oppg 1 - $\delta t = 1nS$")
    plt.xlabel("Tid $\mu S$")
    plt.ylabel("Posisjon, m")
    plt.legend()
    plt.show()
    
        
    fig2 =p.figure(figsize=plt.figaspect(0.5)*1.5)
    ax = fig2.gca(projection='3d')
    ax.plot(R[:,0],R[:,1],R[:,2])    
    p.title("Oppg 1 3D: $\delta t = 1nS$")
    p.show()
                        
    fig3 = plt.figure("Oppgave 2", figsize=(7,7)) 
    plt.plot(t_m, r_m[:,0], "y", label="x-retning")
    plt.plot(t_m, r_m[:,1], "r", label = "y-retning")
    plt.plot(t_m, r_m[:,2], "g",label = "z-retning")
    plt.title("Oppg 2 - $\delta t = 1fS, \vec{V_{1}}$")
    plt.xlabel("Tid $pS$")
    plt.ylabel("Posisjon, m")
    plt.legend()
    plt.show()
    
    fig4 = plt.figure("Oppgave 2", figsize=(7,7)) 
    plt.plot(t_m, v_m[:,0], "y", label="$V_{x}$")
    plt.plot(t_m, v_m[:,1], "r", label = "$V_{y}$")
    plt.plot(t_m, v_m[:,2], "g",label = "$V_{z}$")
    plt.title("Oppg 2 - $\delta t = 1fS \vec{V_{1}}$")
    plt.xlabel("Tid $pS$")
    plt.ylabel("Fart m/s")
    plt.legend()
    plt.show()
    
    fig5 =p.figure(figsize=plt.figaspect(0.5)*1.5)
    ax = fig5.gca(projection='3d')
    ax.plot(r_m[:,0],r_m[:,1],r_m[:,2])
    ax.set_xlabel("X [m]"), ax.set_ylabel("Y [m]"),ax.set_zlabel("Z [m]")   
    p.title("Oppg 2 3D: $\delta t = 1fS$")    
    p.show()
    
    fig6 = p.figure("Oppgave 2,sammenligning analytisk",figsize=(7,7))
    ax = fig6.gca(projection='3d')
    ax.plot(r_m[:,0],r_m[:,1],r_m[:,2])
    #ax.plot(losnto)
    ax.set_xlabel("X[m]"), ax.set_ylabel("Y[m]"),ax.set_zlabel("Z[m]")
    p.title("Oppg 2 - $\delta t = 1fS \vec{V_{2}}$")
    p.show()
    
    solveequation(e,m_e,r_m,v_m,t_m,dt3)
   
            
