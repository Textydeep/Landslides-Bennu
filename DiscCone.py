# -*- coding: utf-8 -*-
"""
Created on Sat May 23 13:31:27 2020

@author: Deepayan Banik
"""

# -*- coding: utf-8 -*-
"""
Created on Sat May 23 12:22:48 2020

@author: Deepayan Banik
"""
import scipy.special
import numpy as np
import mpmath as mp
import math
#from sympy import elliptic_pi
import matplotlib.pyplot as plt
#from scipy.special.orthogonal import p_roots

res = 200
phi = 45 * np.pi / 180
HE = np.cos(phi) # total HE
RAD = np.sin(phi) # total RAD
b = HE
a = 0

coordZ = np.linspace(HE - HE/res/2,HE/res/2,res)
coordR = np.linspace(RAD/res/2 ,RAD - RAD/res/2,res)
MR = np.zeros(res)
MZ = np.zeros(res)
MAG = np.zeros(res)
ETA = np.zeros(res)
POT = np.zeros(res)

def my_R(z, R, Z):
    a = (HE - abs(z)) * np.tan(phi)
    zeta = Z - z
    delta = np.sqrt((a + R)**2 + (zeta)**2)
    k = 2 * np.sqrt(a * R) / delta
    ks=scipy.special.ellipk(k**2) 
    es=scipy.special.ellipe(k**2) 
    return (2 * delta * ((1- k**2 / 2) * ks - es) / R)

def my_Z(z, R, Z):
    a = (HE - abs(z)) * np.tan(phi)
    if (R < a):
        epsilon = 1
    elif (R > a):
        epsilon = 0
    else:
        epsilon = 0.5
    zeta = Z - z
    delta = np.sqrt((a + R)**2 + (zeta)**2)
    k = 2 * np.sqrt(a * R) / delta
    m = 2 * np.sqrt(a * R) / (a + R)
    ks=scipy.special.ellipk(k**2) 
    return (2 * np.pi * np.sign(zeta) * epsilon + 2 * zeta * ((R - a)/(R + a) * mp.ellippi(m**2,k**2) - ks) / delta)

def my_pot(z, R, Z):
    a = (HE - abs(z)) * np.tan(phi)
    if (R < a):
        epsilon = 1
    elif (R > a):
        epsilon = 0
    else:
        epsilon = 0.5
    zeta = Z - z
    delta = np.sqrt((a + R)**2 + (zeta)**2)
    k = 2 * np.sqrt(a * R) / delta
    m = 2 * np.sqrt(a * R) / (a + R)
    ks=scipy.special.ellipk(k**2) 
    es=scipy.special.ellipe(k**2) 
    return 2 * (-np.pi * abs(zeta) * epsilon + delta * es + (a**2 - R**2) * ks / delta + zeta**2 / delta * (- R + a)/(R + a) * mp.ellippi(m**2,k**2))

#def gauss(n,R,Z,fl):
#    G = 0
#    [z,w] = p_roots(n+1)
#    if (fl == 1):
#        for j in range (0,len(z)):
#            G = G + my_R((b - a) * z[j] / 2 + (b + a) / 2, R, Z) * w[j]
#    else:
#        for j in range (0,len(z)):
#            G = G + my_Z((b - a) * z[j] / 2 + (b + a) / 2, R, Z) * w[j]
#    return G * (b - a) / 2

def deep(n,R,Z,fl):
    D = 0
    dz = HE / n
    if (fl == 1):
        for i in range (-n-1,n):
            z = (i + 1) * dz
            D = D + dz * my_R(z,R,Z)
    elif (fl == 2):
        for i in range (-n-1,n):
            z = (i + 1) * dz
            D = D + dz * my_Z(z,R,Z)
    else:
        for i in range (-n-1,n):
            z = (i + 1) * dz
            D = D + dz * my_pot(z,R,Z)
    return D

#print(gauss(my_f,2,-1,1))
for l in range (0,len(coordR)):
    POT[l] = deep(200,coordR[l],coordZ[l],3)
    MR[l] = abs(deep(200,coordR[l],coordZ[l],1))
    MZ[l] = abs(deep(200,coordR[l],coordZ[l],2))
#    MR[l] = abs(deep(20,coordR[l],coordZ[l],3) - deep(20,coordR[l]+0.01,coordZ[l],3))/0.1
#    MZ[l] = abs(deep(20,coordR[l],coordZ[l],3) - deep(20,coordR[l],coordZ[l]+0.01,3))/0.1
    MAG[l] = np.sqrt(MR[l]**2 + MZ[l]**2)
    ETA[l] = np.arctan(MR[l]/MZ[l])
    
f = open("grav45.txt", "w")
g = open("eta45.txt", "w")

ZETA = np.pi/2-phi
AR = np.zeros(res)
MTAN = np.zeros(res)
MNOR = np.zeros(res)
OTAN = np.zeros(res)
ONOR = np.zeros(res)
STAN = np.zeros(40)
SNOR = np.zeros(40)
OMP = np.zeros(res)
VLO = np.zeros(res)
VESC = np.zeros(res)
mu = np.tan(20 * np.pi / 180)
#omega = np.linspace(0,2,40)
#tt=0
#for om in omega:
for l in range (0,len(coordR)):
#    f.write(str(MAG[l]/max(MAG)) + ",")
#    g.write(str(ETA[l]) + ",")
    VESC[l] = np.sqrt(2*POT[l])
    om = 1.0
#    AR[l] = abs(MR[l] / max(MAG) - om * om * coordR[l])
    MTAN[l] = (MAG[l] / max(MAG) * np.sin(ZETA - ETA[l])) + om * om * coordR[l] * np.cos(ZETA)
    MNOR[l] = (-MAG[l] / max(MAG) * np.cos(ZETA - ETA[l])) + om * om * coordR[l] * np.sin(ZETA)
#    ONOR[l] = np.sqrt(-MNOR[l]/np.cos(phi)/coordR[l])
#    OMP[l] = np.sqrt((MTAN[l] - mu * MNOR[l])/coordR[l]/(mu* np.sin(ZETA) - np.cos(ZETA)))
#    OMP[l] = np.sqrt((-MTAN[l] - mu * MNOR[l])/coordR[l]/(mu* np.sin(ZETA) + np.cos(ZETA)))
    VLO[l] = -om * coordR[l] - np.sqrt((om * coordR[l])**2-coordR[l] * MNOR[l]/np.cos(phi))
#    if ((MTAN[l] - mu * MNOR[l])>0):
#        OMP[l] = 0
#    else:
#        OMP[l] = np.sqrt((MTAN[l] - mu * MNOR[l])/coordR[l]/(mu* np.sin(ZETA) - np.cos(ZETA)))
#    if ((MR[l] * np.cos(ZETA) - MZ[l] * np.sin(ZETA))<0):
#        OTAN[l] = 0
#    else:
#    OTAN[l] = np.sqrt(-MTAN[l]/np.sin(phi)/coordR[l])
#    STAN[tt] = coordR[np.argmin(MTAN)] / np.sin(phi)
#    SNOR[tt] = coordR[np.argmin(MNOR)] / np.sin(phi)    
#    tt = tt + 1
#for l in range (0,len(coordR)):
#    if (OTAN[l]!=0):
#        OTAN[l-1] = 0
#        break;
#    else:
#        OTAN[l] = math.nan      
#for l in range (0,len(coordR)):
#    if (OMP[l]!=0):
#        OMP[l-1] = 0
#        break;
#    else:
#        OMP[l] = math.nan
#for m in range (0,len(coordR)):
#    l = len(coordR) - m - 1
#    if (OMP[l]!=0):
#        OMP[l+1] = 0
#        break;
#    else:
#        OMP[l] = math.nan

#f.close
#g.close

fig = plt.figure(figsize=(3.5,1.75))
axh = fig.add_axes([0.35,0.25,0.6,0.5],xlim=(0,1))
#axh = fig.add_axes([0.175,0.175,0.7,0.7])
#axh.plot(omega,STAN)
#axh.plot(omega,SNOR)

plt.plot(coordR/max(coordR),MTAN/max(MAG))
plt.plot(coordR/max(coordR),0*MNOR/max(MAG))
plt.plot(coordR/max(coordR),MTAN/max(MAG) + mu*MNOR/max(MAG))
#plt.plot(coordR/max(coordR),ONOR)

#plt.plot(coordR/max(coordR),MTAN/max(MAG))
#axh.plot(coordR,MR/max(MAG))
#axh.plot(coordR,MZ/max(MAG))
#axh.plot(coordR,coordZ)
plt.xlabel("s")
plt.ylabel("effective acceleration")