# -*- coding: utf-8 -*-
"""
Created on Sat May 23 13:31:27 2020

@author: Deepayan Banik
"""
import scipy.special
import numpy as np
import mpmath as mp
import math
import matplotlib.pyplot as plt
import seaborn as sns

res = 40
zet = 45 * np.pi / 180 # semi apex angle of the cone
HE = np.cos(zet) # total HE assuming slant height as 1
RAD = np.sin(zet) # total RAD
b = HE
a = 0

coordZ = np.linspace(HE - HE/res/2, HE/res/2, res) # centre points of each stacked disk of thickness HE / res decreasing with height
Tau = np.linspace(RAD/res/2, RAD - RAD/res/2, res) # mean radius of each ring of width RAD / res increasing with radius
MR = np.zeros(res) # radial component of gravity
MZ = np.zeros(res) # axial component of gravity
MAG = np.zeros(res) # total magnitude of gravity
ETA = np.zeros(res) # angle made by the gravity vector at any point on the surface of the core
POT = np.zeros(res) # potential at any point on the surface of the core


def my_R(z, R, Z): # z = vertical lower limit of current disc
    a = (HE - abs(z)) * np.tan(zet)
    smh = Z - z
    delta = np.sqrt((a + R)**2 + (smh)**2)
    k = 2 * np.sqrt(a * R) / delta
    ks = scipy.special.ellipk(k**2)
    es = scipy.special.ellipe(k**2)
    return (2 * delta * ((1 - k**2 / 2) * ks - es) / R)


def my_Z(z, R, Z):
    a = (HE - abs(z)) * np.tan(zet)
    if (R < a):
        epsilon = 1
    elif (R > a):
        epsilon = 0
    else:
        epsilon = 0.5
    smh = Z - z
    delta = np.sqrt((a + R)**2 + (smh)**2)
    k = 2 * np.sqrt(a * R) / delta
    m = 2 * np.sqrt(a * R) / (a + R)
    ks = scipy.special.ellipk(k**2)
    return (2 * np.pi * np.sign(smh) * epsilon + 2 * smh * ((R - a)/(R + a) * mp.ellippi(m**2, k**2) - ks) / delta)


def my_pot(z, R, Z):
    a = (HE - abs(z)) * np.tan(zet)
    if (R < a):
        epsilon = 1
    elif (R > a):
        epsilon = 0
    else:
        epsilon = 0.5
    smh = Z - z
    delta = np.sqrt((a + R)**2 + (smh)**2)
    k = 2 * np.sqrt(a * R) / delta
    m = 2 * np.sqrt(a * R) / (a + R)
    ks = scipy.special.ellipk(k**2)
    es = scipy.special.ellipe(k**2)
    return 2 * (-np.pi * abs(smh) * epsilon + delta * es + (a**2 - R**2) * ks / delta + smh**2 / delta * (- R + a)/(R + a) * mp.ellippi(m**2, k**2))

# def gauss(n,R,Z,fl):
#    G = 0
#    [z,w] = p_roots(n+1)
#    if (fl == 1):
#        for j in range (0,len(z)):
#            G = G + my_R((b - a) * z[j] / 2 + (b + a) / 2, R, Z) * w[j]
#    else:
#        for j in range (0,len(z)):
#            G = G + my_Z((b - a) * z[j] / 2 + (b + a) / 2, R, Z) * w[j]
#    return G * (b - a) / 2


def deep(res, R, Z, fl): # deep(resolution, radial surface location, vertical surface location, flag)
    D = 0 # dummy storage for returning magnitude
    dz = HE / res # thickness of discs
    if (fl == 1):
        for i in range(-res, res + 1):
            z = i * dz # adds up contributions from different discs
            D = D + dz * my_R(z, R, Z)
    elif (fl == 2):
        for i in range(-res, res + 1):
            z = i * dz
            D = D + dz * my_Z(z, R, Z)
    else:
        for i in range(-res, res + 1):
            z = i * dz
            D = D + dz * my_pot(z, R, Z)
    return D

om = 2.0
AR = np.zeros(res) # net radial gravity
MTAN = np.zeros(res) # tangential component of gravity magnitude
MNOR = np.zeros(res) # normal component of gravity magnitude\
ETAN = np.zeros(res) # tangential component of effective gravity magnitude
ENOR = np.zeros(res) # normal component of effective gravity magnitude
VESC = np.zeros(res) # escape velocity


# print(gauss(my_f,2,-1,1))
for l in range(0, res):
    POT[l] = deep(res, Tau[l], coordZ[l], 3) # potential
    MR[l] = abs(deep(res, Tau[l], coordZ[l], 1)) # radial magnitude of gravity
    MZ[l] = abs(deep(res, Tau[l], coordZ[l], 2)) # vertical magnitude of gravity
    MAG[l] = np.sqrt(MR[l]**2 + MZ[l]**2) # total magnitude
    ETA[l] = np.arctan(MR[l]/MZ[l]) # angle of gravity vector wrt horizontal

MAG = MAG / MAG[-1]

#%%

for l in range(0, res):
    AR[l] = abs(MR[l] - om * om * Tau[l])
    MTAN[l] = MAG[l] * np.cos(zet + ETA[l])
    MNOR[l] = MAG[l] * np.sin(zet + ETA[l]) # always positive
    VESC[l] = np.sqrt(2*POT[l])
    ETAN[l] = MTAN[l] + om * om * Tau[l] * np.sin(zet)
    ENOR[l] = - MNOR[l] + om * om * Tau[l] * np.cos(zet)

com_ZET = np.pi/2-zet # complement of the semi apex angle

OTAN = np.zeros(res) # angular velocity omega for zero tangential acceleration
ONOR = np.zeros(res) # angular velocity omega for zero normal acceleration
STAN = np.zeros(40) # downslope location for zero tangential acceleration
SNOR = np.zeros(40) # downslope location for zero normal acceleration
OMP = np.zeros((2,res)) # includes the effect of friction
VLO = np.zeros((2,res)) # lift off velocity
mu = np.tan(20 * np.pi / 180)
# omega = np.linspace(0,2,40)
# tt=0
# for om in omega:
for l in range(0, res): # stores gravity field at different locations along the slope, from the pole to the equator

	# take off angular velocity
    ONOR[l] = np.sqrt(MNOR[l] /np.cos(zet)/Tau[l])

	# slope equilibration angular velocity
    if (-MTAN[l]/np.sin(zet)/Tau[l]<0):
        OTAN[l] = math.nan
        if (-MTAN[l+1]/np.sin(zet)/Tau[l+1]>0):
             OTAN[l] = 0
    else:
        OTAN[l] = np.sqrt(-MTAN[l]/np.sin(zet)/Tau[l])

	# friction - static zones
    temp1 = (MTAN[l] - mu * MNOR[l])/Tau[l]/(- mu* np.cos(zet) - np.sin(zet))
    temp2 = (MTAN[l] + mu * MNOR[l])/Tau[l]/(mu* np.cos(zet) - np.sin(zet))
    if (temp1 < 0):
        OMP[0,l] = math.nan
        if ((MTAN[l+1] - mu * MNOR[l+1])/Tau[l+1]/(- mu* np.cos(zet) - np.sin(zet))>0):
             OMP[0,l] = 0
    else:
        OMP[0,l] = np.sqrt(temp1)
    if (temp2 < 0):
        OMP[1,l] = math.nan
        if ((MTAN[l+1] + mu * MNOR[l+1])/Tau[l+1]/(mu* np.cos(zet) - np.sin(zet))>0):
             OMP[1,l] = 0
    else:
        OMP[1,l] = np.sqrt(temp2)

	# lift-off velocitites
    VLO[0,l] = - om * Tau[l] - np.sqrt((om * Tau[l])**2-Tau[l] * ENOR[l]/np.cos(zet))
    VLO[1,l] = - om * Tau[l] + np.sqrt((om * Tau[l])**2-Tau[l] * ENOR[l]/np.cos(zet))

#    STAN[tt] = Tau[np.argmin(MTAN)] / np.sin(zet)
#    SNOR[tt] = Tau[np.argmin(MNOR)] / np.sin(zet)
#    tt = tt + 1

# %%

# colors = sns.color_palette("rocket", 3)
# fig, axs = plt.subplots(1,1, figsize=(5,5))
'''# fig = plt.figure(1, figsize=(5, 5))
#axh = fig.add_axes([0.35,0.25,0.6,0.5],xlim=(0,1))
#axh = fig.add_axes([0.175,0.175,0.7,0.7])
# axh.plot(omega,STAN)
# axh.plot(omega,SNOR)'''

# axs.plot(Tau/max(Tau), ONOR, linestyle ='-', marker ='^', markersize = 10, mfc = 'w', linewidth = 2, color=colors[0])
# axs.plot(Tau/max(Tau), OTAN, linestyle ='-', marker ='o', markersize = 10, mfc = 'w', linewidth = 2, color=colors[1])
# axs.plot(Tau/max(Tau), MNOR, linestyle ='-', marker ='o', markersize = 10, mfc = 'w', linewidth = 2, color=colors[1])
# axs.plot(Tau/max(Tau), MTAN + mu * MNOR, linestyle ='-', marker ='s', markersize = 10, mfc = 'w', linewidth = 2, color=colors[2])
plt.figure()
plt.plot(Tau/max(Tau),VLO[0,:])
plt.plot(Tau/max(Tau),VLO[1,:])
plt.plot(Tau/max(Tau),ENOR)
# plt.plot(Tau/max(Tau),OTAN)
# plt.plot(Tau/max(Tau),ONOR)

# plt.show()

'''# xticks = np.arange(0,1.01,0.25)
# yticks = np.arange(0, 0.26, 0.1)
# axs.xticks(xticks)
# axs.yticks(yticks)
# axs.xlim(0,1.1)
# axs.ylim(-0.1, 0.26)'''
# axs.minorticks_on()
# axs.tick_params(labelsize = 14)
# axs.tick_params(direction = 'in', right = True, top = True)
# axs.tick_params(direction ='in', which = 'minor', length = 5, right = True, top = True)
# axs.tick_params(direction ='in', which = 'major', length = 15, right = True, top = True)
# axs.grid()

'''# axs.legend()
# axs.plot(Tau/max(Tau),MTAN/max(MAG))
# axh.plot(Tau,MR/max(MAG))
# axh.plot(Tau,MZ/max(MAG))
# axh.plot(Tau,coordZ)'''

# axs.set_xlabel("downslope location", fontsize=14)
# axs.set_ylabel("effective acceleration", fontsize=14)
# for axis in ['top','bottom','left','right']:
# 	axs.spines[axis].set_linewidth(2)
# fig.tight_layout()

