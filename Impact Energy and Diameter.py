# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 12:49:16 2021

@author: Deepayan Banik
"""

""" Wherever not mentioned, units are SI"""
import numpy as np
import matplotlib.pyplot as plt
rho_a = 1260 # density of Bennu
D_a = 520 # diameter of Bennu
R_a = D_a / 2 / 1000 # radius in Km

G = 6.67408e-11
f_char = 0.4 # characteristic frequency for Bennu, Quillen
f = 0.41 * f_char # Quillen
E_s = np.pi * G**2 * rho_a **3 * D_a **5 / 108 / f**2 # seismic energy for 1 g_a acc: 366 Joules per impact

rho_p = 2700 # density of impactor
v_p = 10e3 # Rich 2020 Impact produced seismic shaking and regolith growth on asteroids 433 Eros, 2867 Šteins, and 25143 Itokawa
n = 1e-7 # Rich 2020

# min impactor diam for 1g_a acc
D_p = np.cbrt(G**2 * rho_a **3 * D_a **5/ 9/ n / rho_p / v_p**2 / f**2) / 1000 # in Km, Rich 2005

# total number of impactors in the main-belt having a diameter greater than d (in km)
N = 2.8e6 * D_p**(-2.2)  # Bottke 2005
# no of particles greater than D_p that hits Bennu over it's mean collisional lifetime 'mcl_b'
# mcl_b = 1.75e6 # yr (assumed same as Bennu's NEO lifetime) Balouz 2020
mcl_b = 65e6 # yr (assumed same as Bennu's MBA lifetime) Balouz 2020
#mcl_b = 4.5e9 # yr (assumed same as Bennu's total lifetime)
pi = 2.9e-18 # intrinsic collisional probability km^-2 yr^-1 Balouz 2020, Bottke 2005
n_impact = pi * R_a**2 * N * mcl_b # total number of impacts in mean collisional life
print(n_impact)
''' There are approx. 34 impactors that hit Bennu and are able to cause global
regolith mobilization (i.e. of size greater than 37 cm) during its NEO lifetime'''

#%%

initial_core_slant_height = 500
epsilon = 0.01
number_of_landslides = 14#np.linspace(1,250,250)
final_core_slant_height = (1-2*epsilon)**number_of_landslides * initial_core_slant_height
# final_core_slant_tweek = 0*number_of_landslides * initial_core_slant_height + 145
# plt.plot(number_of_landslides,final_core_slant_height)
# plt.plot(number_of_landslides,final_core_slant_tweek)
# plt.xlabel('No. of landslides')
# plt.ylabel('Final core slant height')