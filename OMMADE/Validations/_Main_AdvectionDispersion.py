# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 11:51:16 2017

@author: AJ Tinet & P Collon

Main Program to validate transport by Advection & Dispersion 
Input files are in a Advection_Dispersion directory
"""

import numpy as np
import matplotlib.pyplot as plt
from math import exp,sqrt,pi

import os
os.chdir("..\Codes_OMMADE")
from readData import *
from timeLoops import *
from classDataPoint import *
from classParameters import *
from time import time

colors = ['k','b','r','g','c','m','y','k','b','r','g','c','m','y']

# =================================================================
#
#               THEORETICAL DATA WRITING
#
# =================================================================
os.chdir("..\Validations")

# Analytical solution analytique of advection - dispersion transport in the case
# of a dirac injection of M/A = 900 g/m, with a dispersion coeffcition D =0.05 m²/s 
# with V= 0.01m/s:
# 
# C = Ci * erf(x/2*rac(Dt))

theory = np.zeros((1501,9))
d = 0.0025
u = 0.01
ma = 900

for ix in range(1501):
    for it in range(9):
        x = ix
        t = it*5*3600
        
        if it != 0:
            theory[ix,it] = (ma/sqrt(4*pi*d*t))*exp(-(x-u*t)**2/(4*d*t))
            
# =================================================================
#
#         NUMERICAL SIMULATION AND RESULT SAVING
#
# =================================================================
            
simfile = "Advection_Dispersion\AdvectionDispersion_simulation.txt" 
datafile = "Advection_Dispersion\AdvectionDispersion_parameters.txt"
boundfile = "Advection_Dispersion\AdvectionDispersion_boundary.txt"

print("Initialisation...")

dx,nx,dt,tmax,c0,Xprt,Tprt,scheme = readGeneralData(simfile)
bound = readBound(boundfile)
dataset, C, points = initialise(datafile,dt,dx, c0, nx, scheme)

print("Start Computing...")
t0 = time()
dataobs = timeloop(points, C, dataset, nx, bound, dx, dt, tmax, Xprt, Tprt, scheme)
t1 = time()
print(t1-t0)

np.save("Advection_Dispersion\AdvectionDispersion_Results_C",dataobs[0])

# =================================================================
#
#         PLOT OF COMPARISON THEORY-SIMULATION 
#
# =================================================================

locs = list(Xprt)
ic = -1

rmse= []

for it in range(1,len(Tprt)):

    t = Tprt[it]
    ic += 1
    
    # Plot of analytical results
    plt.plot(Xprt[::25],theory[::25,it],colors[ic]+"o",label="Theory "+str(int(t/3600))+"h")
    
    # Plot of simulation results
    plt.plot(Xprt,dataobs[0][:,it],colors[ic%len(colors)]+"--",label="OM-MADE "+str(int(t/3600))+"h")
    
    somme = 0
    mean = 0
    
    for i in range(len(Xprt)):
        
        somme += (theory[i,it] - dataobs[0][i,it])**2
        mean += theory[i,it]
        
    somme = somme**0.5 / mean
    
    rmse.append(somme)
    
#plt.legend(loc='best')
plt.xlabel("Location (m)")
plt.ylabel("Concentration (g/m3)")
plt.title("Advection dispersion (o Theory    -- OM-MADE)")
plt.show()

plt.plot((1/3600)*Tprt[1:], rmse, "k*")
plt.yscale('log')
plt.xlabel("Time (h)")
plt.ylabel("NRMSE")
plt.show()



