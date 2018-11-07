# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 11:51:16 2017

@author: AJ Tinet & P Collon

Main Program to validate Pure Advection 
Input files are in a PureAdvection directory
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

import os
os.chdir("..\Codes_OMMADE")
from readData import *
from timeLoops import *
from classDataPoint import *
from classParameters import *

colors = ['k','b','r','g','c','m','y','k','b','r','g','c','m','y']

# =================================================================
#
#               THEORETICAL DATA WRITING
#
# =================================================================

# Flow is V=0.01m/s. The injection step corresponds to a distance of 
# 108m starting in x = V*t - 108.
# Results are observed every 5h.

advection_th = np.zeros((1501,9))
v = 0.01
d = 3*3600
dt = 300

for ix in range(1501):
    for it in range(9):
        x = ix
        t = it*5*3600
        
        if v*(t-d) <= x <= v*(t+dt-1):
            advection_th[ix,it] = 350
            
# =================================================================
#
#         NUMERICAL SIMULATION AND RESULT SAVING
#
# =================================================================

os.chdir("..\Validations")
#Read Data in the PureAdvection folder            
simfile = "PureAdvection\PureAdvection_simulation.txt" 
datafile = "PureAdvection\PureAdvection_parameters.txt"
boundfile = "PureAdvection\PureAdvection_boundary.txt"

print("Initialisation...")

dx,nx,dt,tmax,c0,Xprt,Tprt = readGeneralData(simfile)
bound = readBound(boundfile)
dataset, C, points = initialise(datafile,dt,dx, c0, nx)

print("Start Computing...")
dataobs = timeloop(points, C, dataset, nx, bound, dx, dt, tmax, Xprt, Tprt)


np.save("PureAdvection\PureAdvection_Results_C",dataobs[0])

# =================================================================
#
#         PLOT OF COMPARISON THEORY-SIMULATION 
#
# =================================================================

locs = list(Xprt)
ic = -1
rmse = []

for it in range(len(Tprt)):

    t = Tprt[it]
    ic += 1
    
    # Plot of analytical results
    plt.plot(Xprt,advection_th[:,it],colors[ic],label="Theory "+str(int(t/3600))+"h")
    
    # Plot of simulation results
    plt.plot(Xprt,dataobs[0][:,it],colors[ic%len(colors)]+"--",label="OM-MADE "+str(int(t/3600))+"h")
    
    if it > 0:
        somme = 0
        mean = 0
        
        for i in range(len(Xprt)):
            
            somme += (advection_th[i,it] - dataobs[0][i,it])**2
            mean += advection_th[i,it]
            
        somme = somme**0.5 / mean
        
        rmse.append(somme)
    
#plt.legend(loc='best')
plt.xlabel("Location (m)")
plt.ylabel("Concentration (g/m3)")
plt.title("Pure advection (- theory    -- OM-MADE)")
plt.show()

plt.plot((1/3600)*Tprt[1:], rmse, "k*")
plt.yscale('log')
plt.xlabel("Time (h)")
plt.ylabel("NRMSE")
plt.show()


