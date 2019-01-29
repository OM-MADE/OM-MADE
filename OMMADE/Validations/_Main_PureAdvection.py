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
from time import time

import os
os.chdir("..\Codes_OMMADE")
from readData import *
from timeLoops import *
from classDataPoint import *
from classParameters import *

colors = ['k','b','r','g','c','m','y','k','b','r','g','c','m','y']
            
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

dx,nx,dt,tmax,c0,Xprt,Tprt,scheme = readGeneralData(simfile)
bound = readBound(boundfile)
dataset, C, points = initialise(datafile,dt,dx, c0, nx, scheme)

print("Start Computing...")
time0 = time()
dataobs = timeloop(points, C, dataset, nx, bound, dx, dt, tmax, Xprt, Tprt, scheme)
print(time() - time0)

np.save("PureAdvection\PureAdvection_Results_C",dataobs[0])


# =================================================================
#
#               THEORETICAL DATA WRITING
#
# =================================================================

# Flow is V=0.01m/s. The injection step corresponds to a distance of 
# 108m starting in x = V*t - 108.
# Results are observed every 5h.


advection_th = np.zeros((len(Xprt),len(Tprt)))

# Retrieving parameters
flow = dataset[1][1][0] # Flow rate
v = flow / dataset[1][0].getArea() # Velocity

# Describing the injection step
tstart = None
tstop = 0
cinject = max(bound[1][0])
for i in range(len(bound[0])):

    if bound[1][0][i] == cinject and tstart == None:
        tstart = bound[0][i]

    if bound[1][0][i] == 0 and tstart != None:
        break
    
    tstop = bound[0][i]

ix = 0
for x in Xprt:
    it = 0
    for t in Tprt:
        
        if v*(t - tstart - tstop) <= x <= v*(t - tstart + dt):
            advection_th[ix,it] = cinject
    
        it +=1
        
    ix += 1

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


