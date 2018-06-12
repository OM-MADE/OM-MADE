# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 11:51:16 2017

@author: AJ Tinet & P Collon

Main Program to validate Pure Dispersion 
Input files are in a PureDispersion directory
"""

import numpy as np
import matplotlib.pyplot as plt
from math import erfc,sqrt

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

# Analytical solution of th transport by dispersion only in case of 
# a continuous injection of Ci = 350 g/m3, with a dispersion coefficient D= 0.5 m²/s:
# 
# C = Ci * erf(x/2*rac(Dt))

dispersion_th = np.zeros((1501,9))
d = 0.05
dispersion_th[0,0]=350

for ix in range(1501):
    for it in range(9):
        x = ix
        t = it*150*3600
        
        if it != 0:
            dispersion_th[ix,it] = 350*erfc(x/(2*sqrt(d*t)))
            
# =================================================================
#
#         NUMERICAL SIMULATION AND RESULT SAVING
#
# =================================================================
os.chdir("..\Validations")
            
simfile = "PureDispersion\PureDispersion_simulation.txt" 
datafile = "PureDispersion\PureDispersion_parameters.txt"
boundfile = "PureDispersion\PureDispersion_boundary.txt"

print("Initialisation...")

dx,nx,dt,tmax,c0,Xprt,Tprt = readGeneralData(simfile)
bound = readBound(boundfile)
dataset, C, points = initialise(datafile,dt,dx, c0, nx)

print("Start Computing...")
dataobs = timeloop(points, C, dataset, nx, bound, dx, dt, tmax, Xprt, Tprt)


np.save("PureDispersion\PureDispersion_Results_C",dataobs[0])

# =================================================================
#
#         PLOT OF COMPARISON THEORY-SIMULATION 
#
# =================================================================

locs = list(Xprt)
ic = 0

for it in range(1,len(Tprt)):

    t = Tprt[it]
    ic += 1
    
    # Plot of analytical results
    plt.plot(Xprt[::50],dispersion_th[::50,it],colors[ic]+"o",label="Theory "+str(int(t/3600))+"h")
    
    # Plot of simulation results
    plt.plot(Xprt,dataobs[0][:,it],colors[ic%len(colors)]+"--",label="OM-MADE "+str(int(t/3600))+"h")
    
#plt.legend(loc='best')
plt.xlabel("Location (m)")
plt.ylabel("Concentration (g/m3)")
plt.title("Pure dispersion (o Theory    -- OM-MADE)")
plt.show()



