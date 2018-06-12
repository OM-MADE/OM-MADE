# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 11:51:16 2017

@author: AJ Tinet & P Collon

Main Program to validate transport in two mobile zones without exchange 
Input files are in the Comparison_WSADE directory
"""

import numpy as np
import matplotlib.pyplot as plt

import os
os.chdir("..\Codes_OMMADE")
from readData import *
from timeLoops import *
from classDataPoint import *
from classParameters import *

colors = ['k','b','r','g','c','m','y','k','b','r','g','c','m','y']

# =================================================================
#
#               WSADE DATA READING
#
# =================================================================
os.chdir("..\Validations")
# The WSADE software is used to generate reference data

filename = "Comparison_WSADE\WSADE_RefData.txt"
pfile = open(filename,'r')
WSADE_main = np.zeros((101,3))
WSADE_sec = np.zeros((101,3))
WSADE_t = np.zeros(101)

it = -1
for ligne in pfile:
    it += 1
    mots = ligne.split('\t')
    WSADE_t[it] = float(mots[0])
    for ix in range(3):
        WSADE_main[it,ix] = float(mots[2*ix+1])
        WSADE_sec[it,ix] = float(mots[2*ix+2])
        
pfile.close()
            
# =================================================================
#
#         NUMERICAL SIMULATION AND RESULT SAVING
#
# =================================================================
            
simfile = "Comparison_WSADE\WSADE_simulation.txt" 
datafile = "Comparison_WSADE\WSADE_parameters.txt"
boundfile = "Comparison_WSADE\WSADE_boundary.txt"

print("Initialisation...")

dx,nx,dt,tmax,c0,Xprt,Tprt = readGeneralData(simfile)
bound = readBound(boundfile)
dataset, C, points = initialise(datafile,dt,dx, c0, nx)

print("Start Computing...")
dataobs = timeloop(points, C, dataset, nx, bound, dx, dt, tmax, Xprt, Tprt)

ne = len(dataset)-1
for ie in range(ne):
    np.save("Comparison_WSADE\WSADE_Results_C"+str(ie),dataobs[ie])

# =================================================================
#
#         PLOTS OF RESULTS
#
# =================================================================

locs = list(Xprt)
ic = -1

for ix in range(len(Xprt)):

    x = Xprt[ix]
    ic += 2
    
    # Plot of WSADE results
    plt.plot(WSADE_t,WSADE_main[:,ic//2],"o"+colors[ic])
    plt.plot(WSADE_t,WSADE_sec[:,ic//2],"o"+colors[ic+1])
    
    # Plot of simulation results
    plt.plot(Tprt,dataobs[0][ix,:],colors[ic%len(colors)]+"--",label="OM-MADE "+str(x)+"-C0")
    plt.plot(Tprt,dataobs[1][ix,:],colors[(ic+1)%len(colors)]+"--",label="OM-MADE "+str(x)+"-C1")
    

    
plt.legend(loc='best')
plt.xlabel("Time (days)")
plt.ylabel("Concentration (micro-g/l)")
plt.title("Local values (o WSADE    -- OM-MADE)")
plt.show()



