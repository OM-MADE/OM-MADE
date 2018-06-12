# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 11:51:16 2017

@author: AJ Tinet & P Collon

Main Program to validate transport in two mobile zones WITH exchange 
Input files are in the Comparison_DADE directory
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
#               DADE DATA READING
#
# =================================================================
os.chdir("..\Validations")
# The DADE software is used to generate reference data

filename = "Comparison_DADE\DADE_RefData.txt"
pfile = open(filename,'r')
DADE_main = np.zeros((100,3))
DADE_sec = np.zeros((100,3))
DADE_CV = np.zeros((100,3))
DADE_CE = np.zeros((100,3))
DADE_t = np.zeros(100)

it = -1
for ligne in pfile:
    it += 1
    mots = ligne.split('\t')
    DADE_t[it] = float(mots[0])
    for ix in range(3):
        DADE_main[it,ix] = float(mots[4*ix+1])
        DADE_sec[it,ix] = float(mots[4*ix+2])
        DADE_CV[it,ix] = float(mots[4*ix+3])
        DADE_CE[it,ix] = float(mots[4*ix+4])
        
pfile.close()
            
# =================================================================
#
#         NUMERICAL SIMULATION AND RESULT SAVING 
#
# =================================================================
            
simfile = "Comparison_DADE\DADE_simulation.txt" 
datafile = "Comparison_DADE\DADE_parameters.txt"
boundfile = "Comparison_DADE\DADE_boundary.txt"

print("Initialisation...")

dx,nx,dt,tmax,c0,Xprt,Tprt = readGeneralData(simfile)
bound = readBound(boundfile)
dataset, C, points = initialise(datafile,dt,dx, c0, nx)

print("Start Computing...")
dataobs = timeloop(points, C, dataset, nx, bound, dx, dt, tmax, Xprt, Tprt)

ne = len(dataset)-1
for ie in range(ne):
    np.save("Comparison_DADE\DADE_Results_C"+str(ie),dataobs[ie])

# =================================================================
#
#           PLOTS OF RESULTS
#
# =================================================================

locs = list(Xprt)
ic = -1

# Get the flux in both mobile zones
q1 = dataset[1][-1][0]
q2 = dataset[2][-1][0]

for ix in range(len(Xprt)):

    x = Xprt[ix]
    ic += 2
    
    # Plot DADE results
    plt.plot(DADE_t,DADE_CE[:,ic//2],"o"+colors[ic])
    
    # Plot simulation results
    plt.plot(Tprt,(q1*dataobs[0][ix,:]+q2*dataobs[1][ix,:])/(q1+q2),colors[ic%len(colors)]+"--",label="OM-MADE "+str(x))
    

    
plt.legend(loc='best')
plt.xlabel("Time (days)")
plt.ylabel("Concentration (micro-g/l)")
plt.title("Average concentration values (o DADE    -- OM-MADE)")
plt.show()




