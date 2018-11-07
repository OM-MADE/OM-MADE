# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 11:51:16 2017

@author: AJ Tinet & P Collon

Main Program to validate transport in multi-reach transport zone exchanging 
with storage zone 
Input files are in the Comparison_OTIS_App1 directory
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
#               OTIS DATA READING
#
# =================================================================
os.chdir("..\Validations")

filename = "Comparison_OTIS_App1\OTIS_App1_RefData.txt"
pfile = open(filename,'r')
app1_main = np.zeros((160,5))
app1_stock = np.zeros((160,5))
app1_t = np.zeros(160)

it = -1
for ligne in pfile:
    it += 1
    mots = ligne.split('\t')
    app1_t[it] = float(mots[0])
    for ix in range(5):
        app1_main[it,ix] = float(mots[ix+1])
        app1_stock[it,ix] = float(mots[ix+6])
        
pfile.close()
            
# =================================================================
#
#         NUMERICAL SIMULATION AND RESULT SAVING 
#
# =================================================================
            
simfile = "Comparison_OTIS_App1\OTIS_App1_simulation.txt" 
datafile = "Comparison_OTIS_App1\OTIS_App1_parameters.txt"
boundfile = "Comparison_OTIS_App1\OTIS_App1_boundary.txt"

print("Initialisation...")

dx,nx,dt,tmax,c0,Xprt,Tprt = readGeneralData(simfile)
bound = readBound(boundfile)
dataset, C, points = initialise(datafile,dt,dx, c0, nx)

print("Start Computing...")
dataobs = timeloop(points, C, dataset, nx, bound, dx, dt, tmax, Xprt, Tprt)


ne = len(dataset)-1
for ie in range(ne):
    np.save("Comparison_OTIS_App1\OTIS_App1_Results_C"+str(ie),dataobs[ie])

# =================================================================
#
#         PLOTS OF RESULTS 
#
# =================================================================

ic = -1

#Plots of concentration in main channel
for ix in range(len(Xprt)):

    x = Xprt[ix]
    if x!=433:
        ic += 1
        # Plot OTIS Results
        plt.plot(app1_t,app1_main[:,ic],"o"+colors[ic])
        
        # Plot simulation results
        plt.plot(Tprt,dataobs[0][ix,:],colors[ic]+"--",label="OM-MADE "+str(x))
    else:
        plt.plot(Tprt,dataobs[0][ix,:],colors[ic]+"-",label="OM-MADE "+str(x))

    
plt.legend(loc='best')
plt.xlabel("Time (s)")
plt.ylabel("Concentration (mg/l)")
plt.title("Main Chanel (o OTIS    -- OM-MADE)")
plt.show()


#Plots of concentration in storage zone
ic = 1

for ix in range(2,len(Xprt)):

    x = Xprt[ix]
    ic += 1
    if x!=432:
            
        # Plot OTIS Results
        plt.plot(app1_t[::4],app1_stock[::4,ic],"o"+colors[ic])
        
        # Plot simulation results
        plt.plot(Tprt[:],dataobs[1][ix,:],colors[ic]+"--",label="OM-MADE "+str(x))
    else:
        plt.plot(Tprt[:],dataobs[1][ix,:],colors[ic]+"-",label="OM-MADE "+str(x))
        ic-=1

    
plt.legend(loc='best')
plt.xlabel("Time (s)")
plt.ylabel("Concentration (mg/l)")
plt.title("Storage (o OTIS    -- OM-MADE)")
plt.show()





