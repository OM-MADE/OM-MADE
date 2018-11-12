# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 11:51:16 2017

@author: Anne-Julie Tinet and Pauline Collon

Main Program for Furfooz Tracer test 3 (Dewaide et al., 2017, Hydrogeology journal)
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
#               EXPERIMENTAL DATA READING
#
# =================================================================

# Experimental data are stored in an independent file
# (1 file per location). These files store the 3 tracer tests (rhodamine and uranine)

# Each data is stored in a dictionary, with the location for key. 
# The value is a list concentration through time.
# The time_exp is also a dictionary, but it has for value a list of 2 elements:
#(i) the distance from the swallow hole
#(ii) the array of observation times
# Note that "TrouQuiFume" corresponds to "Site 1",
# "PuitsDesVaux" corresponds to "Site 3",
# "GalerieDesSources"corresponds to "Site 4",
# "TrouDeLaLoutre" corresponds to "Site 5", in the paper Dewaide et al., 2017

os.chdir("..\Furfooz_TracerTest3")

time_exp = {"TrouQuiFume":[150,np.zeros(1372)],\
             "GalerieDesSources":[770,np.zeros(12360)]}
rhoda_exp = {"TrouQuiFume":np.zeros(1372),\
             "GalerieDesSources":np.zeros(12360)}
ura_exp = {"TrouQuiFume":np.zeros(1372),\
             "GalerieDesSources":np.zeros(12360)}

file_exp = {"TrouQuiFume":"Tracage3_TQF.txt",\
             "GalerieDesSources":"Tracage3_GDS.txt"}
             
# File Reading

for loc in file_exp:
    filename = file_exp[loc]
    pfile = open(filename,'r')
    
    it = -1
    time_it = 0
    for ligne in pfile:
        it += 1
        mots = ligne.split('\t')
        
# TO READ time exp in the data file, but there is an apporxiamtion cause it is written in h. 
#In reality it was a regular step of 120s ( =0.0333 h)
        time_exp[loc][1][it] = float(mots[0])
#        time_exp[loc][1][it] = time_it
        if loc in ura_exp:
            ura_exp[loc][it] = float(mots[1])
            rhoda_exp[loc][it] = float(mots[2])
            
        else:
            rhoda_exp[loc][it] = float(mots[1])
        time_it+=120    
    pfile.close()
    
    
# =================================================================
#
#         NUMERICAL SIMULATION AND RESULT STORAGE
#
# =================================================================
            
simfile = "Tracage_3_simulation.txt" 
datafile = "Tracage_3_parameters.txt"
boundfile = "Tracage_3_boundary.txt"

print("Initialisation...")

dx,nx,dt,tmax,c0,Xprt,Tprt,scheme = readGeneralData(simfile)
bound = readBound(boundfile)
dataset, C, points = initialise(datafile,dt,dx, c0, nx, scheme)

print("Start Computing...")
dataobs = timeloop(points, C, dataset, nx, bound, dx, dt, tmax, Xprt, Tprt, scheme)

ne = len(dataset)-1
for ie in range(ne):
    np.save("Tracage3_Results_C"+str(ie),dataobs[ie])


# =================================================================
#
#         PLOT THE COMPARAISON FIELD - SIMULATION (rhodamine)
#
# =================================================================

locs = list(Xprt)
ic = -1

for location in time_exp:
    
    x = int(time_exp[location][0])
    ix = locs.index(x)
    ic += 1
    
    # Get the flows
    q1 = dataset[1][-1][0]
    q2 = dataset[2][-1][0]
    
    # Plot experimental data
    plt.plot(time_exp[location][1][::1],rhoda_exp[location][::1],colors[0]+".",label="Field")
                
    # Plot simulation results
    plt.plot(Tprt*(1/3600),(q1*dataobs[0][ix,:] + q2*dataobs[1][ix,:])/(q1+q2),colors[5]+"-",label="OM-MADE")
    #plt.plot(Tprt*(1/3600),dataobs[0][ix,:] ,colors[1]+"--",label="zone 1")
    #plt.plot(Tprt*(1/3600),dataobs[1][ix,:],colors[3]+"--",label="zone 2")
    
    if location == "TrouQuiFume":
        plt.xlim(0,10)
        
    if location == "GalerieDesSources":
        plt.xlim(0,600)
    
    plt.legend(loc='best')
    plt.xlabel("Time (h)")
    plt.ylabel("Concentration (ppb)")
    plt.title(location)
    plt.show()


