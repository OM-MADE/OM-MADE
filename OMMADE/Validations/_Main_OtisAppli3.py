# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 11:51:16 2017

@author: AJ Tinet & P Collon

Main Program to validate transport with first-order degradation rate
Input files are in the Comparison_OTIS_App3 directory.
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
filename = "Comparison_OTIS_App3\OTIS_App3_RefData.txt"
pfile = open(filename,'r')
app3 = np.zeros((303,2))
app3_t = np.zeros(303)

it = -1
for ligne in pfile:
    it += 1
    mots = ligne.split('\t')
    app3_t[it] = float(mots[0])*3600
    for ix in range(2):
        app3[it,ix] = float(mots[ix+1])
        
pfile.close()
            
# =================================================================
#
#          NUMERICAL SIMULATION AND RESULT SAVING
#
# =================================================================
            
simfile = "Comparison_OTIS_App3\OTIS_App3_simulation.txt" 
datafile = "Comparison_OTIS_App3\OTIS_App3_parameters.txt"
boundfile = "Comparison_OTIS_App3\OTIS_App3_boundary.txt"

print("Initialisation...")

dx,nx,dt,tmax,c0,Xprt,Tprt = readGeneralData(simfile)
bound = readBound(boundfile)
dataset, C, points = initialise(datafile,dt,dx, c0, nx)

print("Start Computing...")
dataobs = timeloop(points, C, dataset, nx, bound, dx, dt, tmax, Xprt, Tprt)


np.save("Comparison_OTIS_App3\OTIS_App3_Results_C",dataobs[0])

# =================================================================
#
#          PLOTS OF RESULTS
#
# =================================================================

ix = 0
ic = 0

x = Xprt[ix]
rmse = []


#Plot at location x=100m    

# Plot OTIS Results
plt.plot(app3_t,app3[:,ic],"o"+colors[ic], label="OTIS")
    
# Plot simulation results
plt.plot(Tprt,dataobs[0][ix,:],colors[ic]+"--",label="OM-MADE")
    
  
plt.legend(loc='best')
plt.xlabel("Time (s)")
plt.ylabel("Concentration (g/m3)")
plt.title("X = 100 m")
plt.show()

somme = 0
mean = 0
    
for i in range(len(Tprt)):
    
    somme += (app3[i,ic] - dataobs[0][ix,i])**2
    mean += app3[i,ic]
        

somme = somme**0.5 / mean
    
rmse.append(somme)

ix = 1
ic = 1

x = Xprt[ix]
 
#Plot at location x=2000m    
# Plot OTIS Results
plt.plot(app3_t[::5],app3[::5,ic],"o"+colors[ic], label="OTIS")
    
# Plot simulation results
plt.plot(Tprt,dataobs[0][ix,:],colors[ic]+"--",label="OM-MADE")
    
  
plt.legend(loc='best')
plt.xlabel("Time (s)")
plt.ylabel("Concentration (g/m3)")
plt.title("X = 2000 m")
plt.show()

somme = 0
mean = 0
    
for i in range(len(Tprt)):
    
    somme += (app3[i,ic] - dataobs[0][ix,i])**2
    mean += app3[i,ic]
        
somme /= mean
    
rmse.append(somme)

plt.plot([100, 2000], rmse, "k*")
plt.yscale('log')
plt.xlabel("Distance (m)")
plt.ylabel("NRMSE")
plt.show()







