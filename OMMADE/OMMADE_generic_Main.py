# -*- coding: utf-8 -*-
"""
Created on Wen Oct 25 15:33:51 2017

@author: Pauline Collon and Anne-Julie Tinet

Main Program for an executable version of the program.
Directly export the results in Ascii Files.
Several files are generated : 
1) one per mobile zone, providing timesteps (first column) and corresponding 
    concentrations for each printing locations
2) a global file, providing timesteps (first column) and corresponding average 
    concentrations (a weighted sum of all mobile zone concentrations) for each 
    printing location. This file is called: RESULT_AverageFlow.txt and saved
    in the directory indicated in _INPUFILES.txt
"""

import numpy as np

import os
os.chdir("Codes_OMMADE")
from readData import *
from timeLoops import *
from classDataPoint import *
from classParameters import *

  
# =================================================================
#
#         DATA READING AND NUMERICAL SIMULATION
#
# =================================================================

os.chdir("..")
infile = open("_INPUTFILES.txt",'r')

path = lire_ligne(infile).rstrip('\n')
#Read the input file names in the _INPUTFILES file            
simfile = path + lire_ligne(infile).rstrip('\n')
boundfile = path + lire_ligne(infile).rstrip('\n')
datafile = path + lire_ligne(infile).rstrip('\n')


#Read the Data in the different input files
print("Initialisation...")

dx,nx,dt,tmax,c0,Xprt,Tprt,scheme = readGeneralData(simfile)
bound = readBound(boundfile)
dataset, C, points = initialise(datafile,dt,dx, c0, nx, scheme)

#Computation
print("Start Computing...")
dataobs = timeloop(points, C, dataset, nx, bound, dx, dt, tmax, Xprt, Tprt, scheme)

infile.close()

# =================================================================
#
#         EXPORT SIMULATION RESULTS
#
# =================================================================

#To save in txt format
#One file for each zone, first column is Time-Printed steps
# Following columns are concentration of the zone for each X-printed steps
ne = len(dataset)-1 #number of zones
nx = len (Xprt) #number of printed locations

for ie in range(ne): #for each zone
    res_zone=[]
    res_zone.append(Tprt)
    for ix in range(nx): #for each printed location
        res_zone.append(dataobs[ie][ix])
    
    #transform into Numpy array to facilitate exportation
    res_zone_array=np.array(res_zone)
    np.savetxt(path+"RESULT_FlowZone"+str(ie)+".txt",res_zone_array.transpose())

#Computation and export of average concentration in flow zones (considering a mixing of flow zones)
#Generate One unique file for the all simulation, containing in the first column the time-steps
# in each of the next colums, the average concentrations corresponding to the desired printed locations
aver_conc=[]
aver_conc.append(Tprt)
for ix in range(nx): #for each printed location
    sumQ = 0.
    sumQC = 0.
    for ie in range(ne): #for each zone
        qi = dataset[ie+1][-1][0] #input flow rate of the zone
        qici = qi*dataobs[ie][ix,:]
        sumQ += qi
        sumQC +=qici
    aver_conc.append((sumQC/sumQ))    
    
#transform into Numpy array to facilitate exportation
aver_conc_array=np.array(aver_conc)
np.savetxt(path+"RESULT_AverageFlow.txt",aver_conc_array.transpose())
    
