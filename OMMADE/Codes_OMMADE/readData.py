# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 10:27:49 2017

@author: Anne-Julie Tinet
"""

import numpy as np

from classDataPoint import *
from classParameters import *

def lire_ligne(pfile):
    """Reads the next line from pfile, skipping commentary lines and white lines."""
    
    while True:
        
        ligne = pfile.readline()
        
        if ligne.isspace():
            continue
        
        if ligne[0] == "#":
            continue
        
        return ligne
        
def readGeneralData(filename):
    """Reads simulation data (includes printing data) from a given file."""
    
    pfile = open(filename,'r')
    
    # Read space discretisation information
    dx,xmax = lire_ligne(pfile).split()
    dx = float(dx)
    nx = int(float(xmax)/dx)+1

    # Read time discretisation information
    dt,tmax = lire_ligne(pfile).split()
    dt = float(dt)
    tmax = float(tmax)
    
    # Read initial concentration information
    cini = lire_ligne(pfile).split()
    c0 = [float(c) for c in cini]

    # Read printing information
    xtype,ttype = lire_ligne(pfile).split()
    xtype = int(xtype)
    ttype = int(ttype)
    
    # Read printing location information
    
    # Print all locations
    if xtype == 0:
        Xprt = np.linspace(0,float(xmax),nx)
        
    else:
        # Number of printing locations
        nxprt = int(lire_ligne(pfile))
        
        # Equi-spaced printing locations
        if xtype == 1:
            Xprt = np.linspace(0,float(xmax),nxprt)
        
        # Determined printing locations
        else:
            Xprt = np.zeros(nxprt)
            for ip in range(nxprt):
                Xprt[ip] = float(lire_ligne(pfile))

    # Read printing time information
    
    # Print all times                
    if ttype == 0:
        Tprt = np.arange(0,tmax+dt,dt)
        
    else:
        # Number of printing times
        ntprt = int(lire_ligne(pfile))
        
        # Equi-spaced printing times
        if ttype == 1:
            Tprt = np.linspace(0,tmax,ntprt)
        
        # Determined printing times
        else:
            Tprt = np.zeros(ntprt)
            for ip in range(ntprt):
                Tprt[ip] = float(lire_ligne(pfile))
        
    pfile.close()
    
    return dx,nx,dt,tmax,c0,Xprt,Tprt
    


def readDataset(filename,dx,dt):
    """Generates the dataset from a given file."""
    
    pfile = open(filename,'r')
    
    # Number of flow types and reaches
    ne, nr = lire_ligne(pfile).split()
    ne = int(ne)
    nr = int(nr)
    
    # Initialisation
    dataset = [[]]
    for ie in range(ne):
        dataset.append([])
        
    # Flow rate for each flow type
    flow = [float(Q) for Q in lire_ligne(pfile).split()]
            
    # Allows dealing with evolutive flow rate in space
    qi = flow.copy()
    qf = flow.copy()
    
    # Initialise CFL condition (time step)
    cfl = [dt for Q in flow]
    
    # Read reach dependent parameters
    for ir in range(nr):
        
        # Reach start location and length
        x,lx = lire_ligne(pfile).split()
        dataset[0].append((float(x),float(lx)))
        
        
        # Physical parameters (area, dispersion, degradation rate, lateral flow rate, lateral concentration)
        for ie in range(ne):
            A, D, lam, ql, qlout, cl = lire_ligne(pfile).split()
            dataset[ie+1].append(Parameters(float(A), float(D), float(lam), float(ql), float(qlout), float(cl)))
            
            # Updating CFL condition
            if ir == 0:
                
                qf[ie] += float(lx)*(float(ql) - float(qlout))

                if flow[ie] == 0:
                    cfl[ie] = None
                else:
                    cfl[ie] = min(dx*float(A)/qi[ie], dx*float(A)/qf[ie],cfl[ie])

            else:
                
                qi[ie] = qf[ie]
                qf[ie] += float(lx)*(float(ql) - float(qlout))
                
                if cfl[ie] != None:
                    cfl[ie] = min(cfl[ie],dx*float(A)/float(qi[ie]),dx*float(A)/float(qf[ie]))
                                    
        # Physical parameters (exchange rates)
        for ie in range(ne):
            ligne = lire_ligne(pfile).split()
            alpha = {}

            for je in range(ne):
                if je != ie:
                    alpha[je] = float(ligne[je])
                    
            dataset[ie+1][ir].setAlpha(alpha)
    
    # Flow related parameters (for advection)
    for ie in range(ne):
        if cfl[ie] != None:
            dataset[ie+1].append((float(flow[ie]),0.99999999999999*cfl[ie]))
        else:
            dataset[ie+1].append((float(flow[ie]),None))
        
    pfile.close()
    return dataset
    
    
def readBound(filename):
    """Reads boundary data (Concentrations at x = 0)."""
    
    pfile = open(filename,'r')
    
    # Number of boundary points
    nb = int(lire_ligne(pfile))
    
    # Initialisation
    bound = [[],[]]
    
    # Reading boundary data for the first point
    data = lire_ligne(pfile).split()

    tb = float(data[0])
    cb = [float(c) for c in data[1:]]

    bound[0].append(tb)
    
    for ie in range(len(cb)):
        bound[1].append([cb[ie]])
        
    # Reading boundary data for each point
    for ib in range(1,nb):
        
        data = lire_ligne(pfile).split()

        tb = float(data[0])
        cb = [float(c) for c in data[1:]]

        bound[0].append(tb)
    
        for ie in range(len(cb)):
            bound[1][ie].append(cb[ie]) 
            
    pfile.close()
            
    return bound
    
    
    
