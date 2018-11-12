# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 08:55:44 2017

@author: Anne-Julie Tinet

This module contains the time loops fonction. Advection loop contains the iteration
to perform advection for one time step (the advection time step may be shorter 
than the overall time step due to the CFL condition).

The time loop does the overall loop for each time step, until reaching the total
simulation time. It returns a matrix containing the concentration for each printing
location at each printing time step
"""

import numpy as np
import matplotlib.pyplot as plt

from classParameters import Parameters
from classDataPoint import DataPoint
from readData import readDataset


def initialise(filedat, dt, dx, c0, nx, scheme):
    """Initialisation for the simulation.
    Returns the dataset, the point information and the initial concentration."""
    
    dataset = readDataset(filedat, dx, dt, scheme)
    ne = len(dataset)
    ne = ne-1
    
    C = np.ones(nx*ne)
    for ie in range(ne):
        for ix in range(nx):
            C[ix + ie*nx] = c0[ie]
        
    
    points = []
    
    for ie in range(ne):
        for ix in range(nx):
            
            points.append(DataPoint(dataset, ix, ie, dt, dx))
    
    return dataset, C, points


def advectionCFLloop(points, C, data, nx, ne, clim, dt):
    
    """Perform advection for one overall time step. Modifies the concentration
    vector."""
    
    # Initialisation
    Cprime = np.zeros(nx)
    nr = len(data[0])  # Number of reaches
    
    # For each flow type
    for ie in range(ne):
        
        # Initialisation of advection time loop
        t = 0

        # Skipping storage zones
        if (data[ie+1][nr][1] != None):
            
            # Loop until reaching the overall time step
            while (t<0.999999999*dt):
                
                # Advection time step 
                if t + data[ie+1][nr][1] > dt:
                    corr = (dt - t)/data[ie+1][nr][1]
                else:
                    corr = 1
                    
                t += data[ie+1][nr][1]
            
                # Performs advection 
                for ix in range(nx):
                    Cprime[ix] = points[ix + ie*nx].advectionPoint_explicite(C, nx, clim, corr)
                
                # Updates concentration data
                C[ie*nx:(ie+1)*nx] = Cprime[:]


def timeloop(points, C, data, nx, bound, dx, dt, tmax, Xprt, Tprt, scheme):
    
    """Performs the simulation. Returns the concentration for each flow types at 
    each concentration printing time and printing step.
    
    Returns a list under the form : Data[flow type][printing location, printing time]."""
    
    # Initialisation
    nt = int(tmax/dt)+1
    ne = len(data)-1

    # Xprt corresponds to physical location and Tprt to time and Tiprt to time indexes
    Tiprt = [int(t/dt) for t in Tprt]

    dataobs = []
    for ie in range(ne):
        dataobs.append(np.zeros((min(nx,len(Xprt)),min(nt,len(Tprt)))))
    
    A = np.zeros((nx*ne,nx*ne))
    B = np.zeros(nx*ne)
    BC = np.zeros(nx*ne)
    clim = np.zeros(ne)
    
    I = np.eye(nx*ne)
    
    # Filling the system matrix and second member for implicit scheme
    for point in points:
            
        # Always start with dispersion
        point.dispersionPoint(A, nx)
        point.massloss(A, B, nx)
        point.massexchange(A, nx)
        
        if scheme == 1:
            point.advectionPoint_cranknicholson(A, nx)
            AE = I - A
    
    # Time loop
    for t in range(nt):

        # Updating the system second member with boundary condition
        for ie in range(ne):
            clim[ie] = np.interp(t*dt,bound[0],bound[1][ie])
            BC[ie*nx] = points[ie*nx].Dm_*clim[ie]

            if scheme == 1:
                BC[ie*nx] += points[ie*nx].U_*clim[ie]/2

        if scheme == 0:
            # Perform advection
            advectionCFLloop(points, C, data, nx, ne, clim, dt)
            
            # Perform dispersion, mass exchange and mass loss
            C = np.linalg.solve(A,B+BC+C)
            
        else:
            
            # Perform advection, dispersion, mass exchange and mass loss
            C = np.linalg.solve(I - 0.5*AE, B + BC + np.dot(I + 0.5*AE, C))
        
        # Saving data if corresponding to  a printing time
        # The printing time need to correspond to a simulation time (no time interpolation)
        if t in Tiprt:           
            
            it = Tiprt.index(t)
            
            for ie in range(ne):
                
                i0 = 0
                
                # Saving data if corresponding to  a printing location
                # Data may be interpolated (linear interpolation) if necessary
                
                for ixp in range(0,len(Xprt)):
                    
                    if i0 == nx-1:
                        dataobs[ie][ixp,it] = C[i0 + ie*nx]
                        continue
                    
                    for ix in range(i0+1,nx):
                        check = False
                        
                        if (ix-1)*dx < Xprt[ixp] <= ix*dx:
                            
                            dataobs[ie][ixp,it] = ((ix*dx-Xprt[ixp])*C[ix-1 + ie*nx]+(Xprt[ixp]-(ix-1)*dx)*C[ix + ie*nx])/dx
                            check = True

                        elif ix*dx < Xprt[ixp]:
                            i0 += 1
                            
                        else:
                            dataobs[ie][ixp,it] = C[ix-1 + ie*nx]
                            check = True
                            
                        if check:
                            break
                    
                    

    return dataobs