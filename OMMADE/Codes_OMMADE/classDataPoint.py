# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 11:01:19 2017

@author: Anne-Julie Tinet
"""

import numpy as np
from classParameters import *

class DataPoint:
    """Contains the physical data required for computation for one point 
    (flow type and location) of the domain). 
    Allows the calculation of advection (explicit, Lax-Wendroff scheme),
    dispersion (implicit scheme), mass exchange between two flow types point at
    the same location and mass exchange through either degradation on lateral
    flow."""
    
    def __init__(self, dataset, ix, ie, dt, dx):
        """Creates the point from the general information."""
        
        self.ix_ = ix   # Point location (index)
        self.ie_ = ie   # Zone index
        
        # Generation of the physical parameters from the reach number and dataset
        # All parameters are generated in adimensional form
        
        # Retrieving reach information
        nr = len(dataset[0])
        ir = self.reachNumber(dataset,ix,dx)
        
        # Adimensional flow (constraint by CFL condition)
        area = dataset[ie+1][ir].getArea()
        dtcfl = dataset[ie+1][nr][1]
        
        if dtcfl != None:
            dtcfl = min(dataset[ie+1][nr][1],dt)
            self.U_ = dataset[ie+1][nr][0]*dtcfl/dx/area 
            
        else:
            self.U_ = None
            
        # dispersivity x-dx/2, x+dx/2 - adimensional
        self.Dm_, self.Dp_ = self.initDispersivity(dataset, ix, ie, ir, area, dt, dx) 
        
        # exchange rates - adimensional 
        self.alpha_ = {}
        alp = dataset[ie+1][ir].getAlpha()
        for flowId in alp:

            self.alpha_[flowId] = alp[flowId]*dt/area
        
        # degradation rate - adimensional
        self.lambda_ = dataset[ie+1][ir].getLambda()*dt
        
        # Lateral flow (adimensional) and concentration
        self.ql_, self.cl_ = dataset[ie+1][ir].getLateral()
        self.ql_ *= (dt/area)
        
        
    def reachNumber(self, dataset, ix, dx):
        """Calculates and returns the reach index from the location."""
        
        for ir in range(len(dataset[0])):
            
            x,lx = dataset[0][ir]
            if x <= ix*dx + dx/2 < x+lx:

                return ir
                
        return len(dataset[0])-1
        
    def initDispersivity(self, dataset, ix, ie, ir, area, dt, dx):
        """Calculates and returns the dispersivity at faces x-dx/2 (Dm) and
        x+dx/2 (Dp)."""
        
        nr = len(dataset[0])
        
        # Calculation for dispersive transport only
        if (dataset[ie+1][ir].getD() != 0):
            
            # Initialisation at the point (x) value
            Di = dataset[ie+1][ir].getD()*dt/(dx**2)
            Dp = Dm = Di

            x,lx = dataset[0][ir]
            
            # Weighed average if the points connects two different reaches
            if ((ix+1)*dx >= x+lx) and (ir < nr-1):
                Dp = dataset[ie+1][ir+1].getD()*dt*dataset[ie+1][ir+1].getArea()/(dx**2)/(2*area)
                Dp += Di/2
                
            if (ix*dx < x) and (ir > 0):
                Dm = dataset[ie+1][ir-1].getD()*dt*dataset[ie+1][ir-1].getArea()/(dx**2)/(2*area)
                Dm += Di/2

            return Dm, Dp
            
        return 0., 0.
        
    def advectionPoint(self, C, nx, cin, corr):
        """Performs one (advection) time step using explicit Lax-Wendroff scheme.
        The Lax-Wendroff scheme is second order accurate. 
        Returns and modifies the concentration."""
        
        if corr < 1:
            q = self.U_*corr  # Adimensional flow
        else:
            q = self.U_  # Adimensional flow
            
        i = self.ix_ + self.ie_*nx  # Global index in concentration vector
        clim = cin[self.ie_]        # Boundary condition at x = 0
        
        # If non-advective transport - do nothing
        if q != None:
            
            if self.ix_ == 0:
                return q*(q/2 - 1/2)*C[i+1] + (1-q**2)*C[i] + q*(q/2 + 1/2)*clim
    
            elif self.ix_ == nx-1:
                
                return (1-q)*C[i] + q*C[i-1]
    
            else:
                
                return q*(q/2 - 1/2)*C[i+1] + (1-q**2)*C[i] + q*(q/2 + 1/2)*C[i-1]
    
        return C[self.ix_]


    def dispersionPoint(self, A, nx):
        """Calculates and modifies the system matrix A for implicit scheme considering dispersion. 
        The second member is 0 and therefore not modified. The Laplacian discretisation scheme is second order."""
        
        dp,dm = self.Dp_, self.Dm_
        i = self.ix_ + self.ie_*nx # Global index in concentration vector
        
        if self.ix_ < nx-1:
            A[i,i] = 1+dp+dm
            A[i,i+1] = -dp
    
            if self.ix_ != 0:

                A[i,i-1] = -dm
    
        else:
            A[i,i] = 1+dm
            A[i,i-1] = -dm


    def massloss(self, A, B, nx):
        """Calculates and modifies the system matrix A and second member for implicit scheme 
        considering mass loss from degradation, sedimentation and lateral source/sink."""
        
        deg = self.lambda_
        ql,cl = self.ql_, self.cl_
        i = self.ix_ + self.ie_*nx
        
        A[i,i] += deg + ql
        B[i] += ql*cl


    def massexchange(self, A, nx):
        """Calculates and modifies the system matrix A and second member for implicit scheme 
        considering mass exchange with other zones."""
        
        alp = self.alpha_
        i = self.ix_ + self.ie_*nx
        
        for flowId in alp:
            
            A[i,i] += alp[flowId]
            A[i,self.ix_+flowId*nx] = -alp[flowId]
        

