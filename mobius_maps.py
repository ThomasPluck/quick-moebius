#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  4 21:10:50 2018

@author: thomas
"""

import numpy as np
from repetend import no_repet, look, check

#Takes arguements of 4-tuples parameterizing
#mobius transforms and performs calculations
#for easy application

class MoebiusTransform:
    
    def __init__(self,param):
        #Allow input of parameter tuple or numpy array
        if type(param) == tuple:
            a,b,c,d = param[0],param[1],param[2],param[3]
        
            self.array = np.array([[a,b],
                                   [c,d]])
        elif type(param) == np.ndarray:
            self.array = param
            a,b,c,d = param[0,0],param[0,1],param[1,0],param[1,1]
        #Calculate fixed points
        if c!= 0+0j:
            self.fixed_points = [((a-d)+((a-d)**2+(4*b*c))**0.5)/(2*c),
                                 ((a-d)-((a-d)**2+(4*b*c))**0.5)/(2*c)]
        else:
            if a != d:
                self.fixed_points = [float('inf'),-b/(a-d)]
            elif b == 0+0j:
                self.fixed_points = [float('inf'),0]
            else:
                self.fixed_points = [float('inf'),float('inf')]

        #Classify points as either sinks or sources
        self.sink,self.source = None,None

        #Test stationary points for either being a sink or a source
        delt = 1e-11
        n = 24
        ref = [np.exp(0+1j*i*np.pi/n) for i in range(n)]
        eps = []
        for p in self.fixed_points:
            if p != float('inf'):
                for r in ref:
                    eps.append(abs(((self.array[0,0]*(p+delt*r)+self.array[0,1])/\
                               (self.array[1,0]*(p+delt*r)+self.array[1,1]))-p))
            
    
            if all(e>delt for e in eps):
                if self.sink == None:
                    self.sink = p
            elif all(e<delt for e in eps):
                if self.sink == None:
                    self.source = p
            eps = []
        
        if self.sink != None and self.source == None:
            try:
                self.source = [s for s in self.fixed_points if s != self.sink][0]
            except:
                self.source = self.sink
        if self.source != None and self.sink == None:
            try:
                self.sink = [s for s in self.fixed_points if s != self.source][0]
            except:
                self.sink = self.source
        
        #Classify transform
        
        trace = np.trace(self.array)**2
        w,v = np.linalg.eig(self.array)
        
        if trace == 0:
            self.kind = "Circular"
        elif trace < 4 and trace >= 0 and w[0]**2 == -1:
            self.kind = "Elliptic"
        elif trace == 4 and w[0] == 1:
            self.kind = "Parabolic"
        elif trace > 4 and w[0].real > 0 and not np.iscomplex(w[0]):
            self.kind = "Hyperbolic"
        elif trace < 0 or trace > 4 and abs(w[0]) != 1:
            self.kind = "Loxodromic"
        else:
            self.kind = "Unclassified"
          
    #Transform a complex number with the transform        
    
    def transform(self,num):
        return (self.array[0,0]*num+self.array[0,1])/(self.array[1,0]*num+self.array[1,1])

#Semigroup object which takes a single arguement of a list of
#Moebius transforms. Onboard applications include limitset
#calculation and orbit calculation

class SemiGroup:
    
    #Index tuples convert into numpy
    def __init__(self,gen_list):
        self.gen_list = {str(idx):np.array([v[0],v[1]],[v[2],v[3]]) for idx,v in enumerate(gen_list)}
        self.limit_set = None
        self.orbits = None
    
    def compute_limits(self,depth,breadth):
        for d in depth:
            for k,v in self.gen_list.iteritems():
                for i in range(max(float(self.gen_list.keys()))):
                    if no_repet(i+k):
                        if check(self.gen_list.keys(),i+k):
                            self.gen_list.update({i+k,np.matmul(self.gen_list[i],
                                                                self.gen_list[k])})
        for b in breadth:
            return False
    def compute_orbits(grades):
        return False