# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 17:11:44 2016

1D Stream Channel Profile Model

@author: dav
"""

import numpy as np
import matplotlib as mpl

streamlength = 3000
numnodes = 300
movern = 0.85
n = 1
U = 0.01
K = 5.0E-07
k_a = 6.69  #Hack's constant
h = 1.667  #Hack's law exponent
chanhead_pos = 500 #Channel head fixed position
Scrit_deg = 30 #Critical slope for failure in degrees


def calc_Scrit_radians(Scrit_deg):
    return np.tan(np.radians(Scrit_deg))

def calc_DX(nodes,streamlength):
    return streamlength/numnodes

def calc_m(movern,n):
    return movern*n

# Calculate the simple steady state profile
# Using the simple stream power law
def SimpleSteadyState():
    x_coords = np.linspace(0.0, streamlength, num=numnodes+1)
    z_coords_raw = np.zeros(numnodes+1)
    
    print len(x_coords), len(z_coords_raw)
    
    m = calc_m(movern, n)
    
    print x_coords
    
    Scrit_radians = calc_Scrit_radians(Scrit_deg)
    print Scrit_radians    
    
    DX = calc_DX(numnodes, streamlength)
    print DX
    
    for i in range(len(x_coords)-1, 0, -1):
        if x_coords[i] < chanhead_pos and i != numnodes:
            z_coords_raw[i] = z_coords_raw[i+1] + Scrit_radians*DX
        else:
            z_coords_raw[i] = -( (U / (K*k_a**m) )**(1/n)) * (x_coords[i]**(1-h*m/n)) / (1-h*m/n)
            
    #z_coords_raw[z_coords_raw < channhead_pos]         
            
    # Get the last element of the z coordinate array        
    z_raw_end = z_coords_raw[-1]     
    
    z_coords = z_coords_raw - z_raw_end
    
    print z_coords
            

SimpleSteadyState()    
