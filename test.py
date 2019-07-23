import numpy as np
import galactic_equatorial as ge

print(
'''
Test using Vela pulsar proper motion
(Dodson, R., Legge, D., Reynolds, J. E., et al. 2003, ApJ, 596, 1137):
'''

# In galactic coordinates
# Observed proper motion: -53.6, -21.9
# Correction for the solar motion: 6.7, 4.9
# Correction for the galactic rotation: 5.35, 0
#
# Converted back to RA and Dec: -38.5, 23.1
#
# 
# ''')

# Dehnen & Binney 1998:

U = 10.0 # km/s
V = 5.25 # km/s
W = 7.17 # km/s

# Fich 1989:

R0 = 8.5 # kpc
VR0 = 220 # km/s
VR = 220 # km/s

# # Schonrich,Binney& Dehnen (2010):
# U = 11.1 # km/s
# V = 12.24 # km/s
# W = 7.25 # km/s
#
# # Reid 1993; Carrillo et al. 2018
# R0 = 8.0 # kpc
# VR0 = 230 # km/s
# VR = 230 # km/s

pars = [U, V, W, R0, VR0, VR]

def hms2deg(h, m, s):
    return (h + m / 60 + s / 3600) * 360 / 24

def dms2deg(d, m, s):
    return d + np.sign(d) * (m / 60 + s / 3600)
    
    
ra = [8, 35, 20.55] 
dec = [-45, 10, 34.8]  

ra, dec = hms2deg(*ra), dms2deg(*dec)

print("radec: ", ra, dec) 

l, b = ge.RaDec2LB(ra, dec)
print("lb: ", l, b)

corr = ge.RaDec_corr(ra, dec, 0.287, U, V, W, R0, VR0, VR)

muradec = [-49.68, 29.9]

print("corrections, radec: ", corr)
print("proper motion, radec:", muradec)
print("corrected proper motion, radec: ", muradec - sum(corr))














   