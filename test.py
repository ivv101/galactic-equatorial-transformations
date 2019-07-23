import numpy as np
import galactic_equatorial as ge

print(
'''
Vela pulsar test 
Corrections to the observed Vela proper motion:
    
In galactic coordinates
Observed proper motion: -53.6, -21.9
Correction for the solar motion: 6.7, 4.9
Correction for the galactic rotation: 5.35, 0
    
Converted back to RA and Dec: -38.5, 23.1
    
output:
''')

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
            
l = 263.55196463 # degrees
b = -2.7872552 # degrees 
dist = 0.293 # kpc
    
corr = ge.DGRandLSRcorrections(l, b, dist, *pars).reshape(2,2)
corr = [-53.6, -21.9] - np.array(sum(corr)) 
print(ge.muLB2RaDec(l, b, *corr))

print(ge.DGRandLSRcorrections.__doc__, ge.muLB2RaDec.__doc__)




















   