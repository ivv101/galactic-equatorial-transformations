import numpy as np
import inspect

ra0 = 192.85948 # NGP
dec0 = 27.12825 # NGP
l0 = 32.93192 + 90 # NCP

[ra0, dec0, l0] = np.radians([ra0, dec0, l0])


def RaDec2LB(ra, dec):
    
    '''
    transfroms coordinates from equatorial (ra, dec) to galactic (l, b)
    input/output in degrees
    '''
            
    [ra, dec] = np.radians([ra, dec])
    
    b = np.arcsin(np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0))
        
    sn = (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / np.cos(b)
    cs = np.cos(dec) * np.sin(ra - ra0) / np.cos(b) 
    
    l = np.mod(l0 - np.arctan2(cs, sn), 2 * np.pi)
    
    return np.degrees([l, b])
    
    
def LB2RaDec(l, b):
    
    '''
    transfroms coordinates from galactic (l, b) to equatorial (ra, dec)
    input/output in degrees
    '''
    
    [l, b] = np.radians([l, b])
    
    dec = np.arcsin(np.sin(dec0) * np.sin(b) + np.cos(dec0) * np.cos(b) * np.cos(l0 - l))
        
    sn = (np.cos(dec0) * np.sin(b) - np.sin(dec0) * np.cos(b) * np.cos(l0 - l)) / np.cos(dec)
    cs = np.cos(b) * np.sin(l0 - l) / np.cos(dec) 
    
    ra = np.mod(ra0 + np.arctan2(cs, sn), 2 * np.pi)
    
    return np.degrees([ra, dec]) 
    
    
def muRaDec2LB(ra, dec, murastar, mudec):
    
    'transfroms velocities from equatorial (mu_ra_*, mu_dec) to galactic (mu_l_*, mu_b)'
        
    [ra, dec] = np.radians([ra, dec])
        
    sn = np.cos(dec) * np.sin(dec0) - np.cos(dec0) * np.cos(ra - ra0) * np.sin(dec)
    cs = np.cos(dec0) * np.sin(ra - ra0)
        
    phi = np.arctan2(cs, sn)
    
    rot_mat = np.array([[np.cos(phi), np.sin(phi)], [-np.sin(phi), np.cos(phi)]])
    
    return rot_mat.dot([murastar, mudec]) 
    

def muLB2RaDec(l, b, mulstar, mub):
    
    'transfroms velocities from galactic (mu_l_*, mu_b) to equatorial (mu_ra_*, mu_dec)'
    
    [l, b] = np.radians([l, b])
        
    sn = np.cos(b) * np.sin(dec0) - np.cos(dec0) * np.cos(l0 - l) * np.sin(b)
    cs = -np.cos(dec0) * np.sin(l0 - l)
        
    phi = np.arctan2(cs, sn)
    
    rot_mat = np.array([[np.cos(phi), np.sin(phi)], [-np.sin(phi), np.cos(phi)]])
    
    return rot_mat.dot([mulstar, mub])    
    
    
def DGRandLSRcorrections(l, b, dist, U, V, W, R0, VR0, VR):
    
    '''
    output numpy[mulstarDGR, mubDGR, mulstarLSR, mubLSR]
    '''
    
    [l, b] = np.radians([l, b])
    
    thl = np.arctan2(R0 * np.sin(l), R0 * np.cos(l) - dist * np.cos(b))
    
    mulstarDGR = -VR0 * np.cos(l) + VR * np.cos(thl)
    mubDGR = (VR0 * np.sin(l) - VR * np.sin(thl)) * np.sin(b)
    mulstarLSR = U * np.sin(l) - V * np.cos(l)
    mubLSR = U * np.cos(l) * np.sin(b) + V * np.sin(l) * np.sin(b) - W * np.cos(b)
    
    return np.array([mulstarDGR, mubDGR, mulstarLSR, mubLSR]) / (4.74 * dist)
    
      
# TEST
if __name__ == "__main__":
    
    #import sys
    #print(muRaDec2LB(float(sys.argv[1]), float(sys.argv[2]), 10, 10))
    
    print(inspect.cleandoc('''
    Vela pulsar test
    Corrections to the observed Vela proper motion:

    In galactic coordinates
    Observed proper motion: -53.6, -21.9
    Correction for the solar motion: 6.7, 4.9
    Correction for the galactic rotation: 5.35, 0

    Converted back to RA and Dec: -38.5, 23.1

    output:
    '''))

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

    corr = DGRandLSRcorrections(l, b, dist, *pars).reshape(2,2)
    corr = [-53.6, -21.9] - np.array(sum(corr))
    print(muLB2RaDec(l, b, *corr))