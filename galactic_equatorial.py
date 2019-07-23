import numpy as np

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
    
    '''
    transfroms velocities from equatorial (mu_ra_star, mu_dec) to galactic (mu_l_star, mu_b)
    mulstar and murastar imply multiplication by Cos[b] and Cos[Dec], respectively 
    '''
        
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
    
    
def LB_corr(l, b, dist, U, V, W, R0, VR0, VR):
    
    '''
    Calculates corrections to proper motion (in Galactic coordinates)
    for differential Galactic rotation (DGR)
    and Sun's peculiar motion (LSR)
    See e.g., Verbunt F., Igoshev, A., & Cator E., 2017, A&A, 608, 57.
    
    output numpy[[mu_l_star_DGR, mu_b_DGR], [mu_l_star_LSR, mu_b_LSR]]  
    mulstar implies multiplication by Cos[b]  
    '''
    
    [l, b] = np.radians([l, b])
    
    thl = np.arctan2(R0 * np.sin(l), R0 * np.cos(l) - dist * np.cos(b))
    
    mulstarDGR = -VR0 * np.cos(l) + VR * np.cos(thl)
    mubDGR = (VR0 * np.sin(l) - VR * np.sin(thl)) * np.sin(b)
    mulstarLSR = U * np.sin(l) - V * np.cos(l)
    mubLSR = U * np.cos(l) * np.sin(b) + V * np.sin(l) * np.sin(b) - W * np.cos(b)
    
    return np.array([[mulstarDGR, mubDGR], [mulstarLSR, mubLSR]]) / (4.74 * dist)
    
    
def RaDec_corr(ra, dec, dist, U, V, W, R0, VR0, VR): 
    
    '''
    Calculates corrections to proper motion (in equatorial coordinates)
    for differential Galactic rotation (DGR)
    and Sun's peculiar motion (LSR)
    See e.g., Verbunt F., Igoshev, A., & Cator E., 2017, A&A, 608, 57.
    
    output numpy[[mu_ra_star_DGR, mu_dec_DGR], [mu_ra_star_LSR, mu_dec_LSR]]  
    mu_ra_star implies multiplication by Cos[b]  
    '''
      
    l, b = RaDec2LB(ra, dec)
    lbcorr = LB_corr(l, b, dist, U, V, W, R0, VR0, VR);
    
    return np.array([muLB2RaDec(l, b, *lbcorr[0]), muLB2RaDec(l, b, *lbcorr[1])])
    
    
def hms2deg(h, m, s):
    'transforms hours, min, sec into degrees'
    return (h + m / 60 + s / 3600) * 360 / 24
    

def dms2deg(d, m, s):
    'transforms degrees, min, sec into degrees'
    return d + np.sign(d) * (m / 60 + s / 3600)    