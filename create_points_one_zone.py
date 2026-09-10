import numpy as np
from scipy.stats import norm
import input_params as ip
import matplotlib.pyplot as plt
import fun_time_independent_mh as fi

p_arr_core_min = np.zeros(ip.m * (2*ip.n_kin-2)+1)
gamma_arr_core_min = np.zeros(ip.m * (2*ip.n_kin-2)+1)

psi_arr_core_min = np.zeros( ip.m * (2*ip.n_kin-2)+1)
alpha_arr_core_min = np.zeros(ip.m * (2*ip.n_kin-2)+1)

p_arr_core_core_1 = np.zeros(ip.m * (2*ip.n_kin-2)+1)
gamma_arr_core_core_1 = np.zeros(ip.m * (2*ip.n_kin-2)+1)

psi_arr_core_core_1 = np.zeros( ip.m * (2*ip.n_kin-2)+1)
alpha_arr_core_core_1 = np.zeros(ip.m * (2*ip.n_kin-2)+1)

p_arr_pfr_max = np.zeros(ip.m * (2*ip.n_kin-2)+1)
gamma_arr_pfr_max = np.zeros(ip.m * (2*ip.n_kin-2)+1)

psi_arr_pfr_max = np.zeros(ip.m * (2*ip.n_kin-2)+1)
alpha_arr_pfr_max = np.zeros(ip.m * (2*ip.n_kin-2)+1)

p_arr_pfr_pfr_1  = np.zeros(ip.m * (2*ip.n_kin-2)+1)
gamma_arr_pfr_pfr_1 = np.zeros(ip.m * (2*ip.n_kin-2)+1)

psi_arr_pfr_pfr_1 = np.zeros(ip.m * (2*ip.n_kin-2)+1)
alpha_arr_pfr_pfr_1 = np.zeros(ip.m * (2*ip.n_kin-2)+1)

for m_period in range(ip.m):
    g = -ip.A
    psi_0 = (ip.iota_res - ip.iota_b) / ip.iota_a
    
#==============================================================================
# Min point  
#==============================================================================
   
    p_min = np.sqrt(ip.iota_a) * (ip.psi_t_min-psi_0) 
    E_tilde_min = p_min ** 2 / 2  
    
    kin_tilde_arr_min = np.zeros(ip.n_kin)
    p_arr_min = np.zeros(ip.n_kin)
    p_arr_neg_min = np.zeros(ip.n_kin)
    
    gamma_arr_min = np.zeros(ip.n_kin)
    gamma_arr_neg_min = np.zeros(ip.n_kin) 
    
    p_arr_tot_min = np.zeros(2*ip.n_kin-2)
    gamma_arr_tot_min = np.zeros(2*ip.n_kin-2)
    
    kin_tilde_arr_min = fi.poloidal_spacing(E_tilde_min-2*g, 
                                            E_tilde_min, ip.n_kin,
                                            ip.alpha_pol)
    
    for i_kin_tilde in range(len(kin_tilde_arr_min)):
        gamma_arr_min[i_kin_tilde] = \
            np.arccos(- (E_tilde_min - kin_tilde_arr_min[i_kin_tilde]) / g + 1) \
            + 2 * np.pi * m_period

        p_arr_neg_min[i_kin_tilde] = - \
            np.sqrt(2*kin_tilde_arr_min[i_kin_tilde])

        gamma_arr_neg_min[i_kin_tilde] = \
            -np.arccos(- (E_tilde_min - kin_tilde_arr_min[ i_kin_tilde]) / g + 1) \
            + 2 * np.pi * m_period
    
    # from -pi to 0
    p_arr_tot_min[0:ip.n_kin-1] = p_arr_neg_min[0:ip.n_kin-1]
    gamma_arr_tot_min[0:ip.n_kin-1] = gamma_arr_neg_min[0:ip.n_kin-1]
    
    # from 0 to pi
    p_arr_tot_min[ip.n_kin-1:2*ip.n_kin-2] = p_arr_neg_min[ip.n_kin-1:0:-1]
    gamma_arr_tot_min[ip.n_kin-1:2*ip.n_kin-2] = gamma_arr_min[ip.n_kin-1:0:-1]
    
    
    p_arr_core_min[(2*ip.n_kin-2)*m_period:(2*ip.n_kin-2)*(m_period+1)] = p_arr_tot_min[0:2*ip.n_kin-2]
    gamma_arr_core_min[(2*ip.n_kin-2)*m_period:(2*ip.n_kin-2)*(m_period+1)] = gamma_arr_tot_min[0:2*ip.n_kin-2]
    
    if m_period == ip.m-1:
        p_arr_core_min[(2*ip.n_kin-2)*(m_period+1)] = p_arr_neg_min[0]
        gamma_arr_core_min[(2*ip.n_kin-2)*(m_period+1)] = gamma_arr_min[0]
        
    psi_arr_core_min = fi.psi_from_momentum(p_arr_core_min, ip.iota_a, psi_0)
    alpha_arr_core_min = fi.alpha_from_gamma_mod(gamma_arr_core_min, ip.m)
    
#==============================================================================
# Core 1 point     
#==============================================================================   
 
    E_tilde_core_1 = 2*g + ip.frac_d_1 * 2*g 

    kin_tilde_arr_core_1 = np.zeros(ip.n_kin)
    p_arr_core_1 = np.zeros(ip.n_kin)
    p_arr_neg_core_1 = np.zeros(ip.n_kin)
    
    gamma_arr_core_1 = np.zeros(ip.n_kin)
    gamma_arr_neg_core_1 = np.zeros(ip.n_kin) 
    
    p_arr_tot_core_1 = np.zeros(2*ip.n_kin-2)
    gamma_arr_tot_core_1 = np.zeros(2*ip.n_kin-2)
    
    kin_tilde_arr_core_1 = fi.poloidal_spacing(E_tilde_core_1-2*g, 
                                               E_tilde_core_1, ip.n_kin,
                                               ip.alpha_pol)    
    
    for i_kin_tilde in range(len(kin_tilde_arr_core_1)):
        gamma_arr_core_1[i_kin_tilde] = \
            np.arccos(- (E_tilde_core_1 - kin_tilde_arr_core_1[i_kin_tilde]) / g + 1) \
            + 2 * np.pi * m_period

        p_arr_neg_core_1[i_kin_tilde] = - \
            np.sqrt(2*kin_tilde_arr_core_1[i_kin_tilde])

        gamma_arr_neg_core_1[i_kin_tilde] = \
            -np.arccos(- (E_tilde_core_1 - kin_tilde_arr_core_1[ i_kin_tilde]) / g + 1) \
            + 2 * np.pi * m_period
    
    # from -pi to 0
    p_arr_tot_core_1[0:ip.n_kin-1] = p_arr_neg_core_1[0:ip.n_kin-1]
    gamma_arr_tot_core_1[0:ip.n_kin-1] = gamma_arr_neg_core_1[0:ip.n_kin-1]
    
    # from 0 to pi
    p_arr_tot_core_1[ip.n_kin-1:2*ip.n_kin-2] = p_arr_neg_core_1[ip.n_kin-1:0:-1]
    gamma_arr_tot_core_1[ip.n_kin-1:2*ip.n_kin-2] = gamma_arr_core_1[ip.n_kin-1:0:-1]
    
    
    p_arr_core_core_1[(2*ip.n_kin-2)*m_period:(2*ip.n_kin-2)*(m_period+1)] = p_arr_tot_core_1[0:2*ip.n_kin-2]
    gamma_arr_core_core_1[(2*ip.n_kin-2)*m_period:(2*ip.n_kin-2)*(m_period+1)] = gamma_arr_tot_core_1[0:2*ip.n_kin-2]
    
    if m_period == ip.m-1:
        p_arr_core_core_1[(2*ip.n_kin-2)*(m_period+1)] = p_arr_neg_core_1[0]
        gamma_arr_core_core_1[(2*ip.n_kin-2)*(m_period+1)] = gamma_arr_core_1[0]
        
    psi_arr_core_core_1 = fi.psi_from_momentum(p_arr_core_core_1, ip.iota_a, psi_0)
    alpha_arr_core_core_1 = fi.alpha_from_gamma_mod(gamma_arr_core_core_1, ip.m)

#==============================================================================
# Max point     
#============================================================================== 
    
    p_max = np.sqrt(ip.iota_a) * (ip.psi_t_max-psi_0) 
    E_tilde_max = p_max ** 2 / 2  
    
    kin_tilde_arr_max = np.zeros(ip.n_kin)
    p_arr_max = np.zeros(ip.n_kin)
    p_arr_neg_max = np.zeros(ip.n_kin)
    
    gamma_arr_max = np.zeros(ip.n_kin)
    gamma_arr_neg_max = np.zeros(ip.n_kin) 
    
    p_arr_tot_max = np.zeros(2*ip.n_kin-2)
    gamma_arr_tot_max = np.zeros(2*ip.n_kin-2)
    
    kin_tilde_arr_max = fi.poloidal_spacing(E_tilde_max-2*g, 
                                            E_tilde_max, ip.n_kin,
                                            ip.alpha_pol)
    
    for i_kin_tilde in range(len(kin_tilde_arr_max)):
        p_arr_max[i_kin_tilde] = np.sqrt(2*kin_tilde_arr_max[i_kin_tilde])
        gamma_arr_max[i_kin_tilde] = \
            np.arccos(- (E_tilde_max - kin_tilde_arr_max[i_kin_tilde]) / g + 1) \
            + 2 * np.pi * m_period

        gamma_arr_neg_max[i_kin_tilde] = \
            -np.arccos(- (E_tilde_max - kin_tilde_arr_max[ i_kin_tilde]) / g + 1) \
            + 2 * np.pi * m_period
    
    # from -pi to 0
    p_arr_tot_max[0:ip.n_kin-1] = p_arr_max[0:ip.n_kin-1]
    gamma_arr_tot_max[0:ip.n_kin-1] = gamma_arr_neg_max[0:ip.n_kin-1]
    
    # from 0 to pi
    p_arr_tot_max[ip.n_kin-1:2*ip.n_kin-2] = p_arr_max[ip.n_kin-1:0:-1]
    gamma_arr_tot_max[ip.n_kin-1:2*ip.n_kin-2] = gamma_arr_max[ip.n_kin-1:0:-1]
    
    
    p_arr_pfr_max[(2*ip.n_kin-2)*m_period:(2*ip.n_kin-2)*(m_period+1)] = p_arr_tot_max[0:2*ip.n_kin-2]
    gamma_arr_pfr_max[(2*ip.n_kin-2)*m_period:(2*ip.n_kin-2)*(m_period+1)] = gamma_arr_tot_max[0:2*ip.n_kin-2]
    
    if m_period == ip.m-1:
        p_arr_pfr_max[(2*ip.n_kin-2)*(m_period+1)] = p_arr_max[0]
        gamma_arr_pfr_max[(2*ip.n_kin-2)*(m_period+1)] = gamma_arr_max[0]
        
    psi_arr_pfr_max = fi.psi_from_momentum(p_arr_pfr_max, ip.iota_a, psi_0)
    alpha_arr_pfr_max = fi.alpha_from_gamma_mod(gamma_arr_pfr_max, ip.m) 

#==============================================================================
# PFR 1 point        
#============================================================================== 

    E_tilde_pfr_1 = 2*g + ip.frac_d_1 * 2*g 
    
    kin_tilde_arr_pfr_1 = np.zeros(ip.n_kin)
    p_arr_pfr_1 = np.zeros(ip.n_kin)
    p_arr_neg_pfr_1 = np.zeros(ip.n_kin)
    
    gamma_arr_pfr_1 = np.zeros(ip.n_kin)
    gamma_arr_neg_pfr_1 = np.zeros(ip.n_kin) 
    
    p_arr_tot_pfr_1 = np.zeros(2*ip.n_kin-2)
    gamma_arr_tot_pfr_1 = np.zeros(2*ip.n_kin-2)
    
    kin_tilde_arr_pfr_1 = fi.poloidal_spacing(E_tilde_pfr_1-2*g, 
                                            E_tilde_pfr_1, ip.n_kin,
                                            ip.alpha_pol)
    
    for i_kin_tilde in range(len(kin_tilde_arr_pfr_1)):
        p_arr_pfr_1[i_kin_tilde] = np.sqrt(2*kin_tilde_arr_pfr_1[i_kin_tilde])
        gamma_arr_pfr_1[i_kin_tilde] = \
            np.arccos(- (E_tilde_pfr_1 - kin_tilde_arr_pfr_1[i_kin_tilde]) / g + 1) \
            + 2 * np.pi * m_period

        gamma_arr_neg_pfr_1[i_kin_tilde] = \
            -np.arccos(- (E_tilde_pfr_1 - kin_tilde_arr_pfr_1[ i_kin_tilde]) / g + 1) \
            + 2 * np.pi * m_period
    
    # from -pi to 0
    p_arr_tot_pfr_1[0:ip.n_kin-1] = p_arr_pfr_1[0:ip.n_kin-1]
    gamma_arr_tot_pfr_1[0:ip.n_kin-1] = gamma_arr_neg_pfr_1[0:ip.n_kin-1]
    
    # from 0 to pi
    p_arr_tot_pfr_1[ip.n_kin-1:2*ip.n_kin-2] = p_arr_pfr_1[ip.n_kin-1:0:-1]
    gamma_arr_tot_pfr_1[ip.n_kin-1:2*ip.n_kin-2] = gamma_arr_pfr_1[ip.n_kin-1:0:-1]
    
    
    p_arr_pfr_pfr_1[(2*ip.n_kin-2)*m_period:(2*ip.n_kin-2)*(m_period+1)] = p_arr_tot_pfr_1[0:2*ip.n_kin-2]
    gamma_arr_pfr_pfr_1[(2*ip.n_kin-2)*m_period:(2*ip.n_kin-2)*(m_period+1)] = gamma_arr_tot_pfr_1[0:2*ip.n_kin-2]
    
    if m_period == ip.m-1:
        p_arr_pfr_pfr_1[(2*ip.n_kin-2)*(m_period+1)] = p_arr_pfr_1[0]
        gamma_arr_pfr_pfr_1[(2*ip.n_kin-2)*(m_period+1)] = gamma_arr_pfr_1[0]
        
    psi_arr_pfr_pfr_1 = fi.psi_from_momentum(p_arr_pfr_pfr_1, ip.iota_a, psi_0)
    alpha_arr_pfr_pfr_1 = fi.alpha_from_gamma_mod(gamma_arr_pfr_pfr_1, ip.m) 

'''
                                      
plt.figure()
plt.scatter(gamma_arr_core_min, p_arr_core_min, marker='.', s=ip.s_size, color='red')     

plt.figure()
plt.scatter(gamma_arr_core_core_1, p_arr_core_core_1, marker='.', s=ip.s_size, color='blue')  

plt.figure()
plt.scatter(gamma_arr_pfr_max, p_arr_pfr_max, marker='.', s=ip.s_size, color='black')                                          
 
plt.figure()
plt.scatter(gamma_arr_pfr_pfr_1, p_arr_pfr_pfr_1, marker='.', s=ip.s_size, color='green')        
'''
