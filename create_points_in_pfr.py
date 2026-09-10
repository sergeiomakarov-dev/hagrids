import numpy as np
import input_params as ip
import matplotlib.pyplot as plt
import fun_time_independent_mh as fi
import create_points_in_island as cpi
import math

def sign(x):
    return int(math.copysign(1, x)) if x != 0 else 0

p_arr_pfr = np.zeros((ip.n_E_pfr, ip.m * (2*ip.n_kin-2)+1))
gamma_arr_pfr = np.zeros((ip.n_E_pfr, ip.m * (2*ip.n_kin-2)+1))

psi_arr_pfr = np.zeros((ip.n_E_pfr, ip.m * (2*ip.n_kin-2)+1))
alpha_arr_pfr = np.zeros((ip.n_E_pfr, ip.m * (2*ip.n_kin-2)+1))

p_arr_tmp_pfr = np.zeros((ip.n_E_pfr, 2*ip.n_kin-1))
gamma_arr_tmp_pfr = np.zeros((ip.n_E_pfr, 2*ip.n_kin-1))


for m_period in range(ip.m):
    g = -ip.A
    psi_0 = (ip.iota_res - ip.iota_b) / ip.iota_a
    
    p_0 = np.sqrt(ip.iota_a) * (psi_0)
    p_max = np.sqrt(ip.iota_a) * (ip.psi_t_max-psi_0) 
    Q_max = p_max ** 2 / 2 - g * fi.A_fun(p_max, p_0, ip.k, ip.flag_type_A_fun)

    Q_X = fi.X_O_energy(0.5, g*ip.k/p_0, g, ip.flag_type_A_fun)
    
    Q_arr = fi.radial_spacing(Q_X, Q_max, ip.n_E_pfr, alpha=ip.alpha_pfr_rad, beta=ip.beta_pfr_rad,
                              flag_spacing_fun=ip.flag_rad_spacing_fun)
    Q_arr[0]=Q_X
    
    kin_tilde_arr = np.zeros((ip.n_E_pfr, ip.n_kin))
    p_arr = np.zeros((ip.n_E_pfr, ip.n_kin))
    p_arr_neg = np.zeros((ip.n_E_pfr, ip.n_kin))
    gamma_arr = np.zeros((ip.n_E_pfr, ip.n_kin))
    gamma_arr_neg = np.zeros((ip.n_E_pfr, ip.n_kin))
    
    p_arr_tot = np.zeros((ip.n_E_pfr, 2*ip.n_kin-2))
    gamma_arr_tot = np.zeros((ip.n_E_pfr, 2*ip.n_kin-2))
    
    for i_E_tilde in range(len(Q_arr)):
        if i_E_tilde == 0:
            p_arr[i_E_tilde,:], _, gamma_arr[i_E_tilde,:], \
                _, gamma_arr_neg[i_E_tilde,:], _ = \
                cpi.gamma_p_in_island(Q_arr[i_E_tilde], Q_X, g, p_0, m_period)
        else:
            if ip.flag_poloidal_initialization == 0 and ip.flag_type_A_fun == 0:
                kin_tilde_arr[i_E_tilde,:] = fi.poloidal_spacing(Q_arr[i_E_tilde]-g, 
                                                                 Q_arr[i_E_tilde]+g, ip.n_kin,
                                                                 alpha=ip.alpha_pol, beta=ip.beta_pol,
                                                                                           flag_spacing_fun=ip.flag_pol_spacing_fun)
                for i_kin_tilde in range(len(kin_tilde_arr[i_E_tilde])):
                    if kin_tilde_arr[i_E_tilde,i_kin_tilde] >= 0:  
                        p_arr[i_E_tilde,i_kin_tilde] = np.sqrt(2*kin_tilde_arr[i_E_tilde,i_kin_tilde])
                    else:
                        p_arr[i_E_tilde,i_kin_tilde] = 0.0
                    temp = - (Q_arr[i_E_tilde] - np.abs(kin_tilde_arr[i_E_tilde,i_kin_tilde])) / \
                               (g * fi.A_fun(p_arr[i_E_tilde,i_kin_tilde], p_0, ip.k, ip.flag_type_A_fun))
                    if temp > 1:
                        temp = 1
                    elif temp < -1:
                        temp = -1
                    gamma_arr[i_E_tilde,i_kin_tilde] = \
                        np.arccos (temp) \
                        + 2 * np.pi * m_period
                    
                    gamma_arr_neg[i_E_tilde,i_kin_tilde] = \
                        -np.arccos (temp) \
                        + 2 * np.pi * m_period
                    if i_kin_tilde == len(kin_tilde_arr[i_E_tilde])-1:
                        gamma_arr[i_E_tilde,i_kin_tilde] = 2 * np.pi * m_period   
                        gamma_arr_neg[i_E_tilde,i_kin_tilde] = 2 * np.pi * m_period   
            elif ip.flag_poloidal_initialization == 1:
                cos_gamma_arr = fi.poloidal_spacing(-1.0, 1.0, 
                                                                 ip.n_kin, alpha=ip.alpha_pol, beta=ip.beta_pol,
                                                                                           flag_spacing_fun=ip.flag_pol_spacing_fun)
                gamma_arr[i_E_tilde,:] = np.arccos (cos_gamma_arr[:]) + 2 * np.pi * m_period
                gamma_arr_neg[i_E_tilde,:] = -np.arccos (cos_gamma_arr[:]) + 2 * np.pi * m_period
                for i_gamma in range(len(gamma_arr[i_E_tilde])):
                    a_tmp = 0.5
                    b_tmp = -np.cos(gamma_arr[i_E_tilde,i_gamma]) * g * ip.k / p_0
                    c_tmp = -Q_arr[i_E_tilde] - np.cos(gamma_arr[i_E_tilde,i_gamma]) * g
                    if fi.check_discriminant(a_tmp, b_tmp, c_tmp, ip.flag_type_A_fun) < 0:
                        p_arr[i_E_tilde,i_gamma], p_arr_neg[i_E_tilde,i_gamma] = fi.solut_same_roots(a_tmp, b_tmp, c_tmp, ip.flag_type_A_fun)
                    else:
                        p_arr[i_E_tilde,i_gamma], p_arr_neg[i_E_tilde,i_gamma] = fi.solut_algeb_eq(a_tmp, b_tmp, c_tmp, ip.flag_type_A_fun) 
            else :
                raise ValueError("Invalid value for flag_poloidal_initialization. Must be 0 or 1. The flag_poloidal_initialization == 0 option is not supported together with flag_type_A_fun == 1") 
            
    # from -pi to 0
    p_arr_tot[:, 0:ip.n_kin-1] = p_arr[:, 0:ip.n_kin-1]
    gamma_arr_tot[:, 0:ip.n_kin-1] = gamma_arr_neg[:, 0:ip.n_kin-1]
    
    # from 0 to pi
    p_arr_tot[:, ip.n_kin-1:2*ip.n_kin-2] = p_arr[:, ip.n_kin-1:0:-1]
    gamma_arr_tot[:, ip.n_kin-1:2*ip.n_kin-2] = gamma_arr[:, ip.n_kin-1:0:-1]
    
    
    p_arr_pfr[:, (2*ip.n_kin-2)*m_period:(2*ip.n_kin-2)*(m_period+1)] = p_arr_tot[:, 0:2*ip.n_kin-2]
    gamma_arr_pfr[:, (2*ip.n_kin-2)*m_period:(2*ip.n_kin-2)*(m_period+1)] = gamma_arr_tot[:, 0:2*ip.n_kin-2]
    
    if m_period == 0:
        # from -pi to 0
        p_arr_tmp_pfr[:, 0:ip.n_kin-1] = p_arr[:, 0:ip.n_kin-1]
        gamma_arr_tmp_pfr[:, 0:ip.n_kin-1] = gamma_arr_neg[:, 0:ip.n_kin-1]
    
        # from 0 to pi
        p_arr_tmp_pfr[:, ip.n_kin-1:2*ip.n_kin-1] = p_arr[:, ip.n_kin-1::-1]
        gamma_arr_tmp_pfr[:, ip.n_kin-1:2*ip.n_kin-1] = gamma_arr[:, ip.n_kin-1::-1]
    
    if m_period == ip.m-1:
        p_arr_pfr[:, (2*ip.n_kin-2)*(m_period+1)] = p_arr[:, 0]
        gamma_arr_pfr[:, (2*ip.n_kin-2)*(m_period+1)] = gamma_arr[:, 0]
        
    psi_arr_pfr = fi.psi_from_momentum(p_arr_pfr, ip.iota_a, psi_0)
    alpha_arr_pfr = fi.alpha_from_gamma_mod(gamma_arr_pfr, ip.m) 


#plt.figure()
#plt.scatter(gamma_arr_tmp_pfr, p_arr_tmp_pfr, marker='.', s=ip.s_size, color='red')
#np.savez('pfr_data.npz', gamma_pfr=gamma_arr_tmp_pfr, p_pfr=p_arr_tmp_pfr)