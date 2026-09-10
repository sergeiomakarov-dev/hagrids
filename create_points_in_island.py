import numpy as np
import input_params as ip
import matplotlib.pyplot as plt
import fun_time_independent_mh as fi
import math

def sign(x):
    return int(math.copysign(1, x)) if x != 0 else 0

def gamma_p_in_island(Q, Q_X, g, p_0, m_period):
    kin_tilde_arr = np.zeros((ip.n_kin))
    p_arr = np.zeros((ip.n_kin))
    p_arr_neg = np.zeros((ip.n_kin))
    gamma_arr = np.zeros((ip.n_kin))
    gamma_arr_neg = np.zeros((ip.n_kin))
    gamma_arr_2 = np.zeros((ip.n_kin))
    gamma_arr_neg_2 = np.zeros((ip.n_kin))
    
    if ip.flag_poloidal_initialization == 0 and ip.flag_type_A_fun == 0:
         kin_tilde_arr[:] = fi.poloidal_spacing(0, g+Q, ip.n_kin, alpha=ip.alpha_pol)
         for i_kin_tilde in range(len(kin_tilde_arr)):
             if i_kin_tilde == len(kin_tilde_arr)-1:
                 p_arr[i_kin_tilde] = np.sqrt(
                     2*kin_tilde_arr[i_kin_tilde])
                 gamma_arr[i_kin_tilde] = 2 * np.pi * m_period
                 gamma_arr_neg[i_kin_tilde] = 2 * \
                     np.pi * m_period
                 p_arr_neg[i_kin_tilde] = - \
                     np.sqrt(2*kin_tilde_arr[i_kin_tilde])
                 gamma_arr_2[i_kin_tilde] = 2 * np.pi * m_period
                 gamma_arr_neg_2[i_kin_tilde] = 2 * np.pi * m_period
    
             else:
                 p_arr[i_kin_tilde] = np.sqrt(
                     2*kin_tilde_arr[i_kin_tilde])
                 gamma_arr[i_kin_tilde] = \
                     np.arccos(- (Q - kin_tilde_arr[i_kin_tilde]) / g) \
                     + 2 * np.pi * m_period
    
                 gamma_arr_neg[i_kin_tilde] = \
                     -np.arccos(- (Q - kin_tilde_arr[i_kin_tilde]) / g) \
                     + 2 * np.pi * m_period
    
                 p_arr_neg[i_kin_tilde] = - \
                     np.sqrt(2*kin_tilde_arr[i_kin_tilde])
    
                 gamma_arr_2[i_kin_tilde] = \
                     np.arccos(- (Q - kin_tilde_arr[i_kin_tilde]) / g) \
                     + 2 * np.pi * m_period
    
                 gamma_arr_neg_2[i_kin_tilde] = \
                     -np.arccos(- (Q - kin_tilde_arr[i_kin_tilde]) / g) \
                     + 2 * np.pi * m_period
    
    elif ip.flag_poloidal_initialization == 1:
         if Q == Q_X:
             temp, max_gamma_p_1 = fi.solut_same_roots(
                 0.5, g * ip.k / p_0, -Q + g, ip.flag_type_A_fun)
             max_gamma_1 = np.pi \
                 + 2 * np.pi * m_period
         else:
             max_gamma_p_1 = ip.flag_type_A_fun * \
                 (- (p_0/ip.k) + np.sqrt((p_0/ip.k)
                  ** 2 - 2 * Q))
             max_gamma_1 = np.arccos (- (Q - max_gamma_p_1 ** 2 / 2) /
                                                 (g * fi.A_fun(max_gamma_p_1, p_0, ip.k, ip.flag_type_A_fun)) ) \
                 + 2 * np.pi * m_period
         cos_gamma_arr = fi.poloidal_spacing(np.cos(max_gamma_1), 1.0,
                                             ip.n_kin, alpha=ip.alpha_pol)
         gamma_arr[:] = np.arccos(cos_gamma_arr[:]) + 2 * np.pi * m_period
         gamma_arr_2[:] = np.arccos(cos_gamma_arr[:]) + 2 * np.pi * m_period
         gamma_arr_neg[:] = -np.arccos(cos_gamma_arr[:]) + 2 * np.pi * m_period
         gamma_arr_neg_2[:] = -np.arccos(cos_gamma_arr[:]) + 2 * np.pi * m_period
         for i_gamma in range(len(gamma_arr)):
             a_tmp = 0.5
             b_tmp = - \
                 np.cos(gamma_arr_2[i_gamma]) * g * ip.k / p_0
             c_tmp = -Q - \
                 np.cos(gamma_arr_2[i_gamma]) * g
             if i_gamma == 0:
                 p_arr[i_gamma], p_arr_neg[i_gamma] = fi.solut_same_roots(a_tmp, b_tmp, c_tmp, ip.flag_type_A_fun)
             else:
                 p_arr[i_gamma], p_arr_neg[i_gamma] = fi.solut_algeb_eq(a_tmp, b_tmp, c_tmp, ip.flag_type_A_fun)
    else:
         raise ValueError(
             "Invalid value for flag_poloidal_initialization. Must be 0 or 1. The flag_poloidal_initialization == 0 option is not supported together with flag_type_A_fun == 1")
    return p_arr, p_arr_neg, gamma_arr, gamma_arr_2, gamma_arr_neg, gamma_arr_neg_2


psi_arr_isl = []
alpha_arr_isl = []

for m_period in range(ip.m):
#for m_period in range(1):
    g = -ip.A
    psi_0 = (ip.iota_res - ip.iota_b) / ip.iota_a
    
    p_0 = np.sqrt(ip.iota_a) * (psi_0)
    Q_O = fi.X_O_energy(0.5, -g*ip.k/p_0, -g, ip.flag_type_A_fun)
    Q_X = fi.X_O_energy(0.5, g*ip.k/p_0, g, ip.flag_type_A_fun)
    
    Q_arr = fi.radial_spacing(Q_X, Q_X/ip.min_E_fac+Q_O, ip.n_E, alpha=ip.alpha_island_rad)
    
    cos_gamma_arr = np.zeros(ip.n_kin)
    p_arr = np.zeros((ip.n_E, ip.n_kin))
    p_arr_neg = np.zeros((ip.n_E, ip.n_kin))
    gamma_arr = np.zeros((ip.n_E, ip.n_kin))
    gamma_arr_neg = np.zeros((ip.n_E, ip.n_kin))
    gamma_arr_2 = np.zeros((ip.n_E, ip.n_kin))
    gamma_arr_neg_2 = np.zeros((ip.n_E, ip.n_kin))
    
    p_arr_tot = np.zeros((ip.n_E, 4*ip.n_kin-3))
    gamma_arr_tot = np.zeros((ip.n_E, 4*ip.n_kin-3))
    
    
    for i_E_tilde in range(len(Q_arr)):
        p_arr[i_E_tilde,:], p_arr_neg[i_E_tilde,:], gamma_arr[i_E_tilde,:], \
            gamma_arr_2[i_E_tilde,:], gamma_arr_neg[i_E_tilde,:], gamma_arr_neg_2[i_E_tilde,:] = \
            gamma_p_in_island(Q_arr[i_E_tilde], Q_X, g, p_0, m_period)

    # from -pi to 0
    p_arr_tot[:, 0:ip.n_kin-1] = p_arr_neg[:, 0:ip.n_kin-1]
    gamma_arr_tot[:, 0:ip.n_kin-1] = gamma_arr_neg_2[:, 0:ip.n_kin-1]
    
    # from 0 to pi
    p_arr_tot[:, ip.n_kin-1:2*ip.n_kin-2] = p_arr_neg[:, ip.n_kin-1:0:-1]
    gamma_arr_tot[:, ip.n_kin-1:2*ip.n_kin-2] = gamma_arr_2[:, ip.n_kin-1:0:-1]

    # from pi to 0
    p_arr_tot[:, 2*ip.n_kin-2:3*ip.n_kin-3]= p_arr[:, 0:ip.n_kin-1]
    gamma_arr_tot[:, 2*ip.n_kin-2:3*ip.n_kin-3] = gamma_arr[:, 0:ip.n_kin-1]
    
    # from 0 to -pi
    p_arr_tot[:, 3*ip.n_kin-3:4*ip.n_kin-4]= p_arr[:, ip.n_kin-1:0:-1]
    gamma_arr_tot[:, 3*ip.n_kin-3:4*ip.n_kin-4] = gamma_arr_neg[:, ip.n_kin-1:0:-1]    

    # at -pi
    p_arr_tot[:, 4*ip.n_kin-4]= p_arr[:, 0]
    gamma_arr_tot[:, 4*ip.n_kin-4] = gamma_arr_neg[:, 0]

    psi_arr = fi.psi_from_momentum(p_arr_tot, ip.iota_a, psi_0)
    alpha_arr = fi.alpha_from_gamma_mod(gamma_arr_tot, ip.m) 
    
    psi_arr_isl.append(psi_arr)
    alpha_arr_isl.append(alpha_arr)
 
    if m_period == 0:
        
        p_arr_tmp_island = p_arr_tot
        gamma_arr_tmp_island = gamma_arr_tot

#plt.figure()
#plt.scatter(gamma_arr_tmp_island, p_arr_tmp_island, marker='.', s=ip.s_size, color='red')
#np.savez('island_data.npz', gamma_island=gamma_arr_tmp_island, p_island=p_arr_tmp_island)
