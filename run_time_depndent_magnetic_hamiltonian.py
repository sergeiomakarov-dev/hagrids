#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.special import ellipj, ellipk, ellipkinc
from scipy.integrate import solve_ivp
import input_params as ip

# functions     

def polar_to_cartesian(r, theta):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

def dual_stoermer_verlet(t, p, gamma, momentum_rhs):
    dt = t[1] - t[0] 
    
    for i in range(len(t)-1):
        
        # First half step for momentum
        a = momentum_rhs(gamma[i], t[i])
        p_half = p[i] + 0.5 * a * dt
        
        # Full step for positions
        gamma[i + 1] = gamma[i] + (ip.iota_a * p_half + ip.iota_b)* dt 
        
        # Second half step for momentum
        a_next = momentum_rhs(gamma[i + 1], t[i + 1])
        p[i + 1] = p_half + 0.5 * a_next * dt
        
    return p, gamma
    
def num_solver(t, y1, momentum_rhs, flag_sv = 0, flag_sn = 0, method = 'RK45'):
    """
    Perform dual Störmer-Verlet integration for a given number of steps.
    
    """
    
    # Initialize arrays
    yout = np.zeros((2, len(t)))
    y_wrk = np.zeros((2, len(t)))
    y_wrk[0, 0] = y1[0, 0]
    y_wrk[1, 0] = y1[1, 0]
    
    if flag_sv:
        tmp_p, tmp_gamma = dual_stoermer_verlet(t, y_wrk[0, :], y_wrk[1, :], momentum_rhs) 
        y_wrk[0, :] = tmp_p
        y_wrk[1, :] = tmp_gamma
        
    if flag_sn:
        y_temp = solve_ivp(hamiltonian_time_dependent, [t[0], t[-1]], [y1[0, 0], y1[1, 0]], args=(ip.A1, ip.iota_a, ip.iota_b, ip.m1, ip.iota_res1), t_eval=t, method=method)
        y_wrk = y_temp.y
        
    yout[0, :] = y_wrk[0, :] 
    yout[1, :] = y_wrk[1, :]
        
    return yout 

def momentum_rhs(x,t):
    return ip.m1 * ip.A1 * np.sin(ip.m1 * (x - ip.iota_res1 * t)) + ip.m2 * ip.A2 * np.sin(ip.m2 * (x - ip.iota_res2 * t)) 

def hamiltonian_time_dependent(t, z, A1, iota_a, iota_b, m1, iota_res1):
    p, q = z
    dpdt = m1 * A1 * np.sin(m1 * (q - iota_res1 * t))
    dqdt = iota_a * p + iota_b
    
    return [dpdt, dqdt]

momentum = lambda E_tilde_loc, g_loc, gamma_loc : np.sqrt(2 * (E_tilde_loc - g_loc * (1 - np.cos(gamma_loc))))           
psi_from_momentum = lambda p_loc, iota_a_loc, psi_0_loc : p_loc / np.sqrt(iota_a_loc) + psi_0_loc   
alpha_from_gamma = lambda gamma_loc, m_period_loc, m_loc : (gamma_loc + 2 * np.pi * m_period_loc)/ m_loc
alpha_from_gamma_mod = lambda gamma_loc, m_loc : gamma_loc / m_loc

def psi_theta_to_R_Z(psi,theta):
    r = np.sqrt(psi / (ip.B_0 * np.pi))
    
    # In cylindrical coordinates
    R = ip.R_0 + r * np.cos(theta)
    Z = r * np.sin(theta)
    return R, Z

def from_phase_to_real(phi_span_loc, y_loc,
                       R_lines_loc, Z_lines_loc, n_per_turn = 1):

    R_one_turn =  np.array([])
    Z_one_turn =  np.array([])

    psi = y_loc [0, :]
    theta = y_loc [1, :]

    R, Z = psi_theta_to_R_Z(psi,theta)
    
    i_temp = 0
    for i_phi in range(0,len(phi_span_loc)-ip.points,int(ip.points/n_per_turn)):
        R_one_turn = np.append(R_one_turn, R[i_phi])
        Z_one_turn = np.append(Z_one_turn, Z[i_phi])
        i_temp += 1
    
    R_lines_loc.append(R_one_turn) 
    Z_lines_loc.append(Z_one_turn)    
    
    # In cartisian coordinates
    X, Y = polar_to_cartesian(R, phi_span_loc)
    return R_lines_loc, Z_lines_loc

def loc_iota_calc(phi_span_loc, y_loc, iota_loc, iota_max_loc, iota_min_loc):       
    
    iota_loc_one_turn =  np.array([])

    #psi = y_loc [0, :]
    theta = y_loc [1, :]
    
    i_temp = 0    
    for i_phi in range(0,len(phi_span_loc)-ip.points,ip.points):
        d_thata = theta[i_phi+1] - theta[i_phi]
        d_phi = phi_span_loc[i_phi+1] - phi_span_loc[i_phi]
        
        tmp = d_thata/d_phi
        iota_loc_one_turn = np.append(iota_loc_one_turn, tmp)
        iota_max_loc = max(iota_max_loc, tmp)
        iota_min_loc = min(iota_min_loc, tmp)
        i_temp += 1
        
    iota_loc.append(iota_loc_one_turn) 
        
    return iota_loc, iota_max_loc, iota_min_loc

########## main body
# Phi span and initial alpha and psi points
phi_span_td = np.linspace(0, ip.n_turns * 2 * np.pi, ip.n_turns * ip.points + 1)
psi_t_span = np.linspace(ip.psi_t_min, ip.psi_t_max, ip.n_psi_surf_td)

temp = np.pi / ip.n_lines_p_surf_td
theta_0_surf = np.linspace(-np.pi + temp, np.pi - temp, ip.n_lines_p_surf_td)

theta_0 = np.array([])
psi_t_0 = np.array([])

for i_psi_surf in range(ip.n_psi_surf_td):
    theta_0 = np.append(theta_0, theta_0_surf)
    psi_t_0 = np.append(psi_t_0, psi_t_span[i_psi_surf] * np.ones(ip.n_lines_p_surf_td))

# Line tracing     
R_lines_sv_td = []
Z_lines_sv_td = []

R_lines_sn_td = []
Z_lines_sn_td = []

iota_sv_td = []
iota_sv_td_max = float('-inf')
iota_sv_td_min = float('inf')

for i_lines in range(ip.n_psi_surf_td * ip.n_lines_p_surf_td):
    y1 = np.array([[psi_t_0[i_lines]],[theta_0[i_lines]]])

    if ip.flag_Stoermer_Verlet_td:
        y_sv_td = num_solver(phi_span_td, y1, momentum_rhs, flag_sv = ip.flag_Stoermer_Verlet_td, method = ip.method)
        
        R_lines_sv_td, Z_lines_sv_td = \
        from_phase_to_real(phi_span_td, y_sv_td, 
                           R_lines_sv_td, Z_lines_sv_td)
        
        iota_sv_td, iota_sv_td_max, iota_sv_td_min = \
            loc_iota_calc(phi_span_td, y_sv_td, iota_sv_td,
                          iota_sv_td_max, iota_sv_td_min) 
        
        
    if ip.flag_scipy_num_td:
        y_sn_td = num_solver(phi_span_td, y1, momentum_rhs, flag_sn = ip.flag_scipy_num_td, method = ip.method)
        
        R_lines_sn_td, Z_lines_sn_td = \
        from_phase_to_real(phi_span_td, y_sn_td, 
                           R_lines_sn_td, Z_lines_sn_td)

psi_test = []
theta_test = []
R_test = []  
Z_test = []
if ip.flag_Stoermer_Verlet_td:
    psi_test = np.linspace(ip.psi_t_min, ip.psi_t_max, ip.n_test_points)
    theta_test = np.linspace(0, 0, ip.n_test_points)
    R_test, Z_test = psi_theta_to_R_Z(psi_test, theta_test)
    
for i_lines in range(ip.n_psi_surf_td * ip.n_lines_p_surf_td):
    if ip.flag_Stoermer_Verlet_td and ip.flag_filter:
        indices_to_include = list(range(0, len(R_lines_sv_td[i_lines]), ip.n_filter))
        indices_to_exclude_1 = [x for x in list(range(0, len(R_lines_sv_td[i_lines]))) if x not in indices_to_include]
        R_lines_sv_td[i_lines] = np.delete(R_lines_sv_td[i_lines], indices_to_exclude_1)
        Z_lines_sv_td[i_lines] = np.delete(Z_lines_sv_td[i_lines], indices_to_exclude_1)
        iota_sv_td[i_lines] = np.delete(iota_sv_td[i_lines], indices_to_exclude_1)
        
    if ip.flag_Stoermer_Verlet_td and ip.flag_filter:
        indices_to_exclude = []
        for i_points in range(len(R_lines_sv_td[i_lines])):
            if R_lines_sv_td[i_lines][i_points] < ip.R_min or \
                R_lines_sv_td[i_lines][i_points] > ip.R_max or \
                    Z_lines_sv_td[i_lines][i_points] < ip.Z_min or \
                        Z_lines_sv_td[i_lines][i_points] > ip.Z_max:
                           indices_to_exclude.append(i_points)
        
        R_lines_sv_td[i_lines] = np.delete(R_lines_sv_td[i_lines], indices_to_exclude)
        Z_lines_sv_td[i_lines] = np.delete(Z_lines_sv_td[i_lines], indices_to_exclude)
        iota_sv_td[i_lines] = np.delete(iota_sv_td[i_lines], indices_to_exclude)

R_configs_prep = []  
Z_configs_prep = []

R_configs = []
Z_configs = []

psi_0_configs = []

if ip.flag_config and ip.flag_Stoermer_Verlet_td:
    psi_0_configs = np.linspace(ip.psi_0_min_configs, ip.psi_0_max_configs, ip.n_configs)
    theta_0_configs = np.linspace(ip.theta_0_min_configs, ip.theta_0_max_configs, ip.n_configs)
    psi_0_configs = np.append(psi_0_configs, ip.psi_0_o_point)
    theta_0_configs = np.append(theta_0_configs, ip.theta_0_o_point)
    for i_lines in range(len(psi_0_configs)):
        y1 = np.array([[psi_0_configs[i_lines]],[theta_0_configs[i_lines]]])
        
        if ip.flag_Stoermer_Verlet_td:
            y_configs = num_solver(phi_span_td, y1, momentum_rhs, flag_sv = ip.flag_Stoermer_Verlet_td, method = ip.method)
        
            R_configs_prep, Z_configs_prep = \
                from_phase_to_real(phi_span_td, y_configs, R_configs_prep, Z_configs_prep, n_per_turn = 5)  
                
    R_configs = [np.zeros(len(psi_0_configs)) for _ in range(len(R_configs_prep[0]))]   
    Z_configs = [np.zeros(len(psi_0_configs)) for _ in range(len(R_configs_prep[0]))]  
    psi_0_configs_color = [np.zeros(len(psi_0_configs)) for _ in range(len(R_configs_prep[0]))]        
    for i_lines in range(len(psi_0_configs)):
        wrk1 = R_configs_prep[i_lines]
        wrk2 = Z_configs_prep[i_lines]
        wrk3 = psi_0_configs[i_lines]
        for i_points in range(len(wrk1)):
            R_configs[i_points][i_lines] = wrk1[i_points]
            Z_configs[i_points][i_lines] = wrk2[i_points]
            #psi_0_configs_color[i_points][i_lines] = wrk3[i_points]
            
psi_iota_res1_val = (ip.iota_res1 - ip.iota_b) / ip.iota_a
psi_iota_res2_val = (ip.iota_res2 - ip.iota_b) / ip.iota_a

psi_iota_res1 = np.linspace(psi_iota_res1_val, psi_iota_res1_val, ip.n_points_iota_res)
psi_iota_res2 = np.linspace(psi_iota_res2_val, psi_iota_res2_val, ip.n_points_iota_res)
theta_0_res = np.linspace(-np.pi, np.pi , ip.n_points_iota_res)

R_iota_res1, Z_iota_res1 = psi_theta_to_R_Z(psi_iota_res1, theta_0_res)
R_iota_res2, Z_iota_res2 = psi_theta_to_R_Z(psi_iota_res2, theta_0_res)