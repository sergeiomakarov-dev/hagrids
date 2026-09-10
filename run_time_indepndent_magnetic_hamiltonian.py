#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.special import ellipj, ellipk, ellipkinc
from scipy.integrate import solve_ivp
import input_params as ip
import create_points_in_island as cp

# functions
def pendulum(t, y1):
    '''Classical pendulum solution'''

    
    # The initial pololdal angle gamma_1, the initial momentum p_1, energy E_tilde and corresponding k
    x_point = 0
    
    g = -ip.A
    psi_0 = (ip.iota_res - ip.iota_b) / ip.iota_a
    gamma_1_tot = ip.m * y1[1, 0]
    m_period = np.floor((gamma_1_tot + np.pi) / (2 * np.pi))
    gamma_1 = gamma_1_tot - 2 * np.pi * m_period 
    
    if gamma_1 >= np.pi or gamma_1 < -np.pi:
        raise ValueError("gamma_1 >= np.pi or gamma_1 < -np.pi")
        
    p_1 = np.sqrt(ip.iota_a) * (y1[0, 0]-psi_0) 
    E_tilde = p_1 ** 2 / 2 + g * (1 - np.cos(gamma_1))
    k = np.sqrt(E_tilde / (2 * g))
    
    # The period T
    if k < 1:
        uT = 4 * ellipk(k)
    elif k > 1:
        uT = 2 / k * ellipk(1 / k)
    elif k == 1:
        uT = 4 * ellipk(k)
        
    T = np.sqrt(1 / g) * uT
    
    # The intial time shift with respect to the time where phase = 0
    if k < 1:
        xi_1 = compute_xi_1(gamma_1,k)
        t_phase = ellipkinc(xi_1, k) / (np.sqrt(g))
    elif k > 1:
        xi_1 = gamma_1 / 2
        t_phase = ellipkinc(xi_1, 1/k) / (np.sqrt(g) * k)
    elif k == 1:
        if gamma_1 == -np.pi:
            x_point = 1
        else:
            x_point = 0  
        xi_1 = gamma_1 / 2
        t_phase = ellipkinc(xi_1, 1/k) / (np.sqrt(g) * k)


    # The sign of the initial momentum      
    if p_1 > 0:
        p_p_init = 1
    else:
        p_p_init = -1
        
    # Modified time array
    u = (p_p_init * (t - t[0]) + t_phase) * np.sqrt(g);
    
    # Initialize arrays
    gamma = np.zeros(len(t))
    gamma_period = np.zeros(len(t))
    p = np.zeros(len(t))
    yout = np.zeros((2, len(t)))
    
    # The classical pendulum solution for gamma and p
    if np.isreal(xi_1):
        if k < 1:
            for i in range(len(t)):   
                sn_u, cn_u, dn_u, _ = ellipj(u[i], k)
                gamma[i] = 2 * np.arcsin(k * sn_u)
                
                t_period = np.floor(u[i] / uT)
                if (np.mod(u[i], uT) <= uT / 4) or (np.mod(u[i], uT) >= 3 * uT / 4):
                    p[i] = p_p_init * momentum(E_tilde, g, gamma[i])
                else:
                    p[i] = -1 * p_p_init * momentum(E_tilde, g, gamma[i])

                if np.isreal(p[i]):
                    yout[0, i] = psi_from_momentum(p[i], ip.iota_a, psi_0)
                else:                    
                    yout[0, i] = np.nan
                yout[1, i] = alpha_from_gamma(gamma[i], m_period, ip.m) 
        elif k > 1:
            for i in range(len(t)):
                sn_u, cn_u, dn_u, _ = ellipj(k * u[i], 1 / k)
                gamma_period[i] = 2 * np.arcsin(sn_u)

                t_period = np.floor((u[i] + uT / 2) / uT)
                if np.mod(t_period, 2) == 0:
                    gamma[i] = 2 * np.pi * t_period + gamma_period[i]
                else:
                    gamma[i] = 2 * np.pi * t_period - gamma_period[i]

                p[i] = p_p_init*momentum(E_tilde, g, gamma[i])

                if np.isreal(p[i]):
                    yout[0, i] = psi_from_momentum(p[i], ip.iota_a, psi_0)
                else:
                    yout[0, i] = np.nan
                yout[1, i] = alpha_from_gamma(gamma[i], m_period, ip.m)    
        elif k == 1:
            if x_point == 1:
                for i in range(len(t)):   
                    gamma[i] = gamma_1
                    p[i] = p_1
                    if np.isreal(p[i]):
                        yout[0, i] = psi_from_momentum(p[i], ip.iota_a, psi_0)
                    else:                    
                        yout[0, i] = np.nan
                    yout[1, i] = alpha_from_gamma(gamma[i], m_period, ip.m)
            elif x_point == 0:
                for i in range(len(t)):
                    sn_u, cn_u, dn_u, _ = ellipj(k * u[i], 1 / k)
                    gamma[i] = 2 * np.arcsin(sn_u)
                    p[i] = p_p_init*momentum(E_tilde, g, gamma[i])
                    yout[0, i] = psi_from_momentum(p[i], ip.iota_a, psi_0)
                    yout[1, i] = alpha_from_gamma(gamma[i], m_period, ip.m)   
    else:
        yout[0, :] = np.nan
        yout[1, :] = np.nan
    return yout 

def compute_xi_1(gamma_1, k):
    # Ensure inputs are numpy arrays for element-wise operations
    gamma_1 = np.asarray(gamma_1)
    k = np.asarray(k)
    
    # Compute the argument for arcsin
    argument = np.sin(gamma_1 / 2) / k
    
    # Clip argument to ensure it's within [-1, 1]
    argument = np.clip(argument, -1, 1)
    
    # Compute xi_1 with the clipped argument
    xi_1 = np.arcsin(argument)
    
    return xi_1       

def polar_to_cartesian(r, theta):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

def dual_stoermer_verlet(t, p, gamma, momentum_rhs):
    dt = t[1] - t[0] 
    
    for i in range(len(t)-1):
        
        # First half step for momentum
        a = momentum_rhs(gamma[i])
        p_half = p[i] + 0.5 * a * dt
        
        # Full step for positions
        gamma[i + 1] = gamma[i] + p_half * dt
        
        # Second half step for momentum
        a_next = momentum_rhs(gamma[i + 1])
        p[i + 1] = p_half + 0.5 * a_next * dt
    return p, gamma
    
def num_solver(t, y1, momentum_rhs, flag_sv = 0, flag_sn = 0, method = 'RK45'):
    """
    Perform dual Störmer-Verlet integration for a given number of steps.
    
    """
    
    # Initialize arrays
    gamma = np.zeros(len(t))

    p = np.zeros(len(t))
    yout = np.zeros((2, len(t)))
    
    psi_0 = (ip.iota_res - ip.iota_b) / ip.iota_a
    gamma[0] = ip.m * y1[1, 0]    
    p[0] = np.sqrt(ip.iota_a) * (y1[0, 0]-psi_0) 
    
    if flag_sv:
        p, gamma = dual_stoermer_verlet(t, p, gamma, momentum_rhs)   
        
    if flag_sn:
        y_temp = solve_ivp(hamiltonian_pendulum, [t[0], t[-1]], [p[0], gamma[0]], args=(ip.A,), t_eval=t, method=method)
        p, gamma = y_temp.y
        
    yout[0, :] = psi_from_momentum(p[:], ip.iota_a, psi_0)
    yout[1, :] = alpha_from_gamma_mod(gamma[:], ip.m)   
        
    return yout 

def momentum_rhs(x):
    g = -ip.A
    return -g * np.sin(x)

def hamiltonian_pendulum(t, z, A):
    g = -A
    p, q = z
    dpdt = -g * np.sin(q)
    dqdt = p
    
    return [dpdt, dqdt]

#momentum = lambda E_tilde_loc, g_loc, gamma_loc : np.sqrt(2 * (E_tilde_loc - g_loc * (1 - np.cos(gamma_loc)))) 
def momentum(E_tilde_loc, g_loc, gamma_loc):
    # Ensure inputs are numpy arrays for element-wise operations
    E_tilde_loc = np.asarray(E_tilde_loc)
    g_loc = np.asarray(g_loc)
    gamma_loc = np.asarray(gamma_loc)
    
    # Compute the argument inside the square root
    argument = 2 * (E_tilde_loc - g_loc * (1 - np.cos(gamma_loc)))
    
    # Handle cases where argument is negative and very close to zero
    mask_small_negative = (argument < 0) & (np.abs(argument) < 1e-15)
    argument = np.where(mask_small_negative, 0, argument)  # Replace small negative values with 0
    
    # Compute momentum with the adjusted argument
    result = np.sqrt(argument)
    
    return result 
         
psi_from_momentum = lambda p_loc, iota_a_loc, psi_0_loc : p_loc / np.sqrt(iota_a_loc) + psi_0_loc   
alpha_from_gamma = lambda gamma_loc, m_period_loc, m_loc : (gamma_loc + 2 * np.pi * m_period_loc)/ m_loc
alpha_from_gamma_mod = lambda gamma_loc, m_loc : gamma_loc / m_loc

def psi_theta_to_R_Z(psi,theta):
    r = np.sqrt(psi / (ip.B_0 * np.pi))
    
    # In cylindrical coordinates
    R = ip.R_0 + r * np.cos(theta)
    Z = r * np.sin(theta)
    return R, Z

def R_Z_to_psi_theta(R,Z):
    # In toroidal coordinates
    r = np.hypot(R - ip.R_0, Z)
    theta = np.arctan2(Z, R - ip.R_0)
    
    psi = r**2 * (ip.B_0 * np.pi)
    return psi, theta

def from_phase_to_real(phi_span_loc, y_loc, i_out_loc):
    
    R_one_turn =  np.array([])
    Z_one_turn =  np.array([])
    theta_one_turn =  np.array([])

    psi = y_loc [0, :]
    alpha = y_loc [1, :]

    psi_iota_res = (ip.iota_res - ip.iota_b) / ip.iota_a
    
    # In phase varibles
    gamma = ip.m * alpha
    p_gamma = np.sqrt(ip.iota_a) * (psi - psi_iota_res)
    
    # In the lab frame
    theta = alpha + ip.iota_res * phi_span_loc

    R, Z = psi_theta_to_R_Z(psi,theta)
    
    i_temp = 0
    for i_phi in i_out_loc:
        R_one_turn = np.append(R_one_turn, R[i_phi])
        Z_one_turn = np.append(Z_one_turn, Z[i_phi])
        theta_one_turn = np.append(theta_one_turn, theta[i_phi])
        i_temp += 1    
    
    # In cartisian coordinates
    X, Y = polar_to_cartesian(R, phi_span_loc)
    return R_one_turn, Z_one_turn, theta_one_turn 

def normalize_angle(angle):
    return angle % (2 * np.pi)

########## main body
# Phi span and initial alpha and psi points
phi_span = np.linspace(0, ip.n_turns * 2 * np.pi / ip.nfp, ip.n_turns * ip.points + 1)   
 
alpha_0 = np.array([])
psi_t_0 = np.array([])

if ip.flag_type_init == 0:
    psi_t_span = np.linspace(ip.psi_t_min, ip.psi_t_max, ip.n_psi_surf)

    temp = np.pi / ip.n_lines_p_surf
    alpha_0_surf = np.linspace(-np.pi + temp, np.pi - temp, ip.n_lines_p_surf)
    


    for i_psi_surf in range(ip.n_psi_surf):
        alpha_0 = np.append(alpha_0, alpha_0_surf)
        psi_t_0 = np.append(psi_t_0, psi_t_span[i_psi_surf] * np.ones(ip.n_lines_p_surf))
    
elif ip.flag_type_init == 1:
    R_span = np.linspace(ip.R_start, ip.R_end, ip.n_lines_tot)
    Z_span = np.linspace(ip.Z_const, ip.Z_const, ip.n_lines_tot)
    
    for i_lines in range(ip.n_lines_tot):
        psi_t_0_temp, alpha_0_temp = R_Z_to_psi_theta(R_span[i_lines], Z_span[i_lines])
        alpha_0 = np.append(alpha_0, alpha_0_temp)
        psi_t_0 = np.append(psi_t_0, psi_t_0_temp)
        
elif ip.flag_type_init == 2:
    with open('base0.dat', 'r') as f:
    
        lines = f.readlines()
        
    # Initialize a list to store processed and formatted float values
    processed_data = []
    
    # Process each line
    for line in lines:
        # Strip leading/trailing whitespace from the line
        stripped_line = line.strip()
        
        # Check if the stripped line is not empty and does not start with '#'
        if stripped_line and stripped_line[0] != '#':
            try:
                # Split the line into individual elements
                elements = stripped_line.split()
                
                # Convert each element to a float and format it
                formatted_values = [f"{float(element):.10f}" for element in elements]
                
                # Append the list of formatted float values to processed_data
                processed_data.append(formatted_values)
            
            except ValueError as e:
                # Handle the error (e.g., log it, print it, or raise it)
                print(f"Error converting line to floats: {stripped_line}")
                print(f"Exception: {e}")
    
    R_span = np.array([])  
    Z_span = np.array([])            
    
    # Print the processed and formatted data
    for i in range(len(processed_data)):
    #for i in range(100):
        R_span = np.append(R_span, float(processed_data[i][0]))
        Z_span = np.append(Z_span, float(processed_data[i][1]))
        
    ip.n_lines_tot = len(R_span)
    
    for i_lines in range(ip.n_lines_tot):
        psi_t_0_temp, alpha_0_temp = R_Z_to_psi_theta(R_span[i_lines], Z_span[i_lines])
        alpha_0 = np.append(alpha_0, alpha_0_temp)
        psi_t_0 = np.append(psi_t_0, psi_t_0_temp)
        
elif ip.flag_type_init == 3: 
    m_period_tmp = 4
    
    alpha_arr_temp = cp.alpha_arr_isl[m_period_tmp]
    psi_arr_temp = cp.psi_arr_isl[m_period_tmp]
    
    alpha_0 = alpha_arr_temp.flatten(order='F')   
    psi_t_0 = psi_arr_temp.flatten(order='F')

# Line tracing 
R_line_tmp = np.array([])
Z_line_tmp = np.array([])
theta_line_tmp = np.array([])
theta_line_grad_tmp = np.array([])

R_lines = [[0 for _ in range(ip.n_lines_tot)] for _ in range(ip.n_poincare)]
Z_lines = [[0 for _ in range(ip.n_lines_tot)] for _ in range(ip.n_poincare)]
theta_lines = [[0 for _ in range(ip.n_lines_tot)] for _ in range(ip.n_poincare)] 

gamma_lines_sv = []
p_gamma_lines_sv = []
R_lines_sv = []
Z_lines_sv = []

gamma_lines_sn = []
p_gamma_lines_sn = []
R_lines_sn = []
Z_lines_sn = []

temp=5

for i_lines in range(ip.n_lines_tot):
    y1 = np.array([[psi_t_0[i_lines]],[alpha_0[i_lines]]])
    
    if ip.flag_analyt:
        # Pendulum solution
        y = pendulum(phi_span, y1)
#        if i_lines % 1000 == 0:
#            print(i_lines)
        
        for i_poincare in range(ip.n_poincare):
            #i_out = list(range(i_poincare*int(ip.points/ip.n_poincare),i_poincare*int(ip.points/ip.n_poincare)+len(phi_span)-ip.points,ip.points))
            i_out = list(range(i_poincare,i_poincare+1,1))
            R_line_tmp, Z_line_tmp, theta_line_tmp = \
            from_phase_to_real(phi_span, y, i_out)
            
            theta_line_grad_tmp = normalize_angle(theta_line_tmp) / np.pi * 180
            
            R_lines[i_poincare][i_lines]=R_line_tmp
            Z_lines[i_poincare][i_lines]=Z_line_tmp  
            theta_lines[i_poincare][i_lines]=theta_line_grad_tmp
            
    if ip.flag_Stoermer_Verlet:
        # Störmer-Verlet solution
        y_sv = num_solver(phi_span, y1, momentum_rhs, flag_sv = ip.flag_Stoermer_Verlet)
        
        gamma_lines_sv, p_gamma_lines_sv, R_lines_sv, Z_lines_sv = \
        from_phase_to_real(phi_span, y_sv, gamma_lines_sv, p_gamma_lines_sv,
                           R_lines_sv, Z_lines_sv)
    if ip.flag_scipy_num:
        y_sn = num_solver(phi_span, y1, momentum_rhs, flag_sn = ip.flag_scipy_num, method = ip.method)
        
        gamma_lines_sn, p_gamma_lines_sn, R_lines_sn, Z_lines_sn = \
        from_phase_to_real(phi_span, y_sn, gamma_lines_sn, p_gamma_lines_sn,
                           R_lines_sn, Z_lines_sn)

R_iota_res_test, Z_iota_res_test = psi_theta_to_R_Z(1.0+1/6, 0)