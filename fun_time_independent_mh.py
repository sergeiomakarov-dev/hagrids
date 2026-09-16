#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.special import ellipj, ellipk, ellipkinc
from scipy.integrate import solve_ivp
import input_params as ip
import mpmath
from scipy.optimize import fsolve
import magnetic_feild_calculation as mfc

# functions
def pendulum(t, y1):
    '''Classical pendulum solution'''
    
    # Set arbitrary precision (e.g., 100 decimal places)
    mpmath.mp.dps = 100

    
    # The initial pololdal angle gamma_1, the initial momentum p_1, energy E_tilde and corresponding k
    x_point = 0
    
    g = -ip.A
    psi_0 = (ip.iota_res - ip.iota_b) / ip.iota_a
    gamma_1_tot = ip.m * y1[1, 0]
    m_period = np.floor((gamma_1_tot + np.pi) / (2 * np.pi)) #floating presision problem
    gamma_1 = gamma_1_tot - 2 * np.pi * m_period
    
    t_norm = t * ip.m * np.sqrt(ip.iota_a)
    
    # Adjust for precision errors near -np.pi
    if -np.pi - 1e-15 < gamma_1 < -np.pi:
        gamma_1 = -np.pi
    
    
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
    u = (p_p_init * (t_norm - t_norm[0]) + t_phase) * np.sqrt(g);
    
    # Initialize arrays
    gamma = np.zeros(len(t_norm))
    gamma_period = np.zeros(len(t_norm))
    p = np.zeros(len(t_norm))
    yout = np.zeros((2, len(t_norm)))
    
    # The classical pendulum solution for gamma and p
    if np.isreal(xi_1):
        if k < 1:
            for i in range(len(t_norm)): 
 
                sn_u, cn_u, dn_u, _ = ellipj(np.float64(u[i]), np.float64(k))
                gamma[i] = 2 * np.arcsin(np.float64(k * sn_u))

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
            for i in range(len(t_norm)):
                sn_u, cn_u, dn_u, _ = ellipj(np.float64(k * u[i]), np.float64(1 / k))
                gamma_period[i] = 2 * np.arcsin(np.float64(sn_u))

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
                for i in range(len(t_norm)):   
                    gamma[i] = gamma_1
                    p[i] = p_1
                    if np.isreal(p[i]):
                        yout[0, i] = psi_from_momentum(p[i], ip.iota_a, psi_0)
                    else:                    
                        yout[0, i] = np.nan
                    yout[1, i] = alpha_from_gamma(gamma[i], m_period, ip.m)
            elif x_point == 0:
                for i in range(len(t_norm)):
                    sn_u, cn_u, dn_u, _ = ellipj(np.float64(k * u[i]), np.float64(1 / k))
                    gamma[i] = 2 * np.arcsin(np.float64(sn_u))
                    p[i] = p_p_init*momentum(E_tilde, g, gamma[i])
                    yout[0, i] = psi_from_momentum(p[i], ip.iota_a, psi_0)
                    yout[1, i] = alpha_from_gamma(gamma[i], m_period, ip.m)   
    else:
        yout[0, :] = np.nan
        yout[1, :] = np.nan
    return yout 

def sine_solution(t, y1):
    '''Classical pendulum solution near O-point'''

    
    # The initial pololdal angle gamma_1, the initial momentum p_1, energy E_tilde and corresponding k
    x_point = 0
    
    g = -ip.A
    psi_0 = (ip.iota_res - ip.iota_b) / ip.iota_a
    gamma_1_tot = ip.m * y1[1, 0]
    m_period = np.floor((gamma_1_tot + np.pi) / (2 * np.pi)) #floating presision problem
    gamma_1 = gamma_1_tot - 2 * np.pi * m_period
    
    t_norm = t * ip.m * np.sqrt(ip.iota_a)
    
    # Adjust for precision errors near -np.pi
    if -np.pi - 1e-15 < gamma_1 < -np.pi:
        gamma_1 = -np.pi
    
    
    if gamma_1 >= np.pi or gamma_1 < -np.pi:
        raise ValueError("gamma_1 >= np.pi or gamma_1 < -np.pi")
        
    p_1 = np.sqrt(ip.iota_a) * (y1[0, 0]-psi_0) 
    E_tilde = p_1 ** 2 / 2 + g * gamma_1 ** 2 / 2
    k = np.sqrt(E_tilde / (2 * g))
    
    xi_1 = gamma_1/2
    
    t_phase=np.arcsin(xi_1/ k ) / (np.sqrt(g))

    # The sign of the initial momentum      
    if p_1 > 0:
        p_p_init = 1
    else:
        p_p_init = -1
        
    gamma_max = np.sqrt(2 * E_tilde / g)
        
    uT = 4 * np.arcsin(gamma_max / (2 * k)) 
    
    T = np.sqrt(1 / g) * uT
        
    # Modified time array
    u = (p_p_init * (t_norm - t_norm[0]) + t_phase) * np.sqrt(g);
    
    # Initialize arrays
    gamma = np.zeros(len(t_norm))
    p = np.zeros(len(t_norm))
    yout = np.zeros((2, len(t_norm)))
    
    # The sine solution for gamma and p
    if np.isreal(xi_1):
        if k < 1:
            for i in range(len(t_norm)):   
                gamma[i] = 2 * k * np.sin(u[i])
                               
                
                if (np.mod(u[i], uT) <= uT / 4) or (np.mod(u[i], uT) >= 3 * uT / 4):
                    p[i] = p_p_init * np.sqrt(2 * (E_tilde - g * gamma[i] ** 2 / 2))
                else:
                    p[i] = -1 * p_p_init * np.sqrt(2 * (E_tilde - g * gamma[i] ** 2 / 2))

                if np.isreal(p[i]):
                    yout[0, i] = psi_from_momentum(p[i], ip.iota_a, psi_0)
                else:                    
                    yout[0, i] = np.nan
                yout[1, i] = alpha_from_gamma(gamma[i], m_period, ip.m) 
        else: 
            yout[0, :] = np.nan
            yout[1, :] = np.nan   
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
    x = r * np.sin(theta)
    y = r * np.cos(theta)
    return x, y

def yoshida_integrator(t, p, gamma):
    """
    Perform simulation using the 4th-order Yoshida integrator over a given time span.

    Parameters:
    - t: Array of time points.
    - gamma: Array to store gamma (analogous to position) at each time step.
    - p: Array to store momentum at each time step.
    - momentum_rhs: Function to compute the gradient of potential energy with respect to gamma (dV/dgamma).
    - coordinate_rhs: Function to compute the gradient of kinetic energy with respect to momentum (dT/dp).

    Returns:
    - gamma: Updated array of gamma values at each time step.
    - p: Updated array of momenta at each time step.
    """
    
    # Calculate the time step size from the time array
    dt = t[1] - t[0]
    
    # Yoshida integrator coefficients
    w0 = -2 ** (1 / 3) / (2 - 2 ** (1 / 3))
    w1 = 1 / (2 - 2 ** (1 / 3))
    c1 = w1 / 2
    c2 = (w0 + w1) / 2
    c3 = c2
    c4 = c1
    d1 = w1
    d2 = w0
    d3 = d1

    # Perform the integration over the time span
    for i in range(len(t) - 1):
        # Step 1: Update gamma with the first coefficient c1
        gamma_temp = gamma[i] + c1 * dt * coordinate_rhs(gamma[i], p[i])
        
        # Step 2: Update momentum p with the first coefficient d1
        p_temp = p[i] + d1 * dt * momentum_rhs(gamma_temp, p[i])
        
        # Step 3: Update gamma with the second coefficient c2
        gamma_temp += c2 * dt * coordinate_rhs(gamma_temp, p_temp)
        
        # Step 4: Update momentum p with the second coefficient d2
        p_temp += d2 * dt * momentum_rhs(gamma_temp, p_temp)
        
        # Step 5: Update gamma with the third coefficient c3
        gamma_temp += c3 * dt * coordinate_rhs(gamma_temp, p_temp)
        
        # Step 6: Update momentum p with the second coefficient d2
        p_temp += d3 * dt * momentum_rhs(gamma_temp, p_temp)
        
        # Step 7: Update gamma with the third coefficient c3
        gamma[i + 1] = gamma_temp + c4 * dt * coordinate_rhs(gamma_temp, p_temp)
        
        # Step 8: Update momentum p with the third coefficient d3
        p[i + 1] = p_temp

    return p, gamma

def dual_stoermer_verlet(t, p, gamma):
    dt = t[1] - t[0] 
    
    for i in range(len(t)-1):
        
        # First half step for momentum
        a = momentum_rhs(gamma[i], p[i])
        p_half = p[i] + 0.5 * a * dt
        
        # Full step for positions
        gamma[i + 1] = gamma[i] + coordinate_rhs(gamma[i], p_half) * dt
        
        # Second half step for momentum
        a_next = momentum_rhs(gamma[i + 1], p_half)
        p[i + 1] = p_half + 0.5 * a_next * dt
    return p, gamma

def dual_stoermer_verlet_mod(t, psi, alpha):
    dt = t[1] - t[0] 
    
    for i in range(len(t)-1):
        
        # First half step for momentum
        a = momentum_rhs_mod(alpha[i], t[i])
        psi_half = psi[i] + 0.5 * a * dt
        
        # Full step for positions
        alpha[i + 1] = alpha[i] + coordinate_rhs_mod(psi_half) * dt 
        
        # Second half step for momentum
        a_next = momentum_rhs_mod(alpha[i + 1], t[i + 1])
        psi[i + 1] = psi_half + 0.5 * a_next * dt
        
    return psi, alpha
    
def num_solver(t, y1, flag_sv = 0, flag_sn = 0, flag_sv_mod = 0, flag_yo = 0, method = 'RK45'):
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
    t_norm = t * ip.m * np.sqrt(ip.iota_a)
    
    alpha = np.zeros(len(t))

    psi = np.zeros(len(t))
    
    alpha[0] = y1[1, 0]    
    psi[0] = y1[0, 0] 

    if flag_sv:
        p, gamma = dual_stoermer_verlet(t_norm, p, gamma)
        
    if flag_sn:
        y_temp = solve_ivp(hamiltonian_pendulum, [t_norm[0], t_norm[-1]], [p[0], gamma[0]], args=(ip.A,), t_eval=t_norm, method=method)
        p, gamma = y_temp.y
        
    if flag_yo:
        p, gamma = yoshida_integrator(t_norm, p, gamma)
        
    if flag_sv_mod :
        psi, alpha = dual_stoermer_verlet_mod(t, psi, alpha)
    
    if flag_sv or flag_sn or flag_yo:
        yout[0, :] = psi_from_momentum(p[:], ip.iota_a, psi_0)
        yout[1, :] = alpha_from_gamma_mod(gamma[:], ip.m)   
    elif flag_sv_mod :
        yout[0, :] = psi[:]
        yout[1, :] = alpha[:]  

    
    return yout 

def momentum_rhs(x, p):
    g = -ip.A
    p_0 = np.sqrt(ip.iota_a) * (ip.psi_0) 
    if ip.flag_type_A_fun == 0:
        mom_rhs = -g * np.sin(x)
    elif ip.flag_type_A_fun == 1:
        mom_rhs = -g * np.heaviside((1 + ip.k * p / p_0), 1) * (1 + ip.k * p / p_0) * np.sin(x)
    else:
        raise ValueError("Invalid value for flag_type_A_fun. Must be 0 or 1.")
    return mom_rhs

def coordinate_rhs(x, p):
    g = -ip.A
    p_0 = np.sqrt(ip.iota_a) * (ip.psi_0) 
    if ip.flag_type_A_fun == 0:
        coord_rhs = p
    elif ip.flag_type_A_fun == 1:
        coord_rhs = p - np.heaviside((1 + ip.k * p / p_0), 1) * g * ip.k / p_0 * np.cos(x)
        if (1 + ip.k * p / p_0) == 0:
            raise ValueError("The (1 + ip.k * p / p_0) = 0 point is non-differentiable")
    else:
        raise ValueError("Invalid value for flag_type_A_fun. Must be 0 or 1.")
    return coord_rhs

def momentum_rhs_mod(x,t):
    return ip.m * ip.A * np.sin(ip.m * x)

def coordinate_rhs_mod(p):
    return ip.iota_a * p + ip.iota_b - ip.iota_res

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

# Function to solve theta = vartheta - varepsilon * sin(vartheta) for a given theta and varepsilon
def solve_for_vartheta(theta_value, varepsilon_value):
    # Define the equation to solve
    def equation(vartheta, varepsilon, theta):
        return vartheta - varepsilon * np.sin(vartheta) - theta
    
    # Initial guess for vartheta
    vartheta_initial_guess = theta_value
    
    # Solve the equation numerically
    vartheta_solution = fsolve(equation, vartheta_initial_guess, args=(varepsilon_value, theta_value))
    
    return vartheta_solution

def psi_theta_to_R_Z(psi,theta, flag_vartheta = 1):
    r = np.sqrt(psi / (ip.B_0 * np.pi))
    varepsilon = r / ip.R_0
    # vartheta
    if flag_vartheta:
        vartheta = solve_for_vartheta(theta, varepsilon)
    else:
        vartheta = theta
    
    # In cylindrical coordinates
    R = ip.R_0 + r * np.cos(vartheta)
    Z = r * np.sin(vartheta)
    return R, Z

def R_Z_to_psi_theta(R,Z, flag_vartheta = 1):
    # In toroidal coordinates
    r = np.hypot(R - ip.R_0, Z)
    vartheta = np.arctan2(Z, R - ip.R_0)
    varepsilon = r / ip.R_0
    # vartheta
    if flag_vartheta:
        theta = vartheta - varepsilon * np.sin(vartheta)
    else:
        theta = vartheta
    
    psi = r**2 * (ip.B_0 * np.pi)
    return psi, theta

def R_Z_to_r_vartheta(R,Z):
    # In toroidal coordinates
    r = np.hypot(R - ip.R_0, Z)
    vartheta = np.arctan2(Z, R - ip.R_0)
    return r, vartheta

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

    R, Z = psi_theta_to_R_Z(psi,theta, flag_vartheta = ip.flag_vartheta)
    
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

def transform_theta_to_vartheta(R_in, Z_in):
    psi_temp, theta_temp = R_Z_to_psi_theta(R_in,Z_in, flag_vartheta = 0)
    R_out, Z_out = psi_theta_to_R_Z(psi_temp, theta_temp, flag_vartheta = 1)
    return R_out, Z_out

def make_poincare (phi_span, alpha_0, psi_t_0, n_lines_tot):
    R_line_tmp = np.array([])
    Z_line_tmp = np.array([])
    theta_line_tmp = np.array([])
    theta_line_grad_tmp = np.array([])
    
    R_lines = [[0 for _ in range(n_lines_tot)] for _ in range(ip.n_poincare)]
    Z_lines = [[0 for _ in range(n_lines_tot)] for _ in range(ip.n_poincare)]
    theta_lines = [[0 for _ in range(n_lines_tot)] for _ in range(ip.n_poincare)] 
    phi = [0 for _ in range(ip.n_poincare)] 
    
    for i_lines in range(n_lines_tot):
        y1 = np.array([[psi_t_0[i_lines]],[alpha_0[i_lines]]])
        
        if ip.flag_analyt:
            # Pendulum solution
            y = pendulum(phi_span, y1)
            
        if ip.flag_analyt_approx:    
            y = sine_solution(phi_span, y1)
        
        if ip.flag_Stoermer_Verlet:
            # Pendulum solution
            y = num_solver(phi_span, y1, flag_sv = ip.flag_Stoermer_Verlet)

        if ip.flag_Stoermer_Verlet_mod:
            # Pendulum solution
            y = num_solver(phi_span, y1, flag_sv_mod = ip.flag_Stoermer_Verlet_mod)
        
        if ip.flag_Yoshida:
            # Pendulum solution
            y = num_solver(phi_span, y1, flag_yo = ip.flag_Yoshida)
            
        for i_poincare in range(ip.n_poincare):
            #i_out = list(range(i_poincare*int(ip.points/ip.n_poincare),i_poincare*int(ip.points/ip.n_poincare)+len(phi_span)-ip.points,ip.points))
            #i_out = list(range(i_poincare,i_poincare+1,1))
            i_out = list(range(ip.n_space_points * i_poincare,ip.n_space_points * i_poincare+1,1))
            R_line_tmp, Z_line_tmp, theta_line_tmp = \
            from_phase_to_real(phi_span, y, i_out)
                
            theta_line_grad_tmp = normalize_angle(theta_line_tmp) / np.pi * 180
                
            R_lines[i_poincare][i_lines]=R_line_tmp
            Z_lines[i_poincare][i_lines]=Z_line_tmp  
            theta_lines[i_poincare][i_lines]=theta_line_grad_tmp

    for i_poincare in range(ip.n_poincare):
        i_out = list(range(ip.n_space_points * i_poincare,ip.n_space_points * i_poincare+1,1))
        phi[i_poincare] = phi_span[i_out]
               
    return R_lines, Z_lines, theta_lines, phi

def cubic_from_derivatives_monotonic(p, q, x_vals):
    """
    Return y(x_vals), where y(x) is the cubic polynomial mapping [0,1] to [0,1]
    with y'(0) = p and y'(1) = q. p and q are the cell sizes at the two ends of
    the interval relative to the uniform spacing. Raises ValueError if the
    polynomial is not monotonic on [0,1].
    """
    a = q + p - 2
    b = 3 - q - 2 * p
    c = p

    def y(x):
        return a * x**3 + b * x**2 + c * x

    def dy_dx(x):
        return 3 * a * x**2 + 2 * b * x + c

    # Check monotonicity at many points
    derivative_vals = dy_dx(x_vals)
    if np.any(derivative_vals < 0):
        raise ValueError(f"The cubic with y'(0)={p}, y'(1)={q} is not monotonic on [0,1].")
    y_out = y(x_vals)
    return y_out

def poloidal_spacing(a, b, n, alpha=1, beta=1, flag_spacing_fun=0):
    # Generate linearly spaced points between 0 and 1
    u = np.linspace(0, 1, n)

    if flag_spacing_fun==0:
        # Transform linearly spaced points using the power function
        output_points = a + (b - a) * u ** alpha
    elif flag_spacing_fun==1:
        # Transform linearly spaced points using the cubic polynomial with the slopes alpha and beta at the ends
        tmp = cubic_from_derivatives_monotonic(alpha, beta, u)
        output_points = a + (b - a) * tmp
    return output_points

def radial_spacing(a, b, n, alpha=1, beta=1, flag_spacing_fun=0):
    # Generate linearly spaced points between 0 and 1
    u = np.linspace(0, 1, n)

    if flag_spacing_fun==0:
        # Transform linearly spaced points using the power function
        output_points = a + (b - a) * u ** alpha
    elif flag_spacing_fun==1:
        # Transform linearly spaced points using the cubic polynomial with the slopes alpha and beta at the ends
        tmp = cubic_from_derivatives_monotonic(alpha, beta, u)
        output_points = a + (b - a) * tmp
    return output_points

def calculate_hypotenuse(x, y):
    return np.sqrt(x**2 + y**2)

def calculate_hypotenuse_3d(x, y, z):
    return np.sqrt(x**2 + y**2 + z**2)

def cell_size(R, Z):
    length = 0
    size = np.zeros(len(R)-1)
    length = np.zeros(len(R)-1)
    i=0
    size[i] = calculate_hypotenuse(R[i+1]-R[i], Z[i+1]-Z[i])
    length[i] = size[i] / 2 
    for i in range(1,len(R)-1):
        size[i] = calculate_hypotenuse(R[i+1]-R[i], Z[i+1]-Z[i])
        length[i] = length[i-1] + size[i-1] / 2 + size[i] / 2
    length_tot = length[-1] + size[-1] / 2
    return size, length, length_tot

def generate_points(R1, Z1, R2, Z2, number_of_points, alpha=1):
    # Generate equally spaced points for R and Z
    Rn = radial_spacing(R1, R2, number_of_points, alpha=alpha)
    Zn = radial_spacing(Z1, Z2, number_of_points, alpha=alpha)
    return Rn, Zn

def structure_R_Z_arrays(R_lines, Z_lines, cols):
    R_lines_theta_r = []
    Z_lines_theta_r = []
    
    for i_poincare in range(ip.n_poincare):
        R_lines_temp=R_lines[i_poincare]
        Z_lines_temp=Z_lines[i_poincare]
        
        R_lines_temp_arr=np.array(R_lines_temp)
        Z_lines_temp_arr=np.array(Z_lines_temp)
        
        #cols = ip.n_E_core+ip.n_E+ip.n_E_pfr-2
        rows = int(len(R_lines_temp) / cols)
        
        R_lines_theta_r.append(R_lines_temp_arr.reshape(rows ,cols))
        Z_lines_theta_r.append(Z_lines_temp_arr.reshape(rows ,cols))
    return R_lines_theta_r, Z_lines_theta_r

def comb_pos_neg(R_lines_theta_r_pos, Z_lines_theta_r_pos, \
                 phi_pos,\
                 R_lines_theta_r_neg, Z_lines_theta_r_neg,\
                 phi_neg):
    R_lines_theta_r = [0] * (2 * len(R_lines_theta_r_pos) - 1)
    Z_lines_theta_r = [0] * (2 * len(R_lines_theta_r_pos) - 1)
    phi_lines_theta_r = [0] * (2 * len(R_lines_theta_r_pos) - 1)
    
    R_lines_theta_r[len(R_lines_theta_r_pos)-1:2 * len(R_lines_theta_r_pos) - 1] = \
    R_lines_theta_r_pos[0:]
    
    Z_lines_theta_r[len(R_lines_theta_r_pos)-1:2 * len(R_lines_theta_r_pos) - 1] = \
    Z_lines_theta_r_pos[0:]
    
    phi_lines_theta_r[len(R_lines_theta_r_pos)-1:2 * len(R_lines_theta_r_pos) - 1] = \
    phi_pos[0:]    
    
    R_lines_theta_r[0:len(R_lines_theta_r_pos)-1] = \
    R_lines_theta_r_neg[:0:-1]
    
    Z_lines_theta_r[0:len(R_lines_theta_r_pos)-1] = \
    Z_lines_theta_r_neg[:0:-1]

    phi_lines_theta_r[0:len(R_lines_theta_r_pos)-1] = \
    phi_neg[:0:-1]
    
    return R_lines_theta_r , Z_lines_theta_r, phi_lines_theta_r

def A_fun(p, p_0, k, flag_type_A_fun):
    if flag_type_A_fun == 0:
        f = 1.0
    elif flag_type_A_fun == 1:
        f = (1.0 + k * p / p_0)
    else:
        raise ValueError("Invalid value for flag_type_A_fun. Must be 0 or 1.")
    
    return f

def zero_discriminant(a, b):
    c = b ** 2 / (4 * a)
    return c

def separatrix_energy(p_0, k, g, flag_type_A_fun):
    if flag_type_A_fun == 0:
        c = zero_discriminant(0.5, 0) 
    elif flag_type_A_fun == 1:
        c = zero_discriminant(0.5, 2 * g * k / p_0) 
    else:
        raise ValueError("Invalid value for flag_type_A_fun. Must be 0 or 1.")  
    E_sep = - c + 2 * g
    return E_sep

def X_O_energy(a, b, c, flag_type_A_fun):
    if flag_type_A_fun == 0:
        energy = - solve_discriminant(a, 0.0, c) / (4 * a)
    elif flag_type_A_fun == 1:
        energy = - solve_discriminant(a, b, c) / (4 * a)
    else:
        raise ValueError("Invalid value for flag_type_A_fun. Must be 0 or 1.")
    return energy

solve_discriminant = lambda a, b, c: b**2 - 4 * a * c

def solve_quadratic(a, b, c):
    # Calculate the discriminant
    discriminant = solve_discriminant(a, b, c)
    
    if discriminant < 0:
        raise ValueError(f"Discriminant = {discriminant} < 0")
        
    # Calculate the two solutions
    root1 = (-b + np.sqrt(discriminant)) / (2 * a)
    root2 = (-b - np.sqrt(discriminant)) / (2 * a)
    
    return root1, root2

def solut_algeb_eq(a, b, c, flag_type_A_fun):
    if flag_type_A_fun == 0:
        sol1, sol2 = solve_quadratic(a, 0.0, c) 
    elif flag_type_A_fun == 1:
        sol1, sol2 = solve_quadratic(a, b, c) 
    else:
        raise ValueError("Invalid value for flag_type_A_fun. Must be 0 or 1.")
    
    return sol1, sol2

def solve_quadratic_same_roots(a, b, c):
    
    # Calculate the two solutions
    root1 = -b / (2 * a)
    root2 = -b / (2 * a)
    
    return root1, root2

def solut_same_roots(a, b, c, flag_type_A_fun):
    if flag_type_A_fun == 0:
        sol1, sol2 = solve_quadratic_same_roots(a, 0.0, c) 
    elif flag_type_A_fun == 1:
        sol1, sol2 = solve_quadratic_same_roots(a, b, c) 
    else:
        raise ValueError("Invalid value for flag_type_A_fun. Must be 0 or 1.")
    
    return sol1, sol2

def check_discriminant(a, b, c, flag_type_A_fun):
    if flag_type_A_fun == 0:
        discriminant = solve_discriminant(a, 0.0, c) 
    elif flag_type_A_fun == 1:
        discriminant = solve_discriminant(a, b, c) 
    else:
        raise ValueError("Invalid value for flag_type_A_fun. Must be 0 or 1.")
    
    return discriminant

def mag_field_in_pol_planes(R_lines_theta_r, Z_lines_theta_r, phi_theta_r, flag_vartheta=ip.flag_vartheta):
    Br_lines_theta_r = []
    Bvartheta_lines_theta_r = []
    Bvarphi_lines_theta_r = []
    BR_lines_theta_r = []
    BZ_lines_theta_r = []
    BY_lines_theta_r = []
    BX_lines_theta_r = []
    B_lines_theta_r = []
    for i_poincare in range(len(R_lines_theta_r)): 
        R_lines_temp=R_lines_theta_r[i_poincare]
        Z_lines_temp=Z_lines_theta_r[i_poincare]
        phi_temp=phi_theta_r[i_poincare][0]
        
        Br_lines_theta_r.append(np.zeros(R_lines_temp.shape))
        Bvartheta_lines_theta_r.append(np.zeros(R_lines_temp.shape))
        Bvarphi_lines_theta_r.append(np.zeros(R_lines_temp.shape))
        BR_lines_theta_r.append(np.zeros(R_lines_temp.shape))
        BZ_lines_theta_r.append(np.zeros(R_lines_temp.shape))
        BX_lines_theta_r.append(np.zeros(R_lines_temp.shape))
        BY_lines_theta_r.append(np.zeros(R_lines_temp.shape))
        B_lines_theta_r.append(np.zeros(R_lines_temp.shape))
        
        for i_theta in range(R_lines_temp.shape[0]):
            for i_r in range(R_lines_temp.shape[1]):
                r_temp, vartheta_temp = R_Z_to_r_vartheta(R_lines_temp[i_theta,i_r],Z_lines_temp[i_theta,i_r])
                Br_temp, Bvartheta_temp, Bvarphi_temp = \
                    mfc.magnetic_field_components(r_temp, vartheta_temp, phi_temp, ip.B_0, ip.R_0,
                                                  iota_a=ip.iota_a, iota_b=ip.iota_b,
                                                  m=ip.m, m1=ip.m1, m2=ip.m2, m3=ip.m3,
                                                  iota_res=ip.iota_res, iota_res1=ip.iota_res1, iota_res2=ip.iota_res2, iota_res3=ip.iota_res3,
                                                  A=ip.A, A1=ip.A1, A2=ip.A2, A3=ip.A3, k=ip.k,
                                                  flag_form=ip.flag_Stoermer_Verlet_td or ip.flag_scipy_num_td or ip.flag_Yoshida_td,
                                                  flag_toroildal=flag_vartheta,
                                                  flag_type_A_fun=ip.flag_type_A_fun)
                
                B_temp = calculate_hypotenuse_3d(Br_temp, Bvartheta_temp, Bvarphi_temp)
                
                BR_temp = Br_temp * np.cos(vartheta_temp) -\
                    Bvartheta_temp * np.sin(vartheta_temp) 
                BZ_temp = Br_temp * np.sin(vartheta_temp) +\
                    Bvartheta_temp * np.cos(vartheta_temp)    
                    
                BX_temp = BR_temp * np.sin(phi_temp) +\
                    Bvarphi_temp * np.cos(phi_temp) 
                BY_temp = BR_temp * np.cos(phi_temp) -\
                    Bvarphi_temp * np.sin(phi_temp)
                    
                Br_lines_theta_r[i_poincare][i_theta,i_r]=Br_temp
                Bvartheta_lines_theta_r[i_poincare][i_theta,i_r]=Bvartheta_temp
                Bvarphi_lines_theta_r[i_poincare][i_theta,i_r]=Bvarphi_temp
                BR_lines_theta_r[i_poincare][i_theta,i_r]=BR_temp
                BZ_lines_theta_r[i_poincare][i_theta,i_r]=BZ_temp
                BX_lines_theta_r[i_poincare][i_theta,i_r]=BX_temp
                BY_lines_theta_r[i_poincare][i_theta,i_r]=BY_temp
                B_lines_theta_r[i_poincare][i_theta,i_r]=B_temp
    return B_lines_theta_r, BX_lines_theta_r, BY_lines_theta_r, BZ_lines_theta_r, BR_lines_theta_r, Bvarphi_lines_theta_r, Br_lines_theta_r, Bvartheta_lines_theta_r
