#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.special import ellipj, ellipk, ellipkinc
from scipy.integrate import solve_ivp
import input_params as ip
import fun_time_independent_mh as fi


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
        gamma_temp = gamma[i] + c1 * dt * coordinate_rhs(gamma[i], p[i], t[i])
        t_temp = t[i] + c1 * dt
        
        # Step 2: Update momentum p with the first coefficient d1
        p_temp = p[i] + d1 * dt * momentum_rhs(gamma_temp, p[i], t_temp)
        
        # Step 3: Update gamma with the second coefficient c2
        gamma_temp += c2 * dt * coordinate_rhs(gamma_temp, p_temp, t_temp)
        t_temp += c2 * dt
        
        # Step 4: Update momentum p with the second coefficient d2
        p_temp += d2 * dt * momentum_rhs(gamma_temp, p_temp, t_temp)
        
        # Step 5: Update gamma with the third coefficient c3
        gamma_temp += c3 * dt * coordinate_rhs(gamma_temp, p_temp, t_temp)
        t_temp += c3 * dt
        
        # Step 6: Update momentum p with the second coefficient d2
        p_temp += d3 * dt * momentum_rhs(gamma_temp, p_temp, t_temp)
        
        # Step 7: Update gamma with the third coefficient c3
        gamma[i + 1] = gamma_temp + c4 * dt * coordinate_rhs(gamma_temp, p_temp, t_temp)
        
        # Step 8: Update momentum p with the third coefficient d3
        p[i + 1] = p_temp

    return p, gamma

def dual_stoermer_verlet(t, p, gamma):
    dt = t[1] - t[0] 
    
    for i in range(len(t)-1):
        
        # First half step for momentum
        a = momentum_rhs(p[i], gamma[i], t[i])
        p_half = p[i] + 0.5 * a * dt
        
        # Full step for positions
        gamma[i + 1] = gamma[i] + coordinate_rhs(gamma[i], p_half, t[i]) * dt 
        
        # Second half step for momentum
        a_next = momentum_rhs(gamma[i + 1], p_half, t[i + 1])
        p[i + 1] = p_half + 0.5 * a_next * dt
        
    return p, gamma
    
def num_solver(t, y1, flag_sv = 0, flag_sn = 0, flag_yo = 0, method = 'RK45'):
    """
    Perform dual Störmer-Verlet integration for a given number of steps.
    
    """
    
    # Initialize arrays
    yout = np.zeros((2, len(t)))
    y_wrk = np.zeros((2, len(t)))
    y_wrk[0, 0] = y1[0, 0]
    y_wrk[1, 0] = y1[1, 0]
    
    if flag_sv:
        tmp_p, tmp_gamma = dual_stoermer_verlet(t, y_wrk[0, :], y_wrk[1, :]) 
        y_wrk[0, :] = tmp_p
        y_wrk[1, :] = tmp_gamma
        
    if flag_sn:
        y_temp = solve_ivp(hamiltonian_time_dependent, [t[0], t[-1]], \
                           [y1[0, 0], y1[1, 0]], \
                               args=(ip.A1, ip.A2, ip.A3, ip.iota_a, ip.iota_b, 
                                     ip.m1, ip.m2, ip.m3, ip.iota_res1, ip.iota_res2, ip.iota_res3), \
                                   t_eval=t, method=method)
        y_wrk = y_temp.y
        
    if flag_yo:
        tmp_p, tmp_gamma = yoshida_integrator(t, y_wrk[0, :], y_wrk[1, :])
        y_wrk[0, :] = tmp_p
        y_wrk[1, :] = tmp_gamma
        
    yout[0, :] = y_wrk[0, :] 
    yout[1, :] = y_wrk[1, :]
        
    return yout 

def tanh_heaviside_with_width_and_shift(x, width=1, shift=0):
    return 0.5 * (1 + np.tanh(2 * (x - shift) / width))

def tanh_heaviside_derivative(x, width=1, shift=0):
    u = 2 * (x - shift) / width
    return (1 / width) * (1 - np.tanh(u)**2)

def anti_tanh_heaviside_with_width_and_shift(x, width=1, shift=0):
    return 0.5 * (1 + np.tanh(-2 * (x - shift) / width))

def anti_tanh_heaviside_derivative(x, width=1, shift=0):
    u = -2 * (x - shift) / width
    return -(1 / width) * (1 - np.tanh(u)**2)

def Q(x, p, t):
    Q_out = ((ip.iota_a * (p - ip.psi_0) ** 2)/2 + ip.A * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) *
         np.cos(ip.m1 * (x - ip.iota_res1 * t)))
    return Q_out

def d_x_Q(x, p, t):
    Q_out = ip.iota_a * (p - ip.psi_0) + ip.k * ip.A / ip.psi_0 * np.cos(ip.m1 * (x - ip.iota_res1 * t))
    return Q_out

def d_p_Q(x, p, t):
    Q_out = (-ip.m * ip.A * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) *
         np.sin(ip.m1 * (x - ip.iota_res1 * t)))
    return Q_out

def f_modificator_2(x, p, t):
    f = tanh_heaviside_with_width_and_shift(Q(x, p, t), width=ip.Q_width_2, shift=ip.Q_shift_2)
    return f

def d_x_f_modificator_2(x, p, t):
    f = (tanh_heaviside_derivative(Q(x, p, t), width=ip.Q_width_2, shift=ip.Q_shift_2) *
         d_x_Q(x, p, t))
    return f

def d_p_f_modificator_2(x, p, t):
    f = (tanh_heaviside_derivative(Q(x, p, t), width=ip.Q_width_2, shift=ip.Q_shift_2) *
         d_p_Q(x, p, t))
    return f

def f_modificator(x, p, t):
    f = (anti_tanh_heaviside_with_width_and_shift(p, width=ip.psi_width_2, shift=ip.psi_shift_2) *
         f_modificator_2(x, p, t))
    return f

def d_x_f_modificator(x, p, t):
    f = (anti_tanh_heaviside_with_width_and_shift(p, width=ip.psi_width_2, shift=ip.psi_shift_2) *
         d_x_f_modificator_2(x, p, t))
    return f

def d_p_f_modificator(x, p, t):
    f = (anti_tanh_heaviside_derivative(p, width=ip.psi_width_2, shift=ip.psi_shift_2) *
         f_modificator_2(x, p, t) +
         anti_tanh_heaviside_with_width_and_shift(p, width=ip.psi_width_2, shift=ip.psi_shift_2) *
         d_p_f_modificator_2(x, p, t))
    return f

def momentum_rhs(x, p, t):
    if ip.flag_type_A_fun == 0:
        mom_rhs = (ip.m1 * ip.A1 * np.sin(ip.m1 * (x - ip.iota_res1 * t)) +  
                   ip.m2 * ip.A2 * np.sin(ip.m2 * (x - ip.iota_res2 * t)) +  
                   ip.m3 * ip.A3 * np.sin(ip.m3 * (x - ip.iota_res3 * t)))
    elif ip.flag_type_A_fun == 1:
        if ip.flag_f_modificator == 0:
            mom_rhs = ( np.heaviside(1 + ip.k * ((p - ip.psi_0) / ip.psi_0), 1) * 
                       (ip.m1 * ip.A1 * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) *
                       np.sin(ip.m1 * (x - ip.iota_res1 * t)) +  
                       ip.m2 * ip.A2 * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) * 
                       np.sin(ip.m2 * (x - ip.iota_res2 * t)) +                    
                       ip.m3 * ip.A3 * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) * 
                       np.sin(ip.m3 * (x - ip.iota_res3 * t))))
        elif ip.flag_f_modificator == 1:
            mom_rhs = ( np.heaviside(1 + ip.k * ((p - ip.psi_0) / ip.psi_0), 1) * 
                       (ip.m1 * ip.A1 * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) *
                       np.sin(ip.m1 * (x - ip.iota_res1 * t)) +  
                       ip.m2 * ip.A2 * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) * 
                       f_modificator(x, p, t) *
                       np.sin(ip.m2 * (x - ip.iota_res2 * t)) + 
                       ip.A2 * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) * 
                       d_x_f_modificator(x, p, t) *
                       np.cos(ip.m2 * (x - ip.iota_res2 * t)) +                   
                       ip.m3 * ip.A3 * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) * 
                       np.sin(ip.m3 * (x - ip.iota_res3 * t))))
        elif ip.flag_f_modificator == 2:
            mom_rhs = ( np.heaviside(1 + ip.k * ((p - ip.psi_0) / ip.psi_0), 1) * 
                       (ip.m1 * ip.A1 * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) *
                       np.sin(ip.m1 * (x - ip.iota_res1 * t))))
            for i_p in range(ip.n_perturbation):
                mom_rhs += ( np.heaviside(1 + ip.k * ((p - ip.psi_0) / ip.psi_0), 1) * 
                    (ip.m_arr[i_p] * ip.A_arr[i_p] * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) * 
                    f_modificator(x, p, t) *
                    np.sin(ip.m_arr[i_p] * (x - ip.iota_res_arr[i_p] * t)) + 
                    ip.A_arr[i_p] * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) * 
                    d_x_f_modificator(x, p, t) *
                    np.cos(ip.m_arr[i_p] * (x - ip.iota_res_arr[i_p] * t))))
        elif ip.flag_f_modificator == 3:
            # Perturbation 1 plus the additional perturbations as they are, without the envelope
            mom_rhs = ( np.heaviside(1 + ip.k * ((p - ip.psi_0) / ip.psi_0), 1) * 
                       (ip.m1 * ip.A1 * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) *
                       np.sin(ip.m1 * (x - ip.iota_res1 * t))))
            for i_p in range(ip.n_perturbation):
                mom_rhs += ( np.heaviside(1 + ip.k * ((p - ip.psi_0) / ip.psi_0), 1) * 
                    ip.m_arr[i_p] * ip.A_arr[i_p] * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) * 
                    np.sin(ip.m_arr[i_p] * (x - ip.iota_res_arr[i_p] * t)))
        else:
            raise ValueError("Invalid value for flag_f_modificator. Must be 0, 1, 2 or 3.")
    else:
        raise ValueError("Invalid value for flag_type_A_fun. Must be 0 or 1.")
    return mom_rhs

def coordinate_rhs(x, p, t):
    if ip.flag_type_A_fun == 0:
        coord_rhs = ip.iota_a * p + ip.iota_b
    elif ip.flag_type_A_fun == 1:
        if ip.flag_f_modificator == 0:
            coord_rhs = (ip.iota_a * p + ip.iota_b +
                         np.heaviside(1 + ip.k * ((p - ip.psi_0) / ip.psi_0), 1) *
                         (ip.A1 * (ip.k / ip.psi_0) *
                         np.cos(ip.m1 * (x - ip.iota_res1 * t)) +
                         ip.A2 * (ip.k / ip.psi_0) * 
                         np.cos(ip.m2 * (x - ip.iota_res2 * t)) +
                         ip.A3 * (ip.k / ip.psi_0) * 
                         np.cos(ip.m3 * (x - ip.iota_res3 * t))))
            if (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) == 0:
                raise ValueError("The (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) = 0 point is non-differentiable")
        elif ip.flag_f_modificator == 1:
            coord_rhs = (ip.iota_a * p + ip.iota_b +
                         np.heaviside(1 + ip.k * ((p - ip.psi_0) / ip.psi_0), 1) *
                         (ip.A1 * (ip.k / ip.psi_0) *
                         np.cos(ip.m1 * (x - ip.iota_res1 * t)) +
                         ip.A2 * (ip.k / ip.psi_0) * 
                         f_modificator(x, p, t) *
                         np.cos(ip.m2 * (x - ip.iota_res2 * t)) +
                         ip.A2 * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) *
                         d_p_f_modificator(x, p, t) *
                         np.cos(ip.m2 * (x - ip.iota_res2 * t)) +
                         ip.A3 * (ip.k / ip.psi_0) * 
                         np.cos(ip.m3 * (x - ip.iota_res3 * t))))
            if (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) == 0:
                raise ValueError("The (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) = 0 point is non-differentiable")
        elif ip.flag_f_modificator == 2:
            coord_rhs = (ip.iota_a * p + ip.iota_b +
                         np.heaviside(1 + ip.k * ((p - ip.psi_0) / ip.psi_0), 1) *
                         (ip.A1 * (ip.k / ip.psi_0) *
                         np.cos(ip.m1 * (x - ip.iota_res1 * t))))
            for i_p in range(ip.n_perturbation):
                coord_rhs += (np.heaviside(1 + ip.k * ((p - ip.psi_0) / ip.psi_0), 1) *
                    (ip.A_arr[i_p] * (ip.k / ip.psi_0) * 
                    f_modificator(x, p, t) *
                    np.cos(ip.m_arr[i_p] * (x - ip.iota_res_arr[i_p] * t)) +
                    ip.A_arr[i_p] * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) *
                    d_p_f_modificator(x, p, t) *
                    np.cos(ip.m_arr[i_p] * (x - ip.iota_res_arr[i_p] * t))))
            if (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) == 0:
                raise ValueError("The (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) = 0 point is non-differentiable")
        elif ip.flag_f_modificator == 3:
            # Perturbation 1 plus the additional perturbations as they are, without the envelope
            coord_rhs = (ip.iota_a * p + ip.iota_b +
                         np.heaviside(1 + ip.k * ((p - ip.psi_0) / ip.psi_0), 1) *
                         (ip.A1 * (ip.k / ip.psi_0) *
                         np.cos(ip.m1 * (x - ip.iota_res1 * t))))
            for i_p in range(ip.n_perturbation):
                coord_rhs += (np.heaviside(1 + ip.k * ((p - ip.psi_0) / ip.psi_0), 1) *
                    ip.A_arr[i_p] * (ip.k / ip.psi_0) * 
                    np.cos(ip.m_arr[i_p] * (x - ip.iota_res_arr[i_p] * t)))
            if (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) == 0:
                raise ValueError("The (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) = 0 point is non-differentiable")
        else:
            raise ValueError("Invalid value for flag_f_modificator. Must be 0, 1, 2 or 3.")
    else:
        raise ValueError("Invalid value for flag_type_A_fun. Must be 0 or 1.")
    return coord_rhs

def hamiltonian_time_dependent(t, z, A1, A2, A3, iota_a, iota_b, m1, m2, m3, iota_res1, iota_res2, iota_res3):
    p, q = z
    if ip.flag_type_A_fun == 0:
        dpdt = (m1 * A1 * np.sin(m1 * (q - iota_res1 * t)) + 
                m2 * A2 * np.sin(m2 * (q - iota_res2 * t)) + 
                m3 * A3 * np.sin(m2 * (q - iota_res3 * t)))
        dqdt = iota_a * p + iota_b
    elif ip.flag_type_A_fun == 1:
        dpdt = ( np.heaviside(1 + ip.k * ((p - ip.psi_0) / ip.psi_0), 1) * 
                (m1 * A1 * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) * 
                np.sin(m1 * (q - iota_res1 * t)) + 
                m2 * A2 * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) * 
                np.sin(m2 * (q - iota_res2 * t)) + 
                m3 * A3 * (1 + ip.k * ((p - ip.psi_0) / ip.psi_0)) * 
                np.sin(m3 * (q - iota_res3 * t))))
        dqdt = (ip.iota_a * p + ip.iota_b +
                np.heaviside(1 + ip.k * ((p - ip.psi_0) / ip.psi_0), 1) *
                (ip.A1 * (ip.k / ip.psi_0) *
                np.cos(ip.m1 * (q - ip.iota_res1 * t)) +
                ip.A2 * (ip.k / ip.psi_0) * 
                np.cos(ip.m2 * (q - ip.iota_res2 * t)) +
                ip.A3 * (ip.k / ip.psi_0) * 
                np.cos(ip.m3 * (q - ip.iota_res3 * t))))
    return [dpdt, dqdt]

momentum = lambda E_tilde_loc, g_loc, gamma_loc : np.sqrt(2 * (E_tilde_loc - g_loc * (1 - np.cos(gamma_loc))))           
psi_from_momentum = lambda p_loc, iota_a_loc, psi_0_loc : p_loc / np.sqrt(iota_a_loc) + psi_0_loc   
alpha_from_gamma = lambda gamma_loc, m_period_loc, m_loc : (gamma_loc + 2 * np.pi * m_period_loc)/ m_loc
alpha_from_gamma_mod = lambda gamma_loc, m_loc : gamma_loc / m_loc

def from_phase_to_real_lab_frame(phi_span_loc, y_loc, i_out_loc):

    R_one_turn =  np.array([])
    Z_one_turn =  np.array([])
    theta_one_turn =  np.array([])

    psi = y_loc [0, :]
    theta = y_loc [1, :]

    R, Z = fi.psi_theta_to_R_Z(psi,theta, flag_vartheta = ip.flag_vartheta)
    
    i_temp = 0
    for i_phi in i_out_loc:
        R_one_turn = np.append(R_one_turn, R[i_phi])
        Z_one_turn = np.append(Z_one_turn, Z[i_phi])
        theta_one_turn = np.append(theta_one_turn, theta[i_phi])
        i_temp += 1  
    
    # In cartisian coordinates
    X, Y = fi.polar_to_cartesian(R, phi_span_loc)
    return R_one_turn, Z_one_turn, theta_one_turn

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

def make_poincare (phi_span_td, alpha_0, psi_t_0, n_lines_tot):
    R_line_tmp = np.array([])
    Z_line_tmp = np.array([])
    theta_line_tmp = np.array([])
    theta_line_grad_tmp = np.array([])
    
    R_lines = [[0 for _ in range(n_lines_tot)] for _ in range(ip.n_poincare)]
    Z_lines = [[0 for _ in range(n_lines_tot)] for _ in range(ip.n_poincare)]
    theta_lines = [[0 for _ in range(n_lines_tot)] for _ in range(ip.n_poincare)] 
    phi = [0 for _ in range(ip.n_poincare)] 
    
    for i_lines in range(n_lines_tot):
        if i_lines % 1000 == 0:
            print(f"Trace {i_lines} line from {n_lines_tot}")
        
        y1 = np.array([[psi_t_0[i_lines]],[alpha_0[i_lines]]])
        
        if ip.flag_Stoermer_Verlet_td:
            y_td = num_solver(phi_span_td, y1, flag_sv = ip.flag_Stoermer_Verlet_td)       
                
        if ip.flag_scipy_num_td:
            y_td = num_solver(phi_span_td, y1, flag_sn = ip.flag_scipy_num_td, method = ip.method)
            
        if ip.flag_Yoshida_td:
            y_td = num_solver(phi_span_td, y1, flag_yo = ip.flag_Yoshida_td)

            
        for i_poincare in range(ip.n_poincare):
            #i_out = list(range(i_poincare*int(ip.points/ip.n_poincare),i_poincare*int(ip.points/ip.n_poincare)+len(phi_span)-ip.points,ip.points))
            #i_out = list(range(i_poincare,i_poincare+1,1))
            i_out = list(range(ip.n_space_points * i_poincare,ip.n_space_points * i_poincare+1,1))
            R_line_tmp, Z_line_tmp, theta_line_tmp = \
            from_phase_to_real_lab_frame(phi_span_td, y_td, i_out)
                
            theta_line_grad_tmp = fi.normalize_angle(theta_line_tmp) / np.pi * 180
                
            R_lines[i_poincare][i_lines]=R_line_tmp
            Z_lines[i_poincare][i_lines]=Z_line_tmp  
            theta_lines[i_poincare][i_lines]=theta_line_grad_tmp
            
    for i_poincare in range(ip.n_poincare):
        i_out = list(range(ip.n_space_points * i_poincare,ip.n_space_points * i_poincare+1,1))
        phi[i_poincare] = phi_span_td[i_out]
            
    return R_lines, Z_lines, theta_lines, phi