import numpy as np
from scipy.special import ellipj, ellipk, ellipkinc
from scipy.integrate import solve_ivp
import input_params as ip
import create_points_in_the_island as cp
import fun_time_independent_mh as fi

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
        psi_t_0_temp, alpha_0_temp = fi.R_Z_to_psi_theta(R_span[i_lines], Z_span[i_lines], flag_vartheta = 1)
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
        psi_t_0_temp, alpha_0_temp = fi.R_Z_to_psi_theta(R_span[i_lines], Z_span[i_lines], flag_vartheta = 1)
        alpha_0 = np.append(alpha_0, alpha_0_temp)
        psi_t_0 = np.append(psi_t_0, psi_t_0_temp)
        
elif ip.flag_type_init == 3: 
    m_period_tmp = 4
    
    alpha_arr_temp = cp.alpha_arr_isl[m_period_tmp]
    psi_arr_temp = cp.psi_arr_isl[m_period_tmp]
    
    alpha_0 = alpha_arr_temp.flatten(order='F')   
    psi_t_0 = psi_arr_temp.flatten(order='F')

R_lines, Z_lines, theta_lines = fi.make_poincare (phi_span, alpha_0, psi_t_0)