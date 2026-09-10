#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import warnings

# common input paramenters 
nfp = 1

n_poincare = 34
#n_poincare = 7

alpha_0 = 5
psi_t_min = 0.7
psi_t_min = 0.5
psi_t_min = 0.65
#psi_t_min = 0.96
#psi_t_min = 0.94
#psi_t_min = 1.0050510257216823+0.01
#psi_t_min = 1.3

psi_t_max = 1.5

B_0=1
#R_0=1.75
R_0=5.9

n_turns = 1
points = 165

n_space_points = 1

iota_b = 0.5
iota_a = 0.4

iota_b = 0.7
iota_a = 0.3

flag_filter = 1

R_min = 2.15
R_max = 2.45
Z_min = -0.40
Z_max = 0.00
n_filter = 1

R_min = 1.0
R_max = 2.5
Z_min = -1.0
Z_max = 1.0

#R_min = 2.15
#R_max = 2.45
#Z_min = -0.20
#Z_max = 0.20


R_min = 5.15
R_max = 6.65
Z_min = -1.0
Z_max = 1.0

n_filter = 1

method = 'Radau'
#method = 'RK45'
#method = 'BDF'
#method = 'LSODA'

# time idependent hamiltonian method input parameters 
m = 5
n = 5

A = -6e-3
k = 1.0

flag_type_init = 3

# if flag_type_init = 0
n_psi_surf = 10
n_lines_p_surf = 50

#n_psi_surf = 1
#n_lines_p_surf = 1

# if flag_type_init = 1
R_start = 5.2
R_end = 5.43
#R_start = 6.33
#R_end = 6.61

Z_const = 0.0000

npoints = 100

# if flag_type_init = 3
n_E = 37
n_E_core = 8
n_E_pfr = 12

n_kin = 63

min_E_fac = 250

frac_d_1 = 0.07

    #n_lines_tot = 37
    
alpha_island_rad = 3
alpha_core_rad = 0.4
alpha_pfr_rad = 2.4

alpha_pol = 2

psi_width_2 = 0.015
psi_shift_2 = 1.73

Q_width_2 = 2.3040000000000003e-05
Q_shift_2 = -0.00617088

# if flag_type_init = 5

R_center = 6.327
Z_center = -0.351
box_side_len = 0.01

n_points_1D = 10
plot_turns = 25


#########################################

flag_analyt = 1
flag_analyt_approx = 0
flag_Stoermer_Verlet = 0
flag_Stoermer_Verlet_mod = 0 
flag_scipy_num = 0
flag_Yoshida = 0

flag_type_A_fun = 0

flag_f_modificator = 0

# time dependent hamiltonian method input parameters
m1 = 11
n1 = 10

#m1 = 10
#n1 = 10

m2 = 5
n2 = 5

#m2 = 10
#n2 = 11

m3 = 9
n3 = 10

A1 = -16e-5
A2 = -6e-3
A3 = 0

A1 = 0
#A2 = 0

#A1 = -0.0020829166875
#A2 = -0.0020829166875

n_perturbation = 1

A_list = [-24e-4]
m_list = [11]
n_list = [10]

n_psi_surf_td = 10
n_lines_p_surf_td = 50

#n_psi_surf_td = 1
#n_lines_p_surf_td = 1

flag_Stoermer_Verlet_td = 0
flag_scipy_num_td = 0
flag_Yoshida_td = 0

n_test_points = 50

path = "_lines_theta_r_oz"

# configs input parameters

flag_config = 1


psi_0_min_configs = 1.23
psi_0_max_configs = 1.27

# =============================================================================
# theta_0_min_configs = 0
# theta_0_max_configs = theta_0_min_configs
# =============================================================================
psi_0_min_configs = 1.0050510257216823+0.01
psi_0_max_configs = 1.0050510257216823+0.02

theta_0_o_point = 0
psi_0_o_point = 1.25

n_configs = 10

n_configs = 1000

psi_0_min_configs = 1.0050510257216823+0.01
psi_0_max_configs = 1.0050510257216823+0.02
n_configs = 1

n_turns_one_fl=20

# Plots input parameters

flag_iota_profile = 0
flag_phase_line = 0
flag_animation_phase_line = 0
flag_poincare = 0
flag_poincare_2 = 0
flag_poincare_3 = 0

flag_plot_configs = 0
flag_plot_configs_ani = 0
flag_plot_one_fl = 0
flag_plot_one_fl_1D_plot = 0
flag_plot_field_lines_n_magnetic_field = 0

flag_plot_grid = 1
flag_plot_cell_size = 1

flag_output = 1

flag_sym = 1

flag_vartheta = 1

flag_poloidal_initialization = 0

s_size = 4

flag_iota_res_surf = 1
n_points_iota_res = 300

#-------------------------------------------------------
# Function to read variables from text file
def read_variables(file_path):
    variables = {}
    encountered_variables = set()
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('#') or not line:
                continue  # Skip comments and empty lines
            # Split line into parts to extract variable and value, ignoring comments
            parts = line.split('#', 1)[0].split('=')
            if len(parts) == 2:
                name, value = parts[0].strip(), parts[1].strip()
                if name in encountered_variables:
                    print(f"Warning: Variable '{name}' appears more than once in the file.")
                encountered_variables.add(name)
                variables[name] = value
    return variables

# Function to update variables from the variables dictionary
def update_variables(variables):
    existing_vars = globals().copy()  # Make a copy of the current global variables
    for name, value in variables.items():
        if name in existing_vars:
            globals()[name] = eval(value)  # Use eval to convert string values to integers
        else:
            print(f"Warning: Variable '{name}' not found in the script.")

# Paths to the text file
variables_file_path = 'input_params.dat'

# Read variables from the text file
variables = read_variables(variables_file_path)

# Update the variables in the script with the values from the text file
update_variables(variables)

#-------------------------------------------------------
# Calculations based on input parameters 

# time idependent hamiltonian method input parameters 
iota_res = n/m

psi_0 = (iota_res - iota_b) / iota_a

if flag_type_init == 0:
    n_lines_tot = n_psi_surf * n_lines_p_surf
    
if flag_type_init == 1:
    n_lines_tot = 100

if flag_type_init == 3:
    n_lines_tot = n_E * (4*n_kin-3)
    n_lines_tot_core = n_E_core * (4*n_kin-3)
    n_lines_tot_pfr = n_E_pfr * (4*n_kin-3)
    
# time dependent hamiltonian method input parameters    

iota_res1 = n1 / m1
iota_res2 = n2 / m2
iota_res3 = n3 / m3

A_arr = np.array(A_list)
m_arr = np.array(m_list)
n_arr = np.array(n_list)

iota_res_arr = np.zeros(len(m_arr))
iota_res_arr[:] = n_arr[:] / m_arr[:]
    
n_lines_tot_td = n_psi_surf_td * n_lines_p_surf_td

# other input parameters calculations

theta_0_min_configs = 2 * np.pi / 10 
theta_0_max_configs = theta_0_min_configs


if flag_type_A_fun == 1 and (flag_analyt == 1 or flag_analyt_approx == 1 or flag_Stoermer_Verlet_mod == 1 
                             or flag_scipy_num == 1):
    # Raise a warning message
    warnings.warn(
        "The flag_type_A_fun == 1 option is not supported together with "
        "flag_analyt == 1"
        "flag_analyt_approx == 1, flag_Stoermer_Verlet_mod == 1, "
        "or flag_scipy_num == 1", 
        UserWarning
    )