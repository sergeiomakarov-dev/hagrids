import pickle

import numpy as np
import input_params as ip
import fun_time_independent_mh as fi
import time

path_in = "_fun0_7z_test"
path_out = "_fun0_7z_test_vartheta"

path_in = "_fun1_yotd_7z"
path_out = "_fun1_yotd_7z_vartheta_3"

path_in ="_for_Nassim_no_vartheta"
path_out ="_for_Nassim_with_vartheta"

path_in ="_fun1_yotd_7z_A2_m24em4"
path_out ="_fun1_yotd_7z_vartheta_A2_m24em4"

path_in ="_fun1_yotd_7z_A2_0"
path_out ="_fun1_yotd_7z_vartheta_A2_0"


flag_sym = 0

# Start the timer
start_time = time.time()

#-----------------------------------------------------------------------------    
# End the timer
end_time = time.time()
    
# Calculate the elapsed time
elapsed_time = end_time - start_time  
    
# Display the execution time
print(f"Load the grid from the files. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------  

# Load the objects from the files
with open(f"./output/R_lines_theta_r_core{path_in}.pkl", "rb") as file:
    R_lines_theta_r_core_in = pickle.load(file)
with open(f"./output/Z_lines_theta_r_core{path_in}.pkl", "rb") as file:
    Z_lines_theta_r_core_in = pickle.load(file)
R_lines_theta_r_core_out = [np.zeros(R_lines_theta_r_core_in[0].shape) for _ in range(len(R_lines_theta_r_core_in))]
Z_lines_theta_r_core_out = [np.zeros(Z_lines_theta_r_core_in[0].shape) for _ in range(len(Z_lines_theta_r_core_in))]

with open(f"./output/R_lines_theta_r_pfr{path_in}.pkl", "rb") as file:
    R_lines_theta_r_pfr_in = pickle.load(file)
with open(f"./output/Z_lines_theta_r_pfr{path_in}.pkl", "rb") as file:
    Z_lines_theta_r_pfr_in = pickle.load(file)
R_lines_theta_r_pfr_out = [np.zeros(R_lines_theta_r_pfr_in[0].shape) for _ in range(len(R_lines_theta_r_pfr_in))]
Z_lines_theta_r_pfr_out = [np.zeros(Z_lines_theta_r_pfr_in[0].shape) for _ in range(len(Z_lines_theta_r_pfr_in))]
    
R_lines_theta_r_island_in = [0 for _ in range(ip.m)]
Z_lines_theta_r_island_in = [0 for _ in range(ip.m)]
    
for m_period_tmp in range(ip.m):
    with open(f"./output/R_lines_theta_r_island_{m_period_tmp}{path_in}.pkl", "rb") as file:  
        R_lines_theta_r_island_in[m_period_tmp] = pickle.load(file)
    with open(f"./output/Z_lines_theta_r_island_{m_period_tmp}{path_in}.pkl", "rb") as file:  
        Z_lines_theta_r_island_in[m_period_tmp] = pickle.load(file)
R_lines_theta_r_island_out = [[np.zeros(R_lines_theta_r_island_in[0][0].shape) for _ in range(len(R_lines_theta_r_island_in[0]))] for _ in range(ip.m)]
Z_lines_theta_r_island_out = [[np.zeros(Z_lines_theta_r_island_in[0][0].shape) for _ in range(len(Z_lines_theta_r_island_in[0]))] for _ in range(ip.m)]


B_lines_theta_r_island_out = [[np.zeros(R_lines_theta_r_island_in[0][0].shape) for _ in range(len(R_lines_theta_r_island_in[0]))] for _ in range(ip.m)]

B_lines_theta_r_core_out = [np.zeros(R_lines_theta_r_core_in[0].shape) for _ in range(len(R_lines_theta_r_core_in))]

B_lines_theta_r_pfr_out = [np.zeros(R_lines_theta_r_pfr_in[0].shape) for _ in range(len(R_lines_theta_r_pfr_in))]

# Form a phi array

#phi = [0 for _ in range(len(R_lines_theta_r_core_in))] 

# Toroidal angles of the stored planes, as in island_alinged_grid.py and make_poincare: the field
# lines are traced over n_turns * points steps per direction and every n_space_points-th plane
# is stored, n_poincare planes in total per direction
if flag_sym == 1:
    phi_span_pos = np.linspace(0, ip.n_turns * np.pi / ip.nfp, ip.n_turns * ip.points + 1)
    phi_span_neg = np.linspace(0, -ip.n_turns * np.pi / ip.nfp, ip.n_turns * ip.points + 1)
else:
    phi_span_pos = np.linspace(0, ip.n_turns * 2 * np.pi / ip.nfp, ip.n_turns * ip.points + 1)
    phi_span_neg = []

phi_pos = phi_span_pos[::ip.n_space_points][:ip.n_poincare]
phi_neg = phi_span_neg[::ip.n_space_points][:ip.n_poincare] if flag_sym == 1 else np.array([])

# For flag_sym = 1 the planes traced in the negative direction come first, reversed, and share
# the plane phi = 0 with the planes traced in the positive direction (see comb_pos_neg)
phi_all = np.concatenate((phi_neg[:0:-1], phi_pos))
if len(phi_all) != len(R_lines_theta_r_core_in):
    raise ValueError(f"{len(R_lines_theta_r_core_in)} planes stored in the grid files, but "
                     f"n_poincare = {ip.n_poincare} and flag_sym = {flag_sym} give {len(phi_all)} planes")
phi_theta_r = [np.array([val]) for val in phi_all]

# Perform the trasformation
#-----------------------------------------------------------------------------    
# End the timer
end_time = time.time()
    
# Calculate the elapsed time
elapsed_time = end_time - start_time  
    
# Display the execution time
print(f"Perform the trasformation. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------  

for i_poincare in range(len(R_lines_theta_r_core_in)):
    #-----------------------------------------------------------------------------    
    # End the timer
    end_time = time.time()
        
    # Calculate the elapsed time
    elapsed_time = end_time - start_time  
        
    # Display the execution time
    print(f"Poloidal plain No {i_poincare}. Execution time: {elapsed_time:.4f} seconds")
    #-----------------------------------------------------------------------------  
    for i_pol in range(R_lines_theta_r_core_in[0].shape[0]):
        for i_rad in range(R_lines_theta_r_core_in[0].shape[1]):
            R_temp_out, Z_temp_out = \
            fi.transform_theta_to_vartheta(R_lines_theta_r_core_in[i_poincare][i_pol,i_rad], Z_lines_theta_r_core_in[i_poincare][i_pol,i_rad])
            R_lines_theta_r_core_out[i_poincare][i_pol,i_rad] = R_temp_out[0]
            Z_lines_theta_r_core_out[i_poincare][i_pol,i_rad] = Z_temp_out[0]
    for i_pol in range(R_lines_theta_r_pfr_in[0].shape[0]):
        for i_rad in range(R_lines_theta_r_pfr_in[0].shape[1]):
            R_temp_out, Z_temp_out = \
            fi.transform_theta_to_vartheta(R_lines_theta_r_pfr_in[i_poincare][i_pol,i_rad], Z_lines_theta_r_pfr_in[i_poincare][i_pol,i_rad])
            R_lines_theta_r_pfr_out[i_poincare][i_pol,i_rad] = R_temp_out[0]
            Z_lines_theta_r_pfr_out[i_poincare][i_pol,i_rad] = Z_temp_out[0]
    for m_period_tmp in range(ip.m):
        for i_pol in range(R_lines_theta_r_island_in[m_period_tmp][0].shape[0]):
            for i_rad in range(R_lines_theta_r_island_in[m_period_tmp][0].shape[1]):
                R_temp_out, Z_temp_out = \
                fi.transform_theta_to_vartheta(R_lines_theta_r_island_in[m_period_tmp][i_poincare][i_pol,i_rad], Z_lines_theta_r_island_in[m_period_tmp][i_poincare][i_pol,i_rad])
                R_lines_theta_r_island_out[m_period_tmp][i_poincare][i_pol,i_rad] = R_temp_out[0]
                Z_lines_theta_r_island_out[m_period_tmp][i_poincare][i_pol,i_rad] = Z_temp_out[0]
                
# Calculate the magentic feild
B_lines_theta_r_core_out, \
    _, _, _, _, _, _, _ = fi.mag_field_in_pol_planes(R_lines_theta_r_core_out, \
                                                    Z_lines_theta_r_core_out, 
                                                    phi_theta_r, flag_vartheta=1)
                
B_lines_theta_r_pfr_out, \
    _, _, _, _, _, _, _ = fi.mag_field_in_pol_planes(R_lines_theta_r_pfr_out, \
                                                    Z_lines_theta_r_pfr_out, 
                                                    phi_theta_r, flag_vartheta=1)
for m_period_tmp in range(ip.m):
    B_lines_theta_r_island_out[m_period_tmp], \
        _, _, _, _, _, _, _ = fi.mag_field_in_pol_planes(R_lines_theta_r_island_out[m_period_tmp], \
                                                        Z_lines_theta_r_island_out[m_period_tmp], 
                                                        phi_theta_r, flag_vartheta=1)
            
# Save the objects to the files
#-----------------------------------------------------------------------------    
# End the timer
end_time = time.time()
    
# Calculate the elapsed time
elapsed_time = end_time - start_time  
    
# Display the execution time
print(f"Save the grid to the files. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------  
with open(f"./output/R_lines_theta_r_core{path_out}.pkl", "wb") as file:
    pickle.dump(R_lines_theta_r_core_out, file)
with open(f"./output/Z_lines_theta_r_core{path_out}.pkl", "wb") as file:
    pickle.dump(Z_lines_theta_r_core_out, file)
with open(f"./output/B_lines_theta_r_core{path_out}.pkl", "wb") as file:
    pickle.dump(B_lines_theta_r_core_out, file)

with open(f"./output/R_lines_theta_r_pfr{path_out}.pkl", "wb") as file:
    pickle.dump(R_lines_theta_r_pfr_out, file)
with open(f"./output/Z_lines_theta_r_pfr{path_out}.pkl", "wb") as file:
    pickle.dump(Z_lines_theta_r_pfr_out, file)
with open(f"./output/B_lines_theta_r_pfr{path_out}.pkl", "wb") as file:
    pickle.dump(B_lines_theta_r_pfr_out, file)
    
for m_period_tmp in range(ip.m):
    with open(f"./output/R_lines_theta_r_island_{m_period_tmp}{path_out}.pkl", "wb") as file:
        pickle.dump(R_lines_theta_r_island_out[m_period_tmp], file)
    with open(f"./output/Z_lines_theta_r_island_{m_period_tmp}{path_out}.pkl", "wb") as file:
        pickle.dump(Z_lines_theta_r_island_out[m_period_tmp], file)
    with open(f"./output/B_lines_theta_r_island_{m_period_tmp}{path_out}.pkl", "wb") as file:
        pickle.dump(B_lines_theta_r_island_out[m_period_tmp], file)
        
#-----------------------------------------------------------------------------    
# End the timer
end_time = time.time()
        
# Calculate the elapsed time
elapsed_time = end_time - start_time  
        
# Display the execution time
print(f"End. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------       