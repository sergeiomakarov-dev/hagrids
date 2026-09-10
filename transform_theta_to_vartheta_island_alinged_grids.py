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

if flag_sym == 1:
    phi_span_pos = np.linspace(0, ip.n_turns * np.pi / ip.nfp, len(R_lines_theta_r_core_in))   
    phi_span_neg = np.linspace(0, -ip.n_turns * np.pi / ip.nfp, len(R_lines_theta_r_core_in))
    
    phi_theta_r = [0] * (len(R_lines_theta_r_core_in))
        
    # phi_theta_r[len(phi_span_pos)-1:2 * len(phi_span_pos) - 1] = \
    # [np.array([phi_span_pos[i]]) for i in range(len(phi_span_pos))]

    # phi_theta_r[0:len(phi_span_pos)-1] = \
    # phi_span_neg[:0:-1]
    
    # Transfer updated phi_theta_r slice as separate numpy arrays
    phi_theta_r[len(phi_span_pos)-1:2 * len(phi_span_pos) - 1] = \
        [np.array([val]) for val in phi_span_pos[0:]]

    # Transfer negative phi_span_neg[:0:-1] values as separate numpy arrays
    phi_theta_r[0:len(phi_span_pos)-1] = \
        [np.array([val]) for val in phi_span_neg[:0:-1]]
else:
    phi_span_pos = np.linspace(0, ip.n_turns * 2 * np.pi / ip.nfp, len(R_lines_theta_r_core_in))   
    phi_span_neg = []  
    
    phi_theta_r = [np.array([val]) for val in phi_span_pos[0:]]
        
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