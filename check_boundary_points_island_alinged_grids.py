import pickle

import numpy as np
import input_params as ip

path = "_fun1_yotd_7z"
path = "_7z_rp_resx2"
path ="_fun0_7z_test"

# Load the objects from the files
with open(f"./output/R_lines_theta_r_core{path}.pkl", "rb") as file:
    R_lines_theta_r_core_in = pickle.load(file)
with open(f"./output/Z_lines_theta_r_core{path}.pkl", "rb") as file:
    Z_lines_theta_r_core_in = pickle.load(file)

with open(f"./output/R_lines_theta_r_pfr{path}.pkl", "rb") as file:
    R_lines_theta_r_pfr_in = pickle.load(file)
with open(f"./output/Z_lines_theta_r_pfr{path}.pkl", "rb") as file:
    Z_lines_theta_r_pfr_in = pickle.load(file)
    
R_lines_theta_r_island_in = [0 for _ in range(ip.m)]
Z_lines_theta_r_island_in = [0 for _ in range(ip.m)]
    
for m_period_tmp in range(ip.m):
    with open(f"./output/R_lines_theta_r_island_{m_period_tmp}{path}.pkl", "rb") as file:  
        R_lines_theta_r_island_in[m_period_tmp] = pickle.load(file)
    with open(f"./output/Z_lines_theta_r_island_{m_period_tmp}{path}.pkl", "rb") as file:  
        Z_lines_theta_r_island_in[m_period_tmp] = pickle.load(file)

nan_num = 0

print("Check NaN numbers")
for i_poincare in range(len(R_lines_theta_r_core_in)):
    temp = np.isnan(R_lines_theta_r_core_in[i_poincare])
    nan_indices = np.where(temp)  # Extract indices
    if temp.any():
        print(f"Core: Pol. pl.: {i_poincare}: R Indices of NaN values:", nan_indices)
        nan_num += 1
    temp = np.isnan(Z_lines_theta_r_core_in[i_poincare])
    nan_indices = np.where(temp)  # Extract indices
    if temp.any():
        print(f"Core: Pol. pl.: {i_poincare}: Z Indices of NaN values:", nan_indices)
        nan_num += 1

    temp = np.isnan(R_lines_theta_r_pfr_in[i_poincare])
    nan_indices = np.where(temp)  # Extract indices
    if temp.any():
        print(f"PFR: Pol. pl.: {i_poincare}: R Indices of NaN values:", nan_indices)
        nan_num += 1
    temp = np.isnan(Z_lines_theta_r_pfr_in[i_poincare])
    nan_indices = np.where(temp)  # Extract indices
    if temp.any():
        print(f"PFR: Pol. pl.: {i_poincare}: Z Indices of NaN values:", nan_indices)
        nan_num += 1
        
    for m_period_tmp in range(ip.m): 
        temp = np.isnan(R_lines_theta_r_island_in[m_period_tmp][i_poincare])
        nan_indices = np.where(temp)  # Extract indices
        if temp.any():
            print(f"Island No {m_period_tmp}: Pol. pl.: {i_poincare}: R Indices of NaN values:", nan_indices)
            nan_num += 1
        temp = np.isnan(Z_lines_theta_r_island_in[m_period_tmp][i_poincare])
        nan_indices = np.where(temp)  # Extract indices
        if temp.any():
            print(f"Island No {m_period_tmp}: Pol. pl.: {i_poincare}: Z Indices of NaN values:", nan_indices)
            nan_num += 1
    
if nan_num == 0:
    print("No NaN numbers found")
else:
    print(f"Number of NaN's: {nan_num}")     
 
no_overlap_num = 0

R_b_core_island = [0 for _ in range(len(R_lines_theta_r_core_in))]
Z_b_core_island = [0 for _ in range(len(R_lines_theta_r_core_in))]
print("Check overlap core and island points")
for i_poincare in range(len(R_lines_theta_r_core_in)):
    
    R_b_core_island[i_poincare]=[0 for _ in range(ip.m)]
    Z_b_core_island[i_poincare]=[0 for _ in range(ip.m)]

    for m_period_tmp in range(ip.m):    
        i_rad_begin = 0
        i_rad_end = R_lines_theta_r_core_in[0].shape[1] - 1
        n_pol_island = R_lines_theta_r_island_in[m_period_tmp][i_poincare].shape[0]
        n_pol_check=int((n_pol_island-1)/2)
        R_b_core_island[i_poincare][m_period_tmp] = \
        R_lines_theta_r_island_in[m_period_tmp][i_poincare][0:n_pol_check+1,i_rad_begin] - \
        R_lines_theta_r_core_in[i_poincare][m_period_tmp*n_pol_check:(m_period_tmp+1)*n_pol_check+1,i_rad_end]
        Z_b_core_island[i_poincare][m_period_tmp] = \
        Z_lines_theta_r_island_in[m_period_tmp][i_poincare][0:n_pol_check+1,i_rad_begin] - \
        Z_lines_theta_r_core_in[i_poincare][m_period_tmp*n_pol_check:(m_period_tmp+1)*n_pol_check+1,i_rad_end]
        
        for i_pol in range(R_b_core_island[i_poincare][m_period_tmp].shape[0]): 
            if R_b_core_island[i_poincare][m_period_tmp][i_pol] != 0.0 or Z_b_core_island[i_poincare][m_period_tmp][i_pol] != 0.0:
            
                print(f"Core: Pol. pl.: {i_poincare}, m: {m_period_tmp}, i_pol: {i_pol}. Diff: R = {R_b_core_island[i_poincare][m_period_tmp][i_pol]:.1e}, Z = {Z_b_core_island[i_poincare][m_period_tmp][i_pol]:.1e}")
                no_overlap_num += 1

R_b_pfr_island = [0 for _ in range(len(R_lines_theta_r_core_in))]
Z_b_pfr_island = [0 for _ in range(len(R_lines_theta_r_core_in))]
print("Check overlap PFR and island points")
for i_poincare in range(len(R_lines_theta_r_core_in)):
    
    R_b_pfr_island[i_poincare]=[0 for _ in range(ip.m)]
    Z_b_pfr_island[i_poincare]=[0 for _ in range(ip.m)]

    for m_period_tmp in range(ip.m):    
        i_rad_begin = 0
        i_rad_end = R_lines_theta_r_core_in[0].shape[1] - 1
        n_pol_island = R_lines_theta_r_island_in[m_period_tmp][i_poincare].shape[0]
        n_pol_check=int((n_pol_island-1)/2)
        R_b_pfr_island[i_poincare][m_period_tmp] = \
        R_lines_theta_r_island_in[m_period_tmp][i_poincare][2*n_pol_check+1:n_pol_check-1:-1,i_rad_begin] - \
        R_lines_theta_r_pfr_in[i_poincare][m_period_tmp*n_pol_check:(m_period_tmp+1)*n_pol_check+1,i_rad_begin]
        Z_b_pfr_island[i_poincare][m_period_tmp] = \
        Z_lines_theta_r_island_in[m_period_tmp][i_poincare][2*n_pol_check+1:n_pol_check-1:-1,i_rad_begin] - \
        Z_lines_theta_r_pfr_in[i_poincare][m_period_tmp*n_pol_check:(m_period_tmp+1)*n_pol_check+1,i_rad_begin]
        
        for i_pol in range(R_b_pfr_island[i_poincare][m_period_tmp].shape[0]): 
            if R_b_pfr_island[i_poincare][m_period_tmp][i_pol] != 0.0 or Z_b_pfr_island[i_poincare][m_period_tmp][i_pol] != 0.0:
            
                print(f"PFR: Pol. pl.: {i_poincare}, m: {m_period_tmp}, i_pol: {i_pol}. Diff: R = {R_b_pfr_island[i_poincare][m_period_tmp][i_pol]:.1e}, Z = {Z_b_pfr_island[i_poincare][m_period_tmp][i_pol]:.1e}")
                no_overlap_num += 1                
if no_overlap_num == 0:
    print("All boudary points overlap.")
else:
    print(f"Number of points without overlap: {no_overlap_num}")        
   