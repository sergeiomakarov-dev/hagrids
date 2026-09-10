import pickle

import numpy as np
import matplotlib.pyplot as plt
import input_params as ip

# Load the object from the file
#with open("R_lines.pkl", "rb") as file:
path = "_test_const_ampl"
path_2 = "_test_lin_ampl"
path_2 = "_test_lin_ampl_corr"
#path = "_test_const_ampl_3_poincare"
#path_2 = "_test_lin_ampl_3_poincare"
path = "_fun1_7z_test"
path_2 = "_fun1_7z_test_vartheta"

path = "_fun0_7z_test"
path_2 = "_fun1_7z_test"

path ="_fun1_yotd_7z_A2_0_min_025"
path_2 = "_fun1_yotd_7z_A2_m24em4_min_025"

        
with open(f"./output/R_lines_theta_r_core{path}.pkl", "rb") as file:
    R_lines_theta_r_core = pickle.load(file)
# Load the object from the file
#with open("Z_lines.pkl", "rb") as file:
with open(f"./output/Z_lines_theta_r_core{path}.pkl", "rb") as file:
    Z_lines_theta_r_core = pickle.load(file)

with open(f"./output/R_lines_theta_r_pfr{path}.pkl", "rb") as file:
    R_lines_theta_r_pfr = pickle.load(file)
# Load the object from the file
#with open("Z_lines.pkl", "rb") as file:
with open(f"./output/Z_lines_theta_r_pfr{path}.pkl", "rb") as file:
    Z_lines_theta_r_pfr = pickle.load(file)
    
R_lines_theta_r_island = [0 for _ in range(ip.m)]
Z_lines_theta_r_island = [0 for _ in range(ip.m)]
    
for m_period_tmp in range(ip.m):
    with open(f"./output/R_lines_theta_r_island_{m_period_tmp}{path}.pkl", "rb") as file:  
        R_lines_theta_r_island[m_period_tmp] = pickle.load(file)
    with open(f"./output/Z_lines_theta_r_island_{m_period_tmp}{path}.pkl", "rb") as file:  
        Z_lines_theta_r_island[m_period_tmp] = pickle.load(file)
        
with open(f"./output/R_lines_theta_r_core{path_2}.pkl", "rb") as file:
    R_lines_theta_r_core_2 = pickle.load(file)
# Load the object from the file
#with open("Z_lines.pkl", "rb") as file:
with open(f"./output/Z_lines_theta_r_core{path_2}.pkl", "rb") as file:
    Z_lines_theta_r_core_2 = pickle.load(file)

with open(f"./output/R_lines_theta_r_pfr{path_2}.pkl", "rb") as file:
    R_lines_theta_r_pfr_2 = pickle.load(file)
# Load the object from the file
#with open("Z_lines.pkl", "rb") as file:
with open(f"./output/Z_lines_theta_r_pfr{path_2}.pkl", "rb") as file:
    Z_lines_theta_r_pfr_2 = pickle.load(file)
    
R_lines_theta_r_island_2 = [0 for _ in range(ip.m)]
Z_lines_theta_r_island_2 = [0 for _ in range(ip.m)]
    
for m_period_tmp in range(ip.m):
    with open(f"./output/R_lines_theta_r_island_{m_period_tmp}{path_2}.pkl", "rb") as file:  
        R_lines_theta_r_island_2[m_period_tmp] = pickle.load(file)
    with open(f"./output/Z_lines_theta_r_island_{m_period_tmp}{path_2}.pkl", "rb") as file:  
        Z_lines_theta_r_island_2[m_period_tmp] = pickle.load(file)

for i_poincare in range(len(R_lines_theta_r_core)):
    
    plt.figure()
    
    R_lines_temp=R_lines_theta_r_core[i_poincare]
    Z_lines_temp=Z_lines_theta_r_core[i_poincare]

    for ipol in range(R_lines_temp.shape[0]-1):
        plt.plot(R_lines_temp[ipol,:], Z_lines_temp[ipol,:], color='red')
        if ipol != R_lines_temp.shape[0]-1:
            for irad in range(R_lines_temp.shape[1]):
                plt.plot(R_lines_temp[ipol:ipol+2,irad], Z_lines_temp[ipol:ipol+2,irad], color='red')
    
    R_lines_temp=R_lines_theta_r_pfr[i_poincare]
    Z_lines_temp=Z_lines_theta_r_pfr[i_poincare]
                
    for ipol in range(R_lines_temp.shape[0]-1):
        plt.plot(R_lines_temp[ipol,:], Z_lines_temp[ipol,:], color='red')
        if ipol != R_lines_temp.shape[0]-1:
            for irad in range(R_lines_temp.shape[1]):
                plt.plot(R_lines_temp[ipol:ipol+2,irad], Z_lines_temp[ipol:ipol+2,irad], color='red')   
                
    for m_period_tmp in range(ip.m):
        R_lines_temp=R_lines_theta_r_island[m_period_tmp][i_poincare]
        Z_lines_temp=Z_lines_theta_r_island[m_period_tmp][i_poincare]
        
        for ipol in range(R_lines_temp.shape[0]-1):
            plt.plot(R_lines_temp[ipol,:], Z_lines_temp[ipol,:], color='red')
            if ipol != R_lines_temp.shape[0]-1:
                for irad in range(R_lines_temp.shape[1]):
                    plt.plot(R_lines_temp[ipol:ipol+2,irad], Z_lines_temp[ipol:ipol+2,irad], color='red')   

    R_lines_temp=R_lines_theta_r_core_2[i_poincare]
    Z_lines_temp=Z_lines_theta_r_core_2[i_poincare]

    for ipol in range(R_lines_temp.shape[0]-1):
        plt.plot(R_lines_temp[ipol,:], Z_lines_temp[ipol,:], color='black')
        if ipol != R_lines_temp.shape[0]-1:
            for irad in range(R_lines_temp.shape[1]):
                plt.plot(R_lines_temp[ipol:ipol+2,irad], Z_lines_temp[ipol:ipol+2,irad], color='black')
    
    R_lines_temp=R_lines_theta_r_pfr_2[i_poincare]
    Z_lines_temp=Z_lines_theta_r_pfr_2[i_poincare]
                
    for ipol in range(R_lines_temp.shape[0]-1):
        plt.plot(R_lines_temp[ipol,:], Z_lines_temp[ipol,:], color='black')
        if ipol != R_lines_temp.shape[0]-1:
            for irad in range(R_lines_temp.shape[1]):
                plt.plot(R_lines_temp[ipol:ipol+2,irad], Z_lines_temp[ipol:ipol+2,irad], color='black')   
                
    for m_period_tmp in range(ip.m):
        R_lines_temp=R_lines_theta_r_island_2[m_period_tmp][i_poincare]
        Z_lines_temp=Z_lines_theta_r_island_2[m_period_tmp][i_poincare]
        
        for ipol in range(R_lines_temp.shape[0]-1):
            plt.plot(R_lines_temp[ipol,:], Z_lines_temp[ipol,:], color='black')
            if ipol != R_lines_temp.shape[0]-1:
                for irad in range(R_lines_temp.shape[1]):
                    plt.plot(R_lines_temp[ipol:ipol+2,irad], Z_lines_temp[ipol:ipol+2,irad], color='black')   
        
    plt.axis('equal')      
    plt.xlim(ip.R_min, ip.R_max)
    plt.ylim(ip.Z_min, ip.Z_max)
    plt.show() 
    plt.title(r'Grid', fontsize=14)        
    plt.xlabel(r'R, m', fontsize=14)
    plt.ylabel(r'Z, m', fontsize=14)
    plt.show()