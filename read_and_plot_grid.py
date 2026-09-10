import pickle

import numpy as np
import matplotlib.pyplot as plt
import input_params as ip

# Load the object from the file
#with open("R_lines.pkl", "rb") as file:
with open("./output/R_lines_theta_r_oz.pkl", "rb") as file:
    R_lines_theta_r = pickle.load(file)
# Load the object from the file
#with open("Z_lines.pkl", "rb") as file:
with open("./output/Z_lines_theta_r_oz.pkl", "rb") as file:
    Z_lines_theta_r = pickle.load(file)
    
#for i_poincare in range(ip.n_poincare):
#for i_poincare in range(3):    
for i_poincare in range(len(R_lines_theta_r)):
    R_lines_temp=R_lines_theta_r[i_poincare]
    Z_lines_temp=Z_lines_theta_r[i_poincare]
    plt.figure()
    
    for ipol in range(R_lines_temp.shape[0]-1):
        plt.plot(R_lines_temp[ipol,:], Z_lines_temp[ipol,:], color='black')
        if ipol != R_lines_temp.shape[0]-1:
            for irad in range(R_lines_temp.shape[1]):
                plt.plot(R_lines_temp[ipol:ipol+2,irad], Z_lines_temp[ipol:ipol+2,irad], color='black')
        
    plt.axis('equal')      
    plt.xlim(ip.R_min, ip.R_max)
    plt.ylim(ip.Z_min, ip.Z_max)
    plt.show() 

'''    
    n_pol = int(len(R_lines_temp) / n_rad)
    
    for ipol in range(n_pol-1):
        plt.plot(R_lines_temp[ipol*n_rad:(ipol+1)*n_rad], Z_lines_temp[ipol*n_rad:(ipol+1)*n_rad], color='black') 
        if ipol != n_pol:
            for irad in range(n_rad):
                indices=[ipol*n_rad+irad,(ipol+1)*n_rad+irad]
                plt.plot([R_lines_temp[i] for i in indices], [Z_lines_temp[i] for i in indices], color='black')
'''
