import pickle

import numpy as np
import matplotlib.pyplot as plt
import input_params as ip

# Load the object from the file
with open("./output/R_lines_theta_r_oz.pkl", "rb") as file:
    R_lines = pickle.load(file)
# Load the object from the file
with open("./output/Z_lines_theta_r_oz.pkl", "rb") as file:
    Z_lines = pickle.load(file)
    
for i_poincare in range(len(R_lines)):
#for i_poincare in range(1):    
    R_lines_temp=R_lines[i_poincare]
    Z_lines_temp=Z_lines[i_poincare]
    plt.figure()
    #plt.scatter(R_lines_temp[0:3000], Z_lines_temp[0:3000], marker='.', s=2, color='black')
    plt.scatter(R_lines_temp, Z_lines_temp, marker='.', s=2, color='black') 
    
    plt.xlim(ip.R_min, ip.R_max)
    plt.ylim(ip.Z_min, ip.Z_max)
    plt.show() 
