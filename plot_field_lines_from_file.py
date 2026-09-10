import pickle

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import input_params as ip

path = "_lines_theta_r_oz_test_an_with_island"
path_1 = "_lines_theta_r_oz_test_ti_with_island"
path_2 = "_lines_theta_r_oz_test_ti_mod_with_island"
path_3 = "_lines_theta_r_oz_test_td_with_island"
path_3 ="_lines_theta_r_oz_test_yoshida"
path_4 = "_lines_theta_r_oz_test_td_sn_with_island"

path_1 = "_lines_theta_r_oz_test_td_with_island"
path_2 = "_lines_theta_r_oz_test_yoshida_td"

'''
path = "_lines_theta_r_oz_m1_ti_with_island"
path_1 = "_lines_theta_r_oz_m1_ti_mod_with_island"
path_2 = "_lines_theta_r_oz_m1_td_with_island"
'''
'''
path = "_lines_theta_r_oz_an"
path_1 = "_lines_theta_r_oz_fixed_hresolut"
'''

# Load the object from the file
with open(f"./output/R{path}.pkl", "rb") as file:
    R_lines = pickle.load(file)
# Load the object from the file
with open(f"./output/Z{path}.pkl", "rb") as file:
    Z_lines = pickle.load(file)
    
# Load the object from the file
with open(f"./output/R{path_1}.pkl", "rb") as file:
    R_lines_1 = pickle.load(file)
# Load the object from the file
with open(f"./output/Z{path_1}.pkl", "rb") as file:
    Z_lines_1 = pickle.load(file)
    
# Load the object from the file
with open(f"./output/R{path_2}.pkl", "rb") as file:
    R_lines_2 = pickle.load(file)
# Load the object from the file
with open(f"./output/Z{path_2}.pkl", "rb") as file:
    Z_lines_2 = pickle.load(file)

# Load the object from the file
with open(f"./output/R{path_3}.pkl", "rb") as file:
    R_lines_3 = pickle.load(file)
# Load the object from the file
with open(f"./output/Z{path_3}.pkl", "rb") as file:
    Z_lines_3 = pickle.load(file)
    
# Load the object from the file
with open(f"./output/R{path_4}.pkl", "rb") as file:
    R_lines_4 = pickle.load(file)
# Load the object from the file
with open(f"./output/Z{path_4}.pkl", "rb") as file:
    Z_lines_4 = pickle.load(file)

if ip.flag_sym == 1:
    phi_span_pos = np.linspace(0, ip.n_turns * np.pi / ip.nfp, ip.n_poincare)   
    phi_span_neg = np.linspace(-ip.n_turns * np.pi / ip.nfp, 0, ip.n_poincare) 

    phi_span = np.append(phi_span_neg, phi_span_pos[1:])
else:
    phi_span = np.linspace(0, ip.n_turns * 2 * np.pi / ip.nfp, ip.n_poincare + 1)   

phi_span_list = []
for i in range(len(phi_span))    :
    phi_span_list.append(np.array(phi_span[i]))

test_pol = 3

fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
for i_line in range(ip.n_E_core+ip.n_E+ip.n_E_pfr-2):
    R_one_line = []
    Z_one_line = []
    for i_poincare in range(2 * ip.n_poincare - 1):
        R_one_line.append(R_lines[i_poincare][test_pol,i_line])
        Z_one_line.append(Z_lines[i_poincare][test_pol,i_line])   
    ax.plot(phi_span_list, R_one_line, Z_one_line, color='black')
    
    R_one_line_1 = []
    Z_one_line_1 = []    
    for i_poincare in range(2 * ip.n_poincare - 1):
        R_one_line_1.append(R_lines_1[i_poincare][test_pol,i_line])
        Z_one_line_1.append(Z_lines_1[i_poincare][test_pol,i_line]) 
    ax.plot(phi_span_list, R_one_line_1, Z_one_line_1, color='red')    

    R_one_line_2 = []
    Z_one_line_2 = []    
    for i_poincare in range(2 * ip.n_poincare - 1):
        R_one_line_2.append(R_lines_2[i_poincare][test_pol,i_line])
        Z_one_line_2.append(Z_lines_2[i_poincare][test_pol,i_line]) 
    ax.plot(phi_span_list, R_one_line_2, Z_one_line_2, color='blue')
    '''
    R_one_line_3 = []
    Z_one_line_3 = []    
    for i_poincare in range(2 * ip.n_poincare - 1):
        R_one_line_3.append(R_lines_3[i_poincare][test_pol,i_line])
        Z_one_line_3.append(Z_lines_3[i_poincare][test_pol,i_line]) 
    ax.plot(phi_span_list, R_one_line_3, Z_one_line_3, color='green')
    
    R_one_line_4 = []
    Z_one_line_4 = []    
    for i_poincare in range(2 * ip.n_poincare - 1):
        R_one_line_4.append(R_lines_4[i_poincare][test_pol,i_line])
        Z_one_line_4.append(Z_lines_4[i_poincare][test_pol,i_line]) 
    ax.plot(phi_span_list, R_one_line_4, Z_one_line_4, color='orange')
    '''      
    plt.show() 

