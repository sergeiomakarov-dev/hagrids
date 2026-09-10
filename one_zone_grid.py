import numpy as np
import input_params as ip
import create_points_one_zone as cpo
import fun_time_independent_mh as fi
import fun_time_dependent_mh as fd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pickle
import fun_plot as fp
import time
import magnetic_feild_calculation as mfc

# Define marker styles for the scatter points
markers = ['o', 's', '^', 'D', 'v', 'p', '*', 'h', 'X', '<']  # Extend if needed
#cmap = cm.get_cmap('viridis', ip.nruns)  # Choose colormap with ip.nruns entries

# Start the timer
start_time = time.time()

if ip.flag_sym == 1:
    phi_span_pos = np.linspace(0, ip.n_turns * np.pi / ip.nfp, ip.n_turns * ip.points + 1)   
    phi_span_neg = np.linspace(0, -ip.n_turns * np.pi / ip.nfp, ip.n_turns * ip.points + 1) 
else:
    phi_span_pos = np.linspace(0, ip.n_turns * 2 * np.pi / ip.nfp, ip.n_turns * ip.points + 1)   
    phi_span_neg = []
 
alpha_0 = np.array([])
psi_t_0 = np.array([])

R_oz = np.array([])
Z_oz = np.array([])

R_lines_temp = []
Z_lines_temp = []

R_lines_theta_r = []
Z_lines_theta_r = []

# Br_lines_theta_r = []
# Bvartheta_lines_theta_r = []
# Bvarphi_lines_theta_r = []
# BR_lines_theta_r = []
BZ_lines_theta_r = []
BY_lines_theta_r = []
BX_lines_theta_r = []
B_lines_theta_r = []

R_lines = [0]
Z_lines = [0]
phi = [0]

#R_oz = np.zeros((ip.m*(2*ip.n_kin-2)) * (ip.n_E_core+ip.n_E+ip.n_E_pfr))
#Z_oz = np.zeros((ip.m*(2*ip.n_kin-2)) * (ip.n_E_core+ip.n_E+ip.n_E_pfr))

if ip.flag_type_init == 4: 
    
#-----------------------------------------------------------------------------    
    # End the timer
    end_time = time.time()
    
    # Calculate the elapsed time
    elapsed_time = end_time - start_time  
    
    # Display the execution time
    print(f"Point initialization. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------             
    alpha_arr_temp = cpo.alpha_arr_core_min
    psi_arr_temp = cpo.psi_arr_core_min
            
    R_core_min, Z_core_min = fi.psi_theta_to_R_Z(psi_arr_temp, alpha_arr_temp, flag_vartheta = ip.flag_vartheta)
    
    alpha_arr_temp = cpo.alpha_arr_core_core_1
    psi_arr_temp = cpo.psi_arr_core_core_1
            
    R_core_core_1, Z_core_core_1 = fi.psi_theta_to_R_Z(psi_arr_temp, alpha_arr_temp, flag_vartheta = ip.flag_vartheta)
    
    alpha_arr_temp = cpo.alpha_arr_pfr_pfr_1
    psi_arr_temp = cpo.psi_arr_pfr_pfr_1
            
    R_pfr_pfr_1, Z_pfr_pfr_1 = fi.psi_theta_to_R_Z(psi_arr_temp, alpha_arr_temp, flag_vartheta = ip.flag_vartheta)
    
    alpha_arr_temp = cpo.alpha_arr_pfr_max
    psi_arr_temp = cpo.psi_arr_pfr_max
            
    R_pfr_max, Z_pfr_max = fi.psi_theta_to_R_Z(psi_arr_temp, alpha_arr_temp, flag_vartheta = ip.flag_vartheta)
    
    for ipol in range(ip.m*(2*ip.n_kin-2)+1):
        R_temp, Z_temp = fi.generate_points(R_core_min[ipol], Z_core_min[ipol],\
                                            R_core_core_1[ipol], Z_core_core_1[ipol], ip.n_E_core,\
                                                ip.alpha_core_rad) 
        
        R_oz = np.append(R_oz, R_temp)
        Z_oz = np.append(Z_oz, Z_temp)
        
        R_temp, Z_temp = fi.generate_points(R_core_core_1[ipol], Z_core_core_1[ipol],\
                                            R_pfr_pfr_1[ipol], Z_pfr_pfr_1[ipol], ip.n_E,\
                                                ip.alpha_island_rad)  
        R_oz = np.append(R_oz, R_temp[1:])
        Z_oz = np.append(Z_oz, Z_temp[1:])
        
        R_temp, Z_temp = fi.generate_points(R_pfr_pfr_1[ipol], Z_pfr_pfr_1[ipol],\
                                            R_pfr_max[ipol], Z_pfr_max[ipol], ip.n_E_pfr,\
                                                ip.alpha_pfr_rad)
        R_oz = np.append(R_oz, R_temp[1:])
        Z_oz = np.append(Z_oz, Z_temp[1:])
        
    for i in range(len(R_oz)):
        R_lines_temp.append(np.array(R_oz[i]))
        Z_lines_temp.append(np.array(Z_oz[i]))
    
    psi_t_0, alpha_0 = fi.R_Z_to_psi_theta(R_oz, Z_oz, flag_vartheta = 1)
    n_lines_tot_oz = len(R_oz)
    
    # Trace in the positive direction 
    
#-----------------------------------------------------------------------------    
    # End the timer
    end_time = time.time()
    
    # Calculate the elapsed time
    elapsed_time = end_time - start_time  
    
    # Display the execution time
    print(f"Tracing in the positive direction. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------   
    
    if ip.flag_analyt or ip.flag_Stoermer_Verlet or ip.flag_Stoermer_Verlet_mod or ip.flag_Yoshida:
        R_lines, Z_lines, theta_lines, phi = fi.make_poincare (phi_span_pos, alpha_0, psi_t_0, n_lines_tot_oz)    
    elif ip.flag_Stoermer_Verlet_td or ip.flag_scipy_num_td or ip.flag_Yoshida_td:
        R_lines, Z_lines, theta_lines, phi = fd.make_poincare (phi_span_pos, alpha_0, psi_t_0, n_lines_tot_oz)        
    
    R_lines_theta_r_pos, Z_lines_theta_r_pos = fi.structure_R_Z_arrays(R_lines, Z_lines, ip.n_E_core+ip.n_E+ip.n_E_pfr-2)
    phi_pos = phi
    
    if ip.flag_sym == 1:
        # Trace in the negative direction  
#-----------------------------------------------------------------------------    
        # End the timer
        end_time = time.time()
        
        # Calculate the elapsed time
        elapsed_time = end_time - start_time  
        
        # Display the execution time
        print(f"Tracing in the negative direction. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------  
        
        if ip.flag_analyt or ip.flag_Stoermer_Verlet or ip.flag_Stoermer_Verlet_mod or ip.flag_Yoshida:
            R_lines, Z_lines, theta_lines, phi = fi.make_poincare (phi_span_neg, alpha_0, psi_t_0, n_lines_tot_oz)
        elif ip.flag_Stoermer_Verlet_td or ip.flag_scipy_num_td or ip.flag_Yoshida_td:
            R_lines, Z_lines, theta_lines, phi = fd.make_poincare (phi_span_neg, alpha_0, psi_t_0, n_lines_tot_oz)
        
        R_lines_theta_r_neg, Z_lines_theta_r_neg = fi.structure_R_Z_arrays(R_lines, Z_lines, ip.n_E_core+ip.n_E+ip.n_E_pfr-2)
        phi_neg = phi
        
        # Combine postive and negative 
        
        R_lines_theta_r, Z_lines_theta_r, phi_theta_r = fi.comb_pos_neg(R_lines_theta_r_pos, Z_lines_theta_r_pos,\
                                                                              phi_pos,\
                                                                              R_lines_theta_r_neg, Z_lines_theta_r_neg,\
                                                                              phi_neg)
            
    else:
        R_lines_theta_r = R_lines_theta_r_pos
        Z_lines_theta_r = Z_lines_theta_r_pos
        phi_theta_r = phi_pos

    # Calculate the magentic feild 
    B_lines_theta_r, BX_lines_theta_r, BY_lines_theta_r, BZ_lines_theta_r, _, _, _, _ = fi.mag_field_in_pol_planes(R_lines_theta_r, Z_lines_theta_r, phi_theta_r)

'''          
    R_lines[0] = R_lines_temp
    Z_lines[0] = Z_lines_temp
'''    
#-----------------------------------------------------------------------------    
# End the timer
end_time = time.time()
        
# Calculate the elapsed time
elapsed_time = end_time - start_time  
        
# Display the execution time
print(f"Output and plotting. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------  
if ip.flag_output == 1:
    with open(f"./output/R{ip.path}.pkl", "wb") as file:
        pickle.dump(R_lines_theta_r, file)
    with open(f"./output/Z{ip.path}.pkl", "wb") as file:
        pickle.dump(Z_lines_theta_r, file)
    with open(f"./output/B{ip.path}.pkl", "wb") as file:
        pickle.dump(B_lines_theta_r, file)
'''      
if ip.flag_output == 1:
    with open("R_lines_oz.pkl", "wb") as file:
        pickle.dump(R_lines, file)
    with open("Z_lines_oz.pkl", "wb") as file:
        pickle.dump(Z_lines, file)
       
    plt.figure()
    plt.scatter(R_core_min, Z_core_min, marker='.', s=2, color='black') 
    plt.scatter(R_core_core_1, Z_core_core_1, marker='.', s=2, color='black')
    plt.scatter(R_pfr_pfr_1, Z_pfr_pfr_1, marker='.', s=2, color='black')
    plt.scatter(R_pfr_max, Z_pfr_max, marker='.', s=2, color='black')
    plt.scatter(R_oz, Z_oz, marker='.', s=2, color='red')
    
    plt.xlim(ip.R_min, ip.R_max)
    plt.ylim(ip.Z_min, ip.Z_max)
    plt.show() 
''' 
if ip.flag_plot_grid:
    for i_poincare in range(len(R_lines_theta_r)):    
        R_lines_temp=R_lines_theta_r[i_poincare]
        Z_lines_temp=Z_lines_theta_r[i_poincare]
        plt.figure()
        
        for ipol in range(R_lines_temp.shape[0]-1):
            plt.plot(R_lines_temp[ipol,:], Z_lines_temp[ipol,:], color='tab:blue', linewidth=1)
            if ipol != R_lines_temp.shape[0]-1:
                for irad in range(R_lines_temp.shape[1]):
                    plt.plot(R_lines_temp[ipol:ipol+2,irad], Z_lines_temp[ipol:ipol+2,irad], color='tab:blue', linewidth=1)
            for irad in [0, ip.n_E_core-1,ip.n_E_core+ip.n_E-2,-1]:
                plt.plot(R_lines_temp[ipol:ipol+2,irad], Z_lines_temp[ipol:ipol+2,irad], color='tab:red')
        plt.axis('equal')      
        plt.xlim(ip.R_min, ip.R_max)
        plt.ylim(ip.Z_min, ip.Z_max)
        plt.title(r'2D base mesh; one-zone', fontsize=20)        
        plt.xlabel(r'R, m', fontsize=20)
        plt.ylabel(r'Z, m', fontsize=20)
        plt.tick_params(axis='both', labelsize=16) 
        plt.show()
        plt.tight_layout()

if ip.flag_plot_field_lines_n_magnetic_field:
    fig=plt.figure()
    ax=fig.add_subplot(111,projection='3d')
    for i_theta in range(R_lines_theta_r[0].shape[0]):
        for i_r in range(R_lines_theta_r[0].shape[1]):
            X_one_line = []
            Y_one_line = []
            Z_one_line = []
    
            
            for i_poincare in range(len(R_lines_theta_r)):
                R_temp = R_lines_theta_r[i_poincare][i_theta,i_r]
                Z_temp = Z_lines_theta_r[i_poincare][i_theta,i_r]
                phi_temp = phi_theta_r[i_poincare][0]
                
                X_temp, Y_temp= fi.polar_to_cartesian(R_temp, phi_temp)
                
                BX_temp = BX_lines_theta_r[i_poincare][i_theta,i_r]
                BY_temp = BY_lines_theta_r[i_poincare][i_theta,i_r]
                BZ_temp = BZ_lines_theta_r[i_poincare][i_theta,i_r]
                
                ax.quiver(X_temp, Y_temp, Z_temp, 
                          BX_temp, BY_temp, BZ_temp, 
                          length=0.5, color='r')
                
                X_one_line.append(X_temp)
                Y_one_line.append(Y_temp)
                Z_one_line.append(Z_temp)  
     
            ax.plot(X_one_line, Y_one_line, Z_one_line)
            #plt.axis('equal')
#-----------------------------------------------------------------------------    
# End the timer
end_time = time.time()
        
# Calculate the elapsed time
elapsed_time = end_time - start_time  
        
# Display the execution time
print(f"End. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------  