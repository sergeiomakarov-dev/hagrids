import numpy as np
import input_params as ip
import fun_time_independent_mh as fi
import fun_time_dependent_mh as fd
import matplotlib.pyplot as plt
import pickle
import time

def generate_box(R_center, Z_center, box_side_len, n_points_1D):
    """
    Generate a box of points in R, Z space.
    
    Parameters:
        R_center (float): Center R coordinate of the box.
        Z_center (float): Center Z coordinate of the box.
        box_side_len (float): Length of the box sides.
        n_points_1D (int): Number of points along one dimension.
    
    Returns:
        np.ndarray: Array of shape (n_points_1D * n_points_1D, 2) containing (R, Z) coordinates.
    """
    half_side = box_side_len / 2.0
    R_vals = np.linspace(R_center - half_side, R_center + half_side, n_points_1D)
    Z_vals = np.linspace(Z_center - half_side, Z_center + half_side, n_points_1D)
    
    R_grid, Z_grid = np.meshgrid(R_vals, Z_vals)
    
    points = np.column_stack([R_grid.ravel(), Z_grid.ravel()])
    return points

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

R_lines_island = [0 for _ in range(ip.m)]
Z_lines_island = [0 for _ in range(ip.m)]

R_lines_theta_r_island = [0 for _ in range(ip.m)]
Z_lines_theta_r_island = [0 for _ in range(ip.m)]

# Br_lines_theta_r_island = [0 for _ in range(ip.m)]
# Bvartheta_lines_theta_r_island = [0 for _ in range(ip.m)]
# Bvarphi_lines_theta_r_island = [0 for _ in range(ip.m)]
# BR_lines_theta_r_island = [0 for _ in range(ip.m)]
BZ_lines_theta_r_island = [0 for _ in range(ip.m)]
BY_lines_theta_r_island = [0 for _ in range(ip.m)]
BX_lines_theta_r_island = [0 for _ in range(ip.m)]
B_lines_theta_r_island = [0 for _ in range(ip.m)]

# Br_lines_theta_r_core = []
# Bvartheta_lines_theta_r_core = []
# Bvarphi_lines_theta_r_core = []
# BR_lines_theta_r_core = []
BZ_lines_theta_r_core = []
BY_lines_theta_r_core = []
BX_lines_theta_r_core = []
B_lines_theta_r_core = []

# Br_lines_theta_r_box = []
# Bvartheta_lines_theta_r_box = []
# Bvarphi_lines_theta_r_box = []
# BR_lines_theta_r_box = []
BZ_lines_theta_r_box = []
BY_lines_theta_r_box = []
BX_lines_theta_r_box = []
B_lines_theta_r_box = []


if ip.flag_type_init == 5:                
    n_lines_tot_box = ip.n_points_1D ** 2
    
    points = generate_box(ip.R_center, ip.Z_center, ip.box_side_len, ip.n_points_1D)
elif ip.flag_type_init == 1:
    n_lines_tot_box = ip.npoints
    
    points = np.zeros((n_lines_tot_box, 2))
    
    points[:,0] = np.linspace(ip.R_start, ip.R_end, n_lines_tot_box)
    points[:,1] = np.linspace(ip.Z_const, ip.Z_const, n_lines_tot_box)

if ip.flag_type_init == 5 or ip.flag_type_init == 1: 

    
    for i_lines in range(n_lines_tot_box):
        psi_0_temp, alpha_0_temp = fi.R_Z_to_psi_theta(points[i_lines][0], points[i_lines][1], flag_vartheta = ip.flag_vartheta)
        alpha_0 = np.append(alpha_0, alpha_0_temp)
        psi_t_0 = np.append(psi_t_0, psi_0_temp)   
            
    # alpha_0 = alpha_arr_temp.flatten(order='F')   
    # psi_t_0 = psi_arr_temp.flatten(order='F')
        
    # Trace in the positive direction 
#-----------------------------------------------------------------------------    
    # End the timer
    end_time = time.time()
    
    # Calculate the elapsed time
    elapsed_time = end_time - start_time  
    
    # Display the execution time
    print(f"box region. Tracing in the positive direction. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------     
    if ip.flag_analyt or ip.flag_Stoermer_Verlet or ip.flag_Stoermer_Verlet_mod or ip.flag_Yoshida:             
        R_lines, Z_lines, theta_lines, phi = fi.make_poincare (phi_span_pos, alpha_0, psi_t_0, n_lines_tot_box) 
    elif ip.flag_Stoermer_Verlet_td or ip.flag_scipy_num_td or ip.flag_Yoshida_td:      
        R_lines, Z_lines, theta_lines, phi = fd.make_poincare (phi_span_pos, alpha_0, psi_t_0, n_lines_tot_box) 
      
    if ip.flag_type_init == 5:
        R_lines_theta_r_pos, Z_lines_theta_r_pos  = fi.structure_R_Z_arrays(R_lines, Z_lines, ip.n_points_1D)
    elif ip.flag_type_init == 1:
        R_lines_theta_r_pos, Z_lines_theta_r_pos  = fi.structure_R_Z_arrays(R_lines, Z_lines, n_lines_tot_box)
    phi_pos = phi
           
    if ip.flag_sym == 1:
        # Trace in the negative direction 
        
#-----------------------------------------------------------------------------    
        # End the timer
        end_time = time.time()
        
        # Calculate the elapsed time
        elapsed_time = end_time - start_time  
        
        # Display the execution time
        print(f"box region. Tracing in the negative direction. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------  
        if ip.flag_analyt or ip.flag_Stoermer_Verlet or ip.flag_Stoermer_Verlet_mod or ip.flag_Yoshida:                         
            R_lines, Z_lines, theta_lines, phi = fi.make_poincare (phi_span_neg, alpha_0, psi_t_0, n_lines_tot_box)
        elif ip.flag_Stoermer_Verlet_td or ip.flag_scipy_num_td or ip.flag_Yoshida_td: 
            R_lines, Z_lines, theta_lines, phi = fd.make_poincare (phi_span_neg, alpha_0, psi_t_0, n_lines_tot_box)
        
        if ip.flag_type_init == 5:
            R_lines_theta_r_neg, Z_lines_theta_r_neg = fi.structure_R_Z_arrays(R_lines, Z_lines, ip.n_points_1D)
        elif ip.flag_type_init == 1:
            R_lines_theta_r_neg, Z_lines_theta_r_neg = fi.structure_R_Z_arrays(R_lines, Z_lines, n_lines_tot_box)
        phi_neg = phi
        
        # Combine postive and negative 
            
        R_lines_theta_r, Z_lines_theta_r, phi_theta_r = fi.comb_pos_neg(R_lines_theta_r_pos, Z_lines_theta_r_pos, \
                                                                        phi_pos,\
                                                                        R_lines_theta_r_neg, Z_lines_theta_r_neg,\
                                                                        phi_neg)
                
    else:
        R_lines_theta_r = R_lines_theta_r_pos
        Z_lines_theta_r = Z_lines_theta_r_pos
        phi_theta_r = phi_pos
    
    R_lines_box=R_lines
    Z_lines_box=Z_lines
    
    R_lines_theta_r_box=R_lines_theta_r
    Z_lines_theta_r_box=Z_lines_theta_r 
    
    #fp.plot_poincare (R_lines, Z_lines, n_lines_tot_core)    
    
            
    if ip.flag_plot_grid or ip.flag_plot_cell_size or ip.flag_poincare:
#-----------------------------------------------------------------------------    
        # End the timer
        end_time = time.time()
            
        # Calculate the elapsed time
        elapsed_time = end_time - start_time  
            
        # Display the execution time
        print(f"Plotting. Execution time: {elapsed_time:.4f} seconds")
#----------------------------------------------------------------------------- 


    if ip.flag_poincare :
        plt.figure()
        if ip.flag_type_init == 5:
            for i_poincare in range(len(R_lines_theta_r_box) - ip.plot_turns, len(R_lines_theta_r_box)): 
                plt.scatter(R_lines_theta_r_box[i_poincare][:,:], Z_lines_theta_r_box[i_poincare][:,:], marker='.', s=ip.s_size, color='red')    
            for i_poincare in range(ip.plot_turns): 
                plt.scatter(R_lines_theta_r_box[i_poincare][:,:], Z_lines_theta_r_box[i_poincare][:,:], marker='.', s=ip.s_size, color='blue')
        elif ip.flag_type_init == 1:
            for i_poincare in range(len(R_lines_theta_r_box)): 
                plt.scatter(R_lines_theta_r_box[i_poincare][:,:], Z_lines_theta_r_box[i_poincare][:,:], marker='.', s=ip.s_size, color='red') 

        # for i_poincare in range(len(R_lines_theta_r_box)): 
        #     plt.scatter(R_lines_theta_r_box[i_poincare][:,:], Z_lines_theta_r_box[i_poincare][:,:], marker='.', s=ip.s_size, color='green')
                
        # if ip.flag_analyt or ip.flag_Stoermer_Verlet or ip.flag_Stoermer_Verlet_mod or ip.flag_Yoshida:
        #     title = r"$A = {:.2e}$".format(ip.A)
        # elif ip.flag_Stoermer_Verlet_td or ip.flag_scipy_num_td or ip.flag_Yoshida_td:
        #     title = r"$A_1 = {:.2e}$; $A_2 = {:.2e}$".format(ip.A1, ip.A2)
        # else:
        #     title = r""        


        # Construct the first row with m, n, A_1
        title_lines = [f"$m_1 = {ip.m}$; $n_1 = {ip.n}$; $A_1 = {ip.A1:.2e}$"]
        
        # Add m_2, n_2, A_2, m_3, n_3, A_3... on the same row
        for i in range(len(ip.m_list)):
            title_lines.append(f"$m_{i+2} = {ip.m_list[i]}$; $n_{i+2} = {ip.n_list[i]}$; $A_{i+2} = {ip.A_list[i]:.2e}$")
        
        # Join all lines with newline characters to make it multi-line
        title_text = "\n".join(title_lines)
                
        # Set the title using plt.title() with the constructed multi-line text
        plt.title(title_text, fontsize=14)

        # plt.xlim(6.25, 6.60)  # Limits for x-axis
        # plt.ylim(-0.45, 0.45)  # Limits for y-axis    
        plt.show() 
        
    if ip.flag_output == 1:
#-----------------------------------------------------------------------------    
        # End the timer
        end_time = time.time()
    
        # Calculate the elapsed time
        elapsed_time = end_time - start_time  
    
        # Display the execution time
        print(f"Output. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------
                
        with open(f"./output/R_lines_theta_r_box{ip.path}.pkl", "wb") as file:
            pickle.dump(R_lines_theta_r_box, file)
        with open(f"./output/Z_lines_theta_r_box{ip.path}.pkl", "wb") as file:
            pickle.dump(Z_lines_theta_r_box, file)
        with open(f"./output/B_lines_theta_r_box{ip.path}.pkl", "wb") as file:
            pickle.dump(B_lines_theta_r_box, file)
                                
                
#-----------------------------------------------------------------------------    
# End the timer
end_time = time.time()
        
# Calculate the elapsed time
elapsed_time = end_time - start_time  
        
# Display the execution time
print(f"End. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------         