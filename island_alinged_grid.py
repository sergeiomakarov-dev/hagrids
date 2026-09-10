import numpy as np
import input_params as ip
import create_points_in_island as cpi
import create_points_in_core as cpc
import create_points_in_pfr as cpp
import fun_time_independent_mh as fi
import fun_time_dependent_mh as fd
import matplotlib.pyplot as plt
import pickle
import fun_plot as fp
import time

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

# Br_lines_theta_r_pfr = []
# Bvartheta_lines_theta_r_pfr = []
# Bvarphi_lines_theta_r_pfr = []
# BR_lines_theta_r_pfr = []
BZ_lines_theta_r_pfr = []
BY_lines_theta_r_pfr = []
BX_lines_theta_r_pfr = []
B_lines_theta_r_pfr = []

if ip.flag_type_init == 3: 
    
#-----------------------------------------------------------------------------    
    # End the timer
    end_time = time.time()
    
    # Calculate the elapsed time
    elapsed_time = end_time - start_time  
    
    # Display the execution time
    print(f"Tracing for islands. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------  
    
    n_lines_tot_island = ip.n_E * (4*ip.n_kin-3)
    n_lines_tot_core = ip.n_E_core *  (ip.m * (2*ip.n_kin-2)+1)
    n_lines_tot_pfr = ip.n_E_pfr *  (ip.m * (2*ip.n_kin-2)+1)
    
    for m_period_tmp in range(ip.m):
        
        alpha_arr_temp = cpi.alpha_arr_isl[m_period_tmp]
        psi_arr_temp = cpi.psi_arr_isl[m_period_tmp]
        
        alpha_0 = alpha_arr_temp.flatten(order='F')   
        psi_t_0 = psi_arr_temp.flatten(order='F')
        
        # Trace in the positive direction 
#-----------------------------------------------------------------------------    
        # End the timer
        end_time = time.time()
    
        # Calculate the elapsed time
        elapsed_time = end_time - start_time  
    
        # Display the execution time
        print(f"Island No {m_period_tmp}. Tracing in the positive direction. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------  
        if ip.flag_analyt or ip.flag_Stoermer_Verlet or ip.flag_Stoermer_Verlet_mod or ip.flag_Yoshida:
            R_lines, Z_lines, theta_lines, phi = fi.make_poincare (phi_span_pos, alpha_0, psi_t_0, n_lines_tot_island)
        elif ip.flag_Stoermer_Verlet_td or ip.flag_scipy_num_td or ip.flag_Yoshida_td:
            R_lines, Z_lines, theta_lines, phi = fd.make_poincare (phi_span_pos, alpha_0, psi_t_0, n_lines_tot_island)
        
        R_lines_theta_r_pos, Z_lines_theta_r_pos  = fi.structure_R_Z_arrays(R_lines, Z_lines, ip.n_E)
        phi_pos = phi
        
        if ip.flag_sym == 1:
            # Trace in the negative direction     
#-----------------------------------------------------------------------------    
            # End the timer
            end_time = time.time()
    
            # Calculate the elapsed time
            elapsed_time = end_time - start_time  
    
            # Display the execution time
            print(f"Island No {m_period_tmp}. Tracing in the negative direction. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------      
            if ip.flag_analyt or ip.flag_Stoermer_Verlet or ip.flag_Stoermer_Verlet_mod or ip.flag_Yoshida:        
                R_lines, Z_lines, theta_lines, phi = fi.make_poincare (phi_span_neg, alpha_0, psi_t_0, n_lines_tot_island)
            elif ip.flag_Stoermer_Verlet_td or ip.flag_scipy_num_td or ip.flag_Yoshida_td:
                R_lines, Z_lines, theta_lines, phi = fd.make_poincare (phi_span_neg, alpha_0, psi_t_0, n_lines_tot_island)   
                
            R_lines_theta_r_neg, Z_lines_theta_r_neg = fi.structure_R_Z_arrays(R_lines, Z_lines, ip.n_E)
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
        
        R_lines_island[m_period_tmp]=R_lines
        Z_lines_island[m_period_tmp]=Z_lines
        
        R_lines_theta_r_island[m_period_tmp]= R_lines_theta_r
        Z_lines_theta_r_island[m_period_tmp]= Z_lines_theta_r
        #fp.plot_poincare (R_lines, Z_lines)
        
        # Calculate the magentic feild 
        B_lines_theta_r_island[m_period_tmp], BX_lines_theta_r_island[m_period_tmp], BY_lines_theta_r_island[m_period_tmp], BZ_lines_theta_r_island[m_period_tmp], _, _, _, _ = fi.mag_field_in_pol_planes(R_lines_theta_r, Z_lines_theta_r, phi_theta_r)
                
    alpha_arr_temp = cpc.alpha_arr_core
    psi_arr_temp = cpc.psi_arr_core
            
    alpha_0 = alpha_arr_temp.flatten(order='F')   
    psi_t_0 = psi_arr_temp.flatten(order='F')
        
    # Trace in the positive direction 
    
#-----------------------------------------------------------------------------    
    # End the timer
    end_time = time.time()
    
    # Calculate the elapsed time
    elapsed_time = end_time - start_time  
    
    # Display the execution time
    print(f"Core region. Tracing in the positive direction. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------   
    if ip.flag_analyt or ip.flag_Stoermer_Verlet or ip.flag_Stoermer_Verlet_mod or ip.flag_Yoshida: 
        R_lines, Z_lines, theta_lines, phi = fi.make_poincare (phi_span_pos, alpha_0, psi_t_0, n_lines_tot_core)
    elif ip.flag_Stoermer_Verlet_td or ip.flag_scipy_num_td or ip.flag_Yoshida_td:  
        R_lines, Z_lines, theta_lines, phi = fd.make_poincare (phi_span_pos, alpha_0, psi_t_0, n_lines_tot_core) 
    
    R_lines_theta_r_pos, Z_lines_theta_r_pos = fi.structure_R_Z_arrays(R_lines, Z_lines, ip.n_E_core)
    phi_pos = phi
    
    if ip.flag_sym == 1:
        # Trace in the negative direction 
#-----------------------------------------------------------------------------    
        # End the timer
        end_time = time.time()
        
        # Calculate the elapsed time
        elapsed_time = end_time - start_time  
        
        # Display the execution time
        print(f"Core region. Tracing in the negative direction. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------  
        if ip.flag_analyt or ip.flag_Stoermer_Verlet or ip.flag_Stoermer_Verlet_mod or ip.flag_Yoshida:             
            R_lines, Z_lines, theta_lines, phi = fi.make_poincare (phi_span_neg, alpha_0, psi_t_0, n_lines_tot_core)
        elif ip.flag_Stoermer_Verlet_td or ip.flag_scipy_num_td or ip.flag_Yoshida_td:
            R_lines, Z_lines, theta_lines, phi = fd.make_poincare (phi_span_neg, alpha_0, psi_t_0, n_lines_tot_core)
            
        R_lines_theta_r_neg, Z_lines_theta_r_neg = fi.structure_R_Z_arrays(R_lines, Z_lines, ip.n_E_core)
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
    
    R_lines_core=R_lines
    Z_lines_core=Z_lines

    R_lines_theta_r_core=R_lines_theta_r
    Z_lines_theta_r_core=Z_lines_theta_r      
    
    #fp.plot_poincare (R_lines, Z_lines, n_lines_tot_core)
    
    # Calculate the magentic feild 
    B_lines_theta_r_core, BX_lines_theta_r_core, BY_lines_theta_r_core, BZ_lines_theta_r_core, _, _, _, _ = fi.mag_field_in_pol_planes(R_lines_theta_r, Z_lines_theta_r, phi_theta_r)
                
    alpha_arr_temp = cpp.alpha_arr_pfr
    psi_arr_temp = cpp.psi_arr_pfr
            
    alpha_0 = alpha_arr_temp.flatten(order='F')   
    psi_t_0 = psi_arr_temp.flatten(order='F')
        
    # Trace in the positive direction 
#-----------------------------------------------------------------------------    
    # End the timer
    end_time = time.time()
    
    # Calculate the elapsed time
    elapsed_time = end_time - start_time  
    
    # Display the execution time
    print(f"PFR region. Tracing in the positive direction. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------     
    if ip.flag_analyt or ip.flag_Stoermer_Verlet or ip.flag_Stoermer_Verlet_mod or ip.flag_Yoshida:             
        R_lines, Z_lines, theta_lines, phi = fi.make_poincare (phi_span_pos, alpha_0, psi_t_0, n_lines_tot_pfr) 
    elif ip.flag_Stoermer_Verlet_td or ip.flag_scipy_num_td or ip.flag_Yoshida_td:      
        R_lines, Z_lines, theta_lines, phi = fd.make_poincare (phi_span_pos, alpha_0, psi_t_0, n_lines_tot_pfr) 
        
    R_lines_theta_r_pos, Z_lines_theta_r_pos  = fi.structure_R_Z_arrays(R_lines, Z_lines, ip.n_E_pfr)
    phi_pos = phi
           
    if ip.flag_sym == 1:
        # Trace in the negative direction 
        
#-----------------------------------------------------------------------------    
        # End the timer
        end_time = time.time()
        
        # Calculate the elapsed time
        elapsed_time = end_time - start_time  
        
        # Display the execution time
        print(f"PFR region. Tracing in the negative direction. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------  
        if ip.flag_analyt or ip.flag_Stoermer_Verlet or ip.flag_Stoermer_Verlet_mod or ip.flag_Yoshida:                         
            R_lines, Z_lines, theta_lines, phi = fi.make_poincare (phi_span_neg, alpha_0, psi_t_0, n_lines_tot_pfr)
        elif ip.flag_Stoermer_Verlet_td or ip.flag_scipy_num_td or ip.flag_Yoshida_td: 
            R_lines, Z_lines, theta_lines, phi = fd.make_poincare (phi_span_neg, alpha_0, psi_t_0, n_lines_tot_pfr)
            
        R_lines_theta_r_neg, Z_lines_theta_r_neg = fi.structure_R_Z_arrays(R_lines, Z_lines, ip.n_E_pfr)
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
    
    R_lines_pfr=R_lines
    Z_lines_pfr=Z_lines
    
    R_lines_theta_r_pfr=R_lines_theta_r
    Z_lines_theta_r_pfr=Z_lines_theta_r 
    
    #fp.plot_poincare (R_lines, Z_lines, n_lines_tot_core)    
    
    # Calculate the magentic feild 
    B_lines_theta_r_pfr, BX_lines_theta_r_pfr, BY_lines_theta_r_pfr, BZ_lines_theta_r_pfr, _, _, _, _ = fi.mag_field_in_pol_planes(R_lines_theta_r, Z_lines_theta_r, phi_theta_r)
    
    m_period_tmp = 0
    for irad in range(ip.n_E-1):
        R_lines_tmp = R_lines_island[m_period_tmp][0][0:ip.n_E]
        Z_lines_tmp = Z_lines_island[m_period_tmp][0][0:ip.n_E]
        size_r_island_X  = np.zeros(len(R_lines_tmp)-1)
        length_r_island_X  = np.zeros(len(R_lines_tmp)-1)
        size_r_island_X, length_r_island_X, length_r_island_X_tot = fi.cell_size(R_lines_tmp, Z_lines_tmp)
        
        R_lines_tmp = R_lines_island[m_period_tmp][0][ip.n_E * (ip.n_kin-1):ip.n_E * ip.n_kin]
        Z_lines_tmp = Z_lines_island[m_period_tmp][0][ip.n_E * (ip.n_kin-1):ip.n_E * ip.n_kin]
        size_r_island_O  = np.zeros(len(R_lines_tmp)-1)
        length_r_island_O  = np.zeros(len(R_lines_tmp)-1)
        size_r_island_O, length_r_island_O, length_r_island_O_tot = fi.cell_size(R_lines_tmp, Z_lines_tmp)
    
    for irad in range(ip.n_E_core-1):
        R_lines_tmp = R_lines_core[0][0:ip.n_E_core]
        Z_lines_tmp = Z_lines_core[0][0:ip.n_E_core]
        size_r_core_X  = np.zeros(len(R_lines_tmp)-1)
        length_r_core_X  = np.zeros(len(R_lines_tmp)-1)
        size_r_core_X, length_r_core_X, length_r_core_X_tot = fi.cell_size(R_lines_tmp, Z_lines_tmp)
        
        R_lines_tmp = R_lines_core[0][ip.n_E_core * (ip.n_kin-1):ip.n_E_core * ip.n_kin]
        Z_lines_tmp = Z_lines_core[0][ip.n_E_core * (ip.n_kin-1):ip.n_E_core * ip.n_kin]
        size_r_core_O  = np.zeros(len(R_lines_tmp)-1)
        length_r_core_O  = np.zeros(len(R_lines_tmp)-1)
        size_r_core_O, length_r_core_O, length_r_core_O_tot = fi.cell_size(R_lines_tmp, Z_lines_tmp)
        
    for irad in range(ip.n_E_pfr-1):
        R_lines_tmp = R_lines_pfr[0][0:ip.n_E_pfr]
        Z_lines_tmp = Z_lines_pfr[0][0:ip.n_E_pfr]
        size_r_pfr_X  = np.zeros(len(R_lines_tmp)-1)
        length_r_pfr_X  = np.zeros(len(R_lines_tmp)-1)
        size_r_pfr_X, length_r_pfr_X, length_r_pfr_X_tot = fi.cell_size(R_lines_tmp, Z_lines_tmp)
        
        R_lines_tmp = R_lines_pfr[0][ip.n_E_pfr * (ip.n_kin-1):ip.n_E_pfr * ip.n_kin]
        Z_lines_tmp = Z_lines_pfr[0][ip.n_E_pfr * (ip.n_kin-1):ip.n_E_pfr * ip.n_kin]
        size_r_pfr_O  = np.zeros(len(R_lines_tmp)-1)
        length_r_pfr_O  = np.zeros(len(R_lines_tmp)-1)
        size_r_pfr_O, length_r_pfr_O, length_r_pfr_O_tot = fi.cell_size(R_lines_tmp, Z_lines_tmp)
        
    m_period_tmp = 0
    for irad in range(ip.n_E-1):
        R_lines_tmp = R_lines_island[m_period_tmp][0][ip.n_E-1:ip.n_E-1+ip.n_kin*ip.n_E:ip.n_E]
        Z_lines_tmp = Z_lines_island[m_period_tmp][0][ip.n_E-1:ip.n_E-1+ip.n_kin*ip.n_E:ip.n_E]
        size_p_island_O  = np.zeros(len(R_lines_tmp)-1)
        length_p_island_O  = np.zeros(len(R_lines_tmp)-1)
        size_p_island_O, length_p_island_O, length_p_island_O_tot = fi.cell_size(R_lines_tmp, Z_lines_tmp)
        
    for irad in range(ip.n_E_core-1):
        R_lines_tmp = R_lines_core[0][0:ip.n_kin*ip.n_E_core:ip.n_E_core]
        Z_lines_tmp = Z_lines_core[0][0:ip.n_kin*ip.n_E_core:ip.n_E_core]
        size_p_core_inner  = np.zeros(len(R_lines_tmp)-1)
        length_p_core_inner  = np.zeros(len(R_lines_tmp)-1)
        size_p_core_inner, length_p_core_inner, length_p_core_inner_tot = fi.cell_size(R_lines_tmp, Z_lines_tmp)
        
        R_lines_tmp = R_lines_core[0][ip.n_E_core-1:ip.n_E_core-1+ip.n_kin*ip.n_E_core:ip.n_E_core]
        Z_lines_tmp = Z_lines_core[0][ip.n_E_core-1:ip.n_E_core-1+ip.n_kin*ip.n_E_core:ip.n_E_core]
        size_p_core_sep  = np.zeros(len(R_lines_tmp)-1)
        length_p_core_sep  = np.zeros(len(R_lines_tmp)-1)
        size_p_core_sep, length_p_core_sep, length_p_core_sep_tot = fi.cell_size(R_lines_tmp, Z_lines_tmp)
        
    for irad in range(ip.n_E_pfr-1):
        R_lines_tmp = R_lines_pfr[0][0:ip.n_kin*ip.n_E_pfr:ip.n_E_pfr]
        Z_lines_tmp = Z_lines_pfr[0][0:ip.n_kin*ip.n_E_pfr:ip.n_E_pfr]
        size_p_pfr_sep  = np.zeros(len(R_lines_tmp)-1)
        length_p_pfr_sep  = np.zeros(len(R_lines_tmp)-1)
        size_p_pfr_sep, length_p_pfr_sep, length_p_pfr_sep_tot = fi.cell_size(R_lines_tmp, Z_lines_tmp)
        
        R_lines_tmp = R_lines_pfr[0][ip.n_E_pfr-1:ip.n_E_pfr-1+ip.n_kin*ip.n_E_pfr:ip.n_E_pfr]
        Z_lines_tmp = Z_lines_pfr[0][ip.n_E_pfr-1:ip.n_E_pfr-1+ip.n_kin*ip.n_E_pfr:ip.n_E_pfr]
        size_p_pfr_outer  = np.zeros(len(R_lines_tmp)-1)
        length_p_pfr_outer  = np.zeros(len(R_lines_tmp)-1)
        size_p_pfr_outer, length_p_pfr_outer, length_p_pfr_outer_tot = fi.cell_size(R_lines_tmp, Z_lines_tmp)

            
    if ip.flag_plot_grid or ip.flag_plot_cell_size or ip.flag_poincare:
#-----------------------------------------------------------------------------    
        # End the timer
        end_time = time.time()
            
        # Calculate the elapsed time
        elapsed_time = end_time - start_time  
            
        # Display the execution time
        print(f"Plotting. Execution time: {elapsed_time:.4f} seconds")
#----------------------------------------------------------------------------- 
    if ip.flag_plot_grid:
        plt.figure()
        for m_period_tmp in range(ip.m):
            #plt.scatter(R_lines_island[m_period_tmp][0], Z_lines_island[m_period_tmp][0], marker='.', s=ip.s_size, color='red')
            for ipol in range(4*ip.n_kin-3):
                plt.plot(R_lines_island[m_period_tmp][0][ipol*ip.n_E:(ipol+1)*ip.n_E], Z_lines_island[m_period_tmp][0][ipol*ip.n_E:(ipol+1)*ip.n_E], color='tab:blue', linewidth=1)
                if ipol != 4*ip.n_kin-3-1:
                    for irad in range(ip.n_E):
                        indices=[ipol*ip.n_E+irad,(ipol+1)*ip.n_E+irad]
                        plt.plot([R_lines_island[m_period_tmp][0][i] for i in indices], [Z_lines_island[m_period_tmp][0][i] for i in indices], color='tab:blue', linewidth=1)
                        
                    for irad in [ip.n_E-1]:
                        indices=[ipol*ip.n_E+irad,(ipol+1)*ip.n_E+irad]
                        plt.plot([R_lines_island[m_period_tmp][0][i] for i in indices], [Z_lines_island[m_period_tmp][0][i] for i in indices], color='tab:red')
                
        #plt.scatter(R_lines_core[0], Z_lines_core[0], marker='.', s=ip.s_size, color='red')
        for ipol in range(ip.m*(2*ip.n_kin-2)):
            plt.plot(R_lines_core[0][ipol*ip.n_E_core:(ipol+1)*ip.n_E_core], Z_lines_core[0][ipol*ip.n_E_core:(ipol+1)*ip.n_E_core], color='tab:blue', linewidth=1)
            if ipol != ip.m*(2*ip.n_kin-2):
                for irad in range(ip.n_E_core):
                    indices=[ipol*ip.n_E_core+irad,(ipol+1)*ip.n_E_core+irad]
                    plt.plot([R_lines_core[0][i] for i in indices], [Z_lines_core[0][i] for i in indices], color='tab:blue', linewidth=1)
                    
                for irad in [0,ip.n_E_core-1]:
                    indices=[ipol*ip.n_E_core+irad,(ipol+1)*ip.n_E_core+irad]
                    plt.plot([R_lines_core[0][i] for i in indices], [Z_lines_core[0][i] for i in indices], color='tab:red')

        #plt.scatter(R_lines_pfr[0], Z_lines_pfr[0], marker='.', s=ip.s_size, color='red')
        for ipol in range(ip.m*(2*ip.n_kin-2)):
            plt.plot(R_lines_pfr[0][ipol*ip.n_E_pfr:(ipol+1)*ip.n_E_pfr], Z_lines_pfr[0][ipol*ip.n_E_pfr:(ipol+1)*ip.n_E_pfr], color='tab:blue', linewidth=1)
            if ipol != ip.m*(2*ip.n_kin-2):
                for irad in range(ip.n_E_pfr):
                    indices=[ipol*ip.n_E_pfr+irad,(ipol+1)*ip.n_E_pfr+irad]
                    plt.plot([R_lines_pfr[0][i] for i in indices], [Z_lines_pfr[0][i] for i in indices], color='tab:blue', linewidth=1)
                    
                for irad in [0,ip.n_E_pfr-1]:
                    indices=[ipol*ip.n_E_pfr+irad,(ipol+1)*ip.n_E_pfr+irad]
                    plt.plot([R_lines_pfr[0][i] for i in indices], [Z_lines_pfr[0][i] for i in indices], color='tab:red')
                    
        plt.axis('equal')      
        plt.xlim(ip.R_min, ip.R_max)
        plt.ylim(ip.Z_min, ip.Z_max)                 
        plt.title(r'2D base mesh; multi-zone', fontsize=20)        
        plt.xlabel(r'R, m', fontsize=20)
        plt.ylabel(r'Z, m', fontsize=20)
        plt.tick_params(axis='both', labelsize=16) 
        plt.show()
        plt.tight_layout()

    if ip.flag_plot_cell_size:
        plt.figure() 
        plt.plot(length_r_island_X / length_r_island_X_tot, size_r_island_X / max(size_r_island_X),
                 color='black', label='X-point level') 
        plt.plot(length_r_island_O / length_r_island_O_tot, size_r_island_O / max(size_r_island_O),
                 color='red', label='O-point level')
                 
        plt.title(r'Radial cell size: island; $N = {:d}$; $\alpha = {:.1f}$'.format(ip.n_E,ip.alpha_island_rad), fontsize=14)        
        plt.xlabel(r'Normilised distance from separatrix to O-point', fontsize=14)
        plt.ylabel(r'Normilised radial cell size', fontsize=14)
        plt.legend( fontsize=14)
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.show() 
        
        plt.figure() 
        plt.plot(length_r_core_X / length_r_core_X_tot, size_r_core_X / max(size_r_core_X),
                 color='black', label='X-point level') 
        plt.plot(length_r_core_O / length_r_core_O_tot, size_r_core_O / max(size_r_core_O),
                 color='red', label='O-point level')
        plt.title(r'Radial cell size: core $N = {:d}$; $\alpha = {:.1f}$'.format(ip.n_E_core,ip.alpha_core_rad), fontsize=14)        
        plt.xlabel(r'Normilised distance from core boundary to separatrix', fontsize=14)
        plt.ylabel(r'Normilised radial cell size', fontsize=14)
        plt.legend( fontsize=14)
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.show() 
        
        plt.figure() 
        plt.plot(length_r_pfr_X / length_r_pfr_X_tot, size_r_pfr_X / max(size_r_pfr_X),
                 color='black', label='X-point level') 
        plt.plot(length_r_pfr_O / length_r_pfr_O_tot, size_r_pfr_O / max(size_r_pfr_O),
                 color='red', label='O-point level')
        plt.title(r'Radial cell size: PFR; $N = {:d}$; $\alpha = {:.1f}$'.format(ip.n_E_pfr,ip.alpha_pfr_rad), fontsize=14)        
        plt.xlabel(r'Normilised distance from separatrix to outer boundary', fontsize=14)
        plt.ylabel(r'Normilised radial cell size', fontsize=14)
        plt.legend( fontsize=14)
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.show() 
        
        plt.figure() 
        plt.plot(length_p_core_inner / length_p_core_inner_tot, size_p_core_inner / max(size_p_core_inner),
                 color='black', label='Core inner boudary') 
        plt.plot(length_p_core_sep / length_p_core_sep_tot, size_p_core_sep / max(size_p_core_sep),
                 color='red', label='Separatrix, core side')
        plt.plot(length_p_island_O/ length_p_island_O_tot, size_p_island_O / max(size_p_island_O),
                 color='green', label='O-point, core side')
        plt.plot(length_p_pfr_sep / length_p_pfr_sep_tot, size_p_pfr_sep / max(size_p_pfr_sep),
                 color='blue', label='Separatrix, PFR side')
        plt.plot(length_p_pfr_outer / length_p_pfr_outer_tot, size_p_pfr_outer / max(size_p_pfr_outer),
                 color='magenta', label='PFR outer boudary') 
        plt.title(r'Poloidal cell size; $N = {:d}$; $\alpha = {:.1f}$'.format(ip.n_kin,ip.alpha_pol), fontsize=14)        
        plt.xlabel(r'Normilised distance from X-point to O-point level', fontsize=14)
        plt.ylabel(r'Normilised poloidal cell size', fontsize=14)
        plt.legend( fontsize=14)
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.show()       

    if ip.flag_poincare :
        plt.figure()
        for i_poincare in range(len(R_lines_theta_r_core)):
            for m_period_tmp in range(ip.m):
                plt.scatter(R_lines_theta_r_island[m_period_tmp][i_poincare][:,1:ip.n_E], Z_lines_theta_r_island[m_period_tmp][i_poincare][:,1:ip.n_E], marker='.', s=ip.s_size, color='red')

            plt.scatter(R_lines_theta_r_core[i_poincare][:,0:ip.n_E_core-1], Z_lines_theta_r_core[i_poincare][:,0:ip.n_E_core-1], marker='.', s=ip.s_size, color='blue')     
            plt.scatter(R_lines_theta_r_pfr[i_poincare][:,1:], Z_lines_theta_r_pfr[i_poincare][:,1:], marker='.', s=ip.s_size, color='green')
                
        if ip.flag_analyt or ip.flag_Stoermer_Verlet or ip.flag_Stoermer_Verlet_mod or ip.flag_Yoshida:
            title = r"$A = {:.2e}$".format(ip.A)
        elif ip.flag_Stoermer_Verlet_td or ip.flag_scipy_num_td or ip.flag_Yoshida_td:
            title = r"$A_1 = {:.2e}$; $A_2 = {:.2e}$".format(ip.A1, ip.A2)
        else:
            title = r""        
        plt.title(title, fontsize=14)  
        plt.xlim(ip.R_min, ip.R_max)
        plt.ylim(ip.Z_min, ip.Z_max)
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
        for m_period_tmp in range(ip.m):
            with open(f"./output/R_lines_theta_r_island_{m_period_tmp}{ip.path}.pkl", "wb") as file:
                pickle.dump(R_lines_theta_r_island[m_period_tmp], file)
            with open(f"./output/Z_lines_theta_r_island_{m_period_tmp}{ip.path}.pkl", "wb") as file:
                pickle.dump(Z_lines_theta_r_island[m_period_tmp], file)
            with open(f"./output/B_lines_theta_r_island_{m_period_tmp}{ip.path}.pkl", "wb") as file:
                pickle.dump(B_lines_theta_r_island[m_period_tmp], file)

        with open(f"./output/R_lines_theta_r_core{ip.path}.pkl", "wb") as file:
            pickle.dump(R_lines_theta_r_core, file)
        with open(f"./output/Z_lines_theta_r_core{ip.path}.pkl", "wb") as file:
            pickle.dump(Z_lines_theta_r_core, file)
        with open(f"./output/B_lines_theta_r_core{ip.path}.pkl", "wb") as file:
            pickle.dump(B_lines_theta_r_core, file)
                
        with open(f"./output/R_lines_theta_r_pfr{ip.path}.pkl", "wb") as file:
            pickle.dump(R_lines_theta_r_pfr, file)
        with open(f"./output/Z_lines_theta_r_pfr{ip.path}.pkl", "wb") as file:
            pickle.dump(Z_lines_theta_r_pfr, file)
        with open(f"./output/B_lines_theta_r_pfr{ip.path}.pkl", "wb") as file:
            pickle.dump(B_lines_theta_r_pfr, file)
                                
                
#-----------------------------------------------------------------------------    
# End the timer
end_time = time.time()
        
# Calculate the elapsed time
elapsed_time = end_time - start_time  
        
# Display the execution time
print(f"End. Execution time: {elapsed_time:.4f} seconds")
#-----------------------------------------------------------------------------         