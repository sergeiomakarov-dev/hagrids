#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.cm as cm
import input_params as ip
from run_time_indepndent_magnetic_hamiltonian import ( phi_span, 
R_lines, Z_lines,
R_lines_sv, Z_lines_sv,
R_lines_sn, Z_lines_sn,
R_iota_res_test, Z_iota_res_test)
from run_time_depndent_magnetic_hamiltonian import (psi_t_0,
R_lines_sv_td, Z_lines_sv_td,
R_lines_sn_td, Z_lines_sn_td,
iota_sv_td, iota_sv_td_max, iota_sv_td_min,
psi_test, R_test, Z_test, R_configs, Z_configs,
psi_0_configs,
R_iota_res1, Z_iota_res1,
R_iota_res2, Z_iota_res2)


if ip.flag_iota_profile :
    psi_span = np.linspace(ip.psi_t_min, ip.psi_t_max, 10)
    iota_span = ip.iota_a * psi_span + ip.iota_b
    
    plt.figure()
    # Plot
    plt.plot(psi_span, iota_span, color='black', label='iota')# label='Runge-Kutta')
    plt.title(r'Iota profile', fontsize=14)        
    plt.xlabel(r'$\psi$', fontsize=14)
    plt.ylabel(r'$\iota$', fontsize=14)
    plt.legend( fontsize=14)
    plt.show()    
    
if ip.flag_poincare :
    for i_poincare in range(ip.n_poincare):
        plt.figure()
        # Plot
        if ip.flag_analyt:
            for i_lines in range(ip.n_lines_tot):
                plt.scatter(R_lines[i_poincare][i_lines], Z_lines[i_poincare][i_lines], marker='.', s=ip.s_size, color='red')
            #plt.scatter(R_iota_res_test, Z_iota_res_test, marker='.', s=ip.s_size*50, color='blue')
    
        
        if ip.flag_Stoermer_Verlet:
            #plt.figure()
            # Plot
            for i_lines in range(ip.n_lines_tot):
                plt.scatter(R_lines_sv[i_lines], Z_lines_sv[i_lines], marker='.', s=ip.s_size, color='blue')
            
    
        if ip.flag_scipy_num:
            #plt.figure()
            # Plot
            for i_lines in range(ip.n_lines_tot):
                plt.scatter(R_lines_sn[i_lines], Z_lines_sn[i_lines], marker='.', s=ip.s_size, color='black')
                
        if ip.flag_Stoermer_Verlet_td:
            #plt.figure()
            # Plot
            for i_lines in range(ip.n_lines_tot_td):
                temp = 1-(psi_t_0[i_lines] - ip.psi_t_min)/(ip.psi_t_max - ip.psi_t_min)
                color = cm.plasma(temp)
                plt.scatter(R_lines_sv_td[i_lines], Z_lines_sv_td[i_lines], marker='.', s=ip.s_size, color=color)  
                
            plt.scatter(R_iota_res_test, Z_iota_res_test, marker='.', s=ip.s_size*50, color='red')
            
            #temp = 1-(psi_test - ip.psi_t_min)/(ip.psi_t_max - ip.psi_t_min)
            #color = cm.plasma(temp)
            #plt.scatter(R_test, Z_test-0.001, marker='.', s=ip.s_size*100, color=color)
            
            if ip.flag_iota_res_surf:
                plt.plot(R_iota_res1, Z_iota_res1, color='gray')
                plt.plot(R_iota_res2, Z_iota_res2, color='gray')
            
            color_bar = plt.colorbar(label=r'Initial $\psi$')
            plt.set_cmap('plasma')
            color_bar.set_ticks([0, 1])
            color_bar.set_ticklabels([ip.psi_t_max, ip.psi_t_min])
    
            
        if ip.flag_scipy_num_td:
            # Plot
            for i_lines in range(ip.n_lines_tot_td):
                plt.scatter(R_lines_sn_td[i_lines], Z_lines_sn_td[i_lines], marker='.', s=ip.s_size, color='black')
                
        #plt.axis('equal')  
        if ip.flag_analyt:
            title = r"$A = {:.2e}$".format(ip.A)
        elif ip.flag_Stoermer_Verlet_td:
            title = r"$A_1 = {:.2e}$; $A_2 = {:.2e}$".format(ip.A1, ip.A2)
        else:
            title = r""        
        plt.title(title, fontsize=14)  
        plt.xlim(ip.R_min, ip.R_max)
        plt.ylim(ip.Z_min, ip.Z_max)
        plt.show()  
            
if ip.flag_poincare_2 :
    plt.figure()
    # Plot
    if ip.flag_Stoermer_Verlet_td:
        #plt.figure()
        # Plot
        for i_lines in range(ip.n_lines_tot_td):
            for i_points in range(len(R_lines_sv_td[i_lines])):
                temp = 1-(iota_sv_td[i_lines][i_points] - iota_sv_td_min)/(iota_sv_td_max - iota_sv_td_min)
                color = cm.plasma(temp)
                #plt.scatter(R_lines_sv_td[i_lines][i_points], Z_lines_sv_td[i_lines][i_points], marker='.', s=ip.s_size, color=color)
                plt.plot(R_lines_sv_td[i_lines][i_points], Z_lines_sv_td[i_lines][i_points], marker='o', markersize=1, color=color)
                
        #color_bar = plt.colorbar(label=r'local iota')                
        #plt.set_cmap('plasma')
        #color_bar.set_ticks([0, 1])
        #color_bar.set_ticklabels([iota_sv_td_max, iota_sv_td_min])
    #plt.axis('equal')  
    if ip.flag_analyt:
        title = r"$A = {:.2e}$".format(ip.A)
    elif ip.flag_Stoermer_Verlet_td:
        title = r"$A_1 = {:.2e}$; $A_2 = {:.2e}$".format(ip.A1, ip.A2)
    else:
        title = r""
    plt.title(title, fontsize=14)  
    plt.xlim(ip.R_min, ip.R_max)
    plt.ylim(ip.Z_min, ip.Z_max)
    plt.show()    

if ip.flag_poincare_3 :
    plt.figure()
    # Plot
    if ip.flag_Stoermer_Verlet_td:
        #plt.figure()
        # Plot
        for i_lines in range(ip.n_lines_tot_td):
            if len(R_lines_sv_td[i_lines]) > 0:
                temp = 1-(iota_sv_td[i_lines][0] - iota_sv_td_min)/(iota_sv_td_max - iota_sv_td_min)
                color = cm.plasma(temp)
                plt.scatter(R_lines_sv_td[i_lines], Z_lines_sv_td[i_lines], marker='.', s=ip.s_size, color=color)
            #plt.plot(R_lines_sv_td[i_lines][i_points], Z_lines_sv_td[i_lines][i_points], marker='o', markersize=1, color=color)
                
        #color_bar = plt.colorbar(label=r'local iota')                
        #plt.set_cmap('plasma')
        #color_bar.set_ticks([0, 1])
        #color_bar.set_ticklabels([iota_sv_td_max, iota_sv_td_min])
    #plt.axis('equal')  
    if ip.flag_analyt:
        title = r"$A = {:.2e}$".format(ip.A)
    elif ip.flag_Stoermer_Verlet_td:
        title = r"$A_1 = {:.2e}$; $A_2 = {:.2e}$".format(ip.A1, ip.A2)
    else:
        title = r""
    plt.title(title, fontsize=14)  
    plt.xlim(ip.R_min, ip.R_max)
    plt.ylim(ip.Z_min, ip.Z_max)
    plt.show()    

# Function to update the plot for each frame
def update(frame):
    scat.set_offsets(np.column_stack((R_configs[frame], Z_configs[frame])))
    text.set_text(f'field periods = {frame}')  # Update text
    return scat,text

if ip.flag_plot_configs :
    #plt.figure()
    # Plot
            
    if ip.flag_Stoermer_Verlet_td:
        fig, ax = plt.subplots()
        # Plot
        for i_lines in range(ip.n_lines_tot_td):
            ax.scatter(R_lines_sv_td[i_lines], Z_lines_sv_td[i_lines], marker='.', s=1, color='gray')  
            
        if ip.flag_plot_configs_ani:
            temp = 1-(psi_0_configs - ip.psi_0_min_configs)/(ip.psi_0_max_configs - ip.psi_0_min_configs)
            color = cm.plasma(temp)
            scat = ax.scatter(R_configs[0], Z_configs[0], marker='.', s=ip.s_size*10, color=color) 
            text = ax.text( 0.9, 0.9, '', transform=ax.transAxes, ha='center')
            ani = FuncAnimation(fig, update, frames=len(R_configs), interval=1000)
        if ip.flag_plot_one_fl:
            for i_points in range(ip.n_turns_one_fl):
                temp = i_points/ip.n_turns_one_fl
                color = cm.plasma(temp)
                scat = ax.scatter(R_configs[i_points][0], Z_configs[i_points][0], marker='.', s=ip.s_size*10, color=color) 
        
        if ip.flag_analyt:
            title = r"$A = {:.2e}$".format(ip.A)
        elif ip.flag_Stoermer_Verlet_td:
            title = r"$A_1 = {:.2e}$; $A_2 = {:.2e}$".format(ip.A1, ip.A2)
        else:
            title = r""
        plt.title(title, fontsize=14)      
        plt.xlim(ip.R_min, ip.R_max)
        plt.ylim(ip.Z_min, ip.Z_max)
        plt.show()   
        
if ip.flag_plot_one_fl_1D_plot :
    plt.figure()
    # Plot
    plt.plot(R_configs[:][0], Z_configs[i_points][0], color='black', label='iota')# label='Runge-Kutta')
    plt.title(r'Iota profile', fontsize=14)        
    plt.xlabel(r'$\psi$', fontsize=14)
    plt.ylabel(r'$\iota$', fontsize=14)
    plt.legend( fontsize=14)
    plt.show()  