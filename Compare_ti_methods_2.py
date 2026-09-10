import numpy as np
import input_params as ip
import create_points_one_zone as cpo
import fun_time_independent_mh as fi
import fun_time_dependent_mh as fd
import matplotlib.pyplot as plt
import pickle
import fun_plot as fp
import time
import create_points_in_island as cpi

psi_t_0 = 1.1
alpha_0 = 0.00

end = ip.n_turns * np.pi / ip.nfp
ponts = 10

plot_from = 0

m_period_tmp=0

alpha_arr_temp = cpi.alpha_arr_isl[m_period_tmp]
psi_arr_temp = cpi.psi_arr_isl[m_period_tmp]
        
alpha_0 = alpha_arr_temp.flatten(order='F')   
psi_t_0 = psi_arr_temp.flatten(order='F')

phi_span = np.linspace(0, end * 100, ponts * 100)
phi_span_long = np.linspace(0, end * 100, ponts * 100)
phi_span_10 = np.linspace(0, end, ponts * 10)
phi_span_neg = np.linspace(0, -end, ponts)

#y1 = np.array([[psi_t_0],[alpha_0]])
#y_sv=fi.num_solver(phi_span, y1, flag_sv = 1)

#y_sv_mod=fi.num_solver(phi_span, y1, flag_sv_mod = 1)
plt.figure()
for i_lines in range(len(psi_t_0)):
    y1 = np.array([[psi_t_0[i_lines]],[alpha_0[i_lines]]])
    y_yo = fi.num_solver(phi_span, y1, flag_yo = 1)
    
    #y_an=fi.pendulum(phi_span, y1)
    
    #y_sn=fi.num_solver(phi_span, y1, flag_sn = 1, method = 'Radau')
    
    #y_aa=fi.sine_solution(phi_span, y1)
    
    
    
    #plt.plot(y_sv[1, plot_from :], y_sv[0, plot_from :], label='sv', color='blue', marker='o')
    #plt.plot(y_sv_mod[1, plot_from :], y_sv_mod[0, plot_from :], label='sv_mod', color='orange', marker='o')
    #plt.plot(y_an[1, plot_from :], y_an[0, plot_from :], label='an', color='red', marker='o')
    #plt.plot(y_sn[1, plot_from :], y_sn[0, plot_from :], label='sn', color='black', marker='o')
    #plt.plot(y_yo_long[1, plot_from :], y_yo_long[0, plot_from :], label='yo_long', color='green')
    plt.plot(y_yo[1, plot_from], y_yo[0, plot_from], label='yo', color='blue', marker='o')
    plt.plot(y_yo[1, plot_from :], y_yo[0, plot_from :], label='yo', color='blue')
    #plt.plot(y_yo_10 [1, plot_from :], y_yo_10 [0, plot_from :], label='yo x 10', color='blue')
    #plt.plot(-y_yo_neg[1, plot_from :], y_yo_neg[0, plot_from :], label='yo_neg', color='red')
    #plt.plot(y_aa[1, plot_from :], y_aa[0, plot_from :], label='aa', color='magenta', marker='o')
    
    #plt.legend()