import numpy as np
import input_params as ip
import create_points_one_zone as cpo
import fun_time_independent_mh as fi
import fun_time_dependent_mh as fd
import matplotlib.pyplot as plt
import pickle
import fun_plot as fp
import time

psi_t_0 = 1.28
alpha_0 = 0.00

end = 1000
ponts = 1000

plot_from = 900

phi_span = np.linspace(0, end, ponts)

y1 = np.array([[psi_t_0],[alpha_0]])
y_sv=fi.num_solver(phi_span, y1, flag_sv = 1)

y_sv_mod=fi.num_solver(phi_span, y1, flag_sv_mod = 1)

y_yo = fi.num_solver(phi_span, y1, flag_yo = 1)

y_an=fi.pendulum(phi_span, y1)

y_sn=fi.num_solver(phi_span, y1, flag_sn = 1, method = 'Radau')

y_aa=fi.ellipse_solution(phi_span, y1)


plt.figure()
plt.plot(y_sv[1, plot_from :], y_sv[0, plot_from :], label='sv', color='blue', marker='o')
plt.plot(y_sv_mod[1, plot_from :], y_sv_mod[0, plot_from :], label='sv_mod', color='orange', marker='o')
plt.plot(y_an[1, plot_from :], y_an[0, plot_from :], label='an', color='red', marker='o')
plt.plot(y_sn[1, plot_from :], y_sn[0, plot_from :], label='sn', color='black', marker='o')
plt.plot(y_yo[1, plot_from :], y_yo[0, plot_from :], label='yo', color='green', marker='o')
plt.plot(y_aa[1, plot_from :], y_aa[0, plot_from :], label='aa', color='magenta', marker='o')

plt.legend()