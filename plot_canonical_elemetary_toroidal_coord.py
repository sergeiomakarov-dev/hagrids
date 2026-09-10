import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import input_params as ip
import fun_time_independent_mh as fi

n_skip = 10

n_r = 4 * n_skip + 1
n_theta = 8 * n_skip + 1
psi_t_max = 5

r_max = np.sqrt(psi_t_max / (ip.B_0 * np.pi))


r_arr = np.linspace(0, r_max, n_r) 
psi_arr = r_arr**2 * (ip.B_0 * np.pi)
theta_arr = np.linspace(-np.pi, np.pi, n_theta) 

R_rad_arr_theta = [0 for _ in range(n_r)]
Z_rad_arr_theta = [0 for _ in range(n_r)]

R_pol_arr_theta = [np.zeros(n_r) for _ in range(n_theta)]
Z_pol_arr_theta = [np.zeros(n_r)  for _ in range(n_theta)]

for i_rad in range(n_r):
    R_rad_arr_theta[i_rad], Z_rad_arr_theta[i_rad] = fi.psi_theta_to_R_Z(psi_arr[i_rad], theta_arr, flag_vartheta = 1)
    
for i_pol in range(n_theta):
    for i_rad in range(n_r):
       R_pol_arr_theta[i_pol][i_rad], Z_pol_arr_theta[i_pol][i_rad] = fi.psi_theta_to_R_Z(psi_arr[i_rad], theta_arr[i_pol], flag_vartheta = 1)

R_rad_arr_vartheta = [0 for _ in range(n_r)]
Z_rad_arr_vartheta = [0 for _ in range(n_r)]

R_pol_arr_vartheta = [np.zeros(n_r) for _ in range(n_theta)]
Z_pol_arr_vartheta = [np.zeros(n_r)  for _ in range(n_theta)]

for i_rad in range(n_r):
    R_rad_arr_vartheta[i_rad], Z_rad_arr_vartheta[i_rad] = fi.psi_theta_to_R_Z(psi_arr[i_rad], theta_arr, flag_vartheta = 0)
    
for i_pol in range(n_theta):
    for i_rad in range(n_r):
       R_pol_arr_vartheta[i_pol][i_rad], Z_pol_arr_vartheta[i_pol][i_rad] = fi.psi_theta_to_R_Z(psi_arr[i_rad], theta_arr[i_pol], flag_vartheta = 0)
       
plt.figure()
# Plot

for i_rad in range(0, n_r, n_skip):
    plt.plot(R_rad_arr_vartheta[i_rad][n_skip:n_theta], Z_rad_arr_vartheta[i_rad][n_skip:n_theta], color='red')
    
for i_pol in range(0, n_theta, n_skip):
    plt.plot(R_pol_arr_vartheta[i_pol], Z_pol_arr_vartheta[i_pol], color='red')   
      
v_scale = 5
i_theta = 5 * n_skip 
# Plot the vector
plt.quiver(
    R_pol_arr_vartheta[i_theta][n_r-1], Z_pol_arr_vartheta[i_theta][n_r-1],  # Starting point
    v_scale * (R_pol_arr_vartheta[i_theta][n_r-1]-R_pol_arr_vartheta[i_theta][n_r-2]), 
    v_scale * (Z_pol_arr_vartheta[i_theta][n_r-1]-Z_pol_arr_vartheta[i_theta][n_r-2]),           # Components of the vector
    angles='xy', scale_units='xy', scale=1, color='red',
    headwidth=5, headlength=5,      # Enlarged arrowhead
    width=0.003                     # Thinner arrow shaft
)

shiftx = 0.1
shifty = -0.05
# Add text to label the vector
plt.text(
    R_pol_arr_vartheta[i_theta][n_r-1] + shiftx, Z_pol_arr_vartheta[i_theta][n_r-1] + shifty,  # Position of the text
    r"$r$",  # LaTeX text
    fontsize=12, color='red'
)

v_scale = 3
# Plot the vector
plt.quiver(
    R_rad_arr_theta[n_r-1][n_theta-1], Z_rad_arr_theta[n_r-1][n_theta-1],  # Starting point
    v_scale * (R_rad_arr_vartheta[n_r-1][n_theta-1]-R_rad_arr_vartheta[n_r-1][n_theta-2]), 
    v_scale * (Z_rad_arr_vartheta[n_r-1][n_theta-1]-Z_rad_arr_vartheta[n_r-1][n_theta-2]),           # Components of the vector
    angles='xy', scale_units='xy', scale=1, color='red',
    headwidth=5, headlength=5,      # Enlarged arrowhead
    width=0.003                     # Thinner arrow shaft
)

shiftx = 0.05
shifty = -0.3
# Add text to label the vector
plt.text(
    R_rad_arr_theta[n_r-1][n_theta-1] + shiftx, Z_rad_arr_theta[n_r-1][n_theta-1] + shifty,  # Position of the text
    r"$\vartheta$",  # LaTeX text
    fontsize=12, color='red'
)




for i_rad in range(0, n_r, n_skip):
    plt.plot(R_rad_arr_theta[i_rad][n_skip:n_theta], Z_rad_arr_theta[i_rad][n_skip:n_theta], color='blue')
    
for i_pol in range(0, n_theta, n_skip):
    plt.plot(R_pol_arr_theta[i_pol], Z_pol_arr_theta[i_pol], color='blue') 
    
v_scale = 5
i_theta = 5 * n_skip 
# Plot the vector
plt.quiver(
    R_pol_arr_theta[i_theta][n_r-1], Z_pol_arr_theta[i_theta][n_r-1],  # Starting point
    v_scale * (R_pol_arr_theta[i_theta][n_r-1]-R_pol_arr_theta[i_theta][n_r-2]), 
    v_scale * (Z_pol_arr_theta[i_theta][n_r-1]-Z_pol_arr_theta[i_theta][n_r-2]),           # Components of the vector
    angles='xy', scale_units='xy', scale=1, color='blue',
    headwidth=5, headlength=5,      # Enlarged arrowhead
    width=0.003                     # Thinner arrow shaft
)

shiftx = 0.1
shifty = 0.1
# Add text to label the vector
plt.text(
    R_pol_arr_theta[i_theta][n_r-1] + shiftx, Z_pol_arr_theta[i_theta][n_r-1] + shifty,  # Position of the text
    r"$\psi$",  # LaTeX text
    fontsize=12, color='blue'
)

v_scale = 2
# Plot the vector
plt.quiver(
    R_rad_arr_theta[n_r-1][n_theta-1], Z_rad_arr_theta[n_r-1][n_theta-1],  # Starting point
    v_scale * (R_rad_arr_vartheta[n_r-1][n_theta-1]-R_rad_arr_vartheta[n_r-1][n_theta-2]), 
    v_scale * (Z_rad_arr_vartheta[n_r-1][n_theta-1]-Z_rad_arr_vartheta[n_r-1][n_theta-2]),           # Components of the vector
    angles='xy', scale_units='xy', scale=1, color='blue',
    headwidth=5, headlength=5,      # Enlarged arrowhead
    width=0.003                     # Thinner arrow shaft
)

shiftx = 0.05
shifty = -0.15
# Add text to label the vector
plt.text(
    R_rad_arr_theta[n_r-1][n_theta-1] + shiftx, Z_rad_arr_theta[n_r-1][n_theta-1] + shifty,  # Position of the text
    r"$\theta$",  # LaTeX text
    fontsize=12, color='blue'
)

v_scale = 15
i_theta = 6 * n_skip 
# Plot the vector
plt.quiver(
    R_pol_arr_theta[i_theta][n_r-1], Z_pol_arr_theta[i_theta][n_r-1],  # Starting point
    v_scale * (R_pol_arr_theta[i_theta][n_r-1]-R_pol_arr_theta[i_theta][n_r-2]), 
    v_scale * (Z_pol_arr_theta[i_theta][n_r-1]-Z_pol_arr_theta[i_theta][n_r-2]),           # Components of the vector
    angles='xy', scale_units='xy', scale=1, color='k',
    linestyle='--',
    headwidth=5, headlength=5,      # Enlarged arrowhead
    width=0.003                     # Thinner arrow shaft
)

shiftx = -0.2
shifty = 0.15
# Add text to label the vector
plt.text(
    R_pol_arr_theta[i_theta][n_r-1] + shiftx, Z_pol_arr_theta[i_theta][n_r-1] + shifty,  # Position of the text
    r"$\vec{e}_\psi$",  # LaTeX text
    fontsize=12, color='k'
)

v_scale = 0.39
i_theta = 6 * n_skip 
# Plot the vector
plt.quiver(
    R_pol_arr_theta[i_theta][n_r-1], Z_pol_arr_theta[i_theta][n_r-1],  # Starting point
    v_scale * (R_pol_arr_theta[i_theta][n_r-1]-ip.R_0), 
    v_scale * (Z_pol_arr_theta[i_theta][n_r-1]),           # Components of the vector
    angles='xy', scale_units='xy', scale=1, color='k',
    linestyle='--',
    headwidth=5, headlength=5,      # Enlarged arrowhead
    width=0.003                     # Thinner arrow shaft
)

shiftx = 0.05
shifty = 0.15
# Add text to label the vector
plt.text(
    R_pol_arr_theta[i_theta][n_r-1] + shiftx, Z_pol_arr_theta[i_theta][n_r-1] + shifty,  # Position of the text
    r"$\vec{e}^\psi=\nabla\psi$",  # LaTeX text
    fontsize=12, color='k'
)

title = r""
plt.axis('equal') 
plt.title(title, fontsize=14)
plt.xlim(ip.R_min, ip.R_max)
plt.ylim(ip.Z_min - 1.0, ip.Z_max + 1.0)

plt.xlabel("R, [m]")
plt.ylabel("Z, [m]")
plt.show()

