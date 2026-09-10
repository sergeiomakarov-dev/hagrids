import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.cm as cm
import input_params as ip

def plot_poincare (R_lines, Z_lines, n_lines_tot):
    
    if ip.flag_poincare :
        for i_poincare in range(ip.n_poincare):
            plt.figure()
            # Plot
            if ip.flag_analyt:
                for i_lines in range(n_lines_tot):
                    plt.scatter(R_lines[i_poincare][i_lines], Z_lines[i_poincare][i_lines], marker='.', s=ip.s_size, color='red')
                            
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