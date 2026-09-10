import numpy as np
import matplotlib.pyplot as plt


# Load data
island = np.load('island_data.npz')
gamma_island = island['gamma_island']
p_island = island['p_island']

core = np.load('core_data.npz')
gamma_core = core['gamma_core']
p_core = core['p_core']

pfr = np.load('pfr_data.npz')
gamma_pfr = pfr['gamma_pfr']
p_pfr = pfr['p_pfr']

# Plot
plt.figure()
plt.scatter(gamma_core,   p_core,   marker='.', s=5, color='tab:blue',  label='Core')
plt.scatter(gamma_island, p_island, marker='.', s=5, color='tab:red',   label='Island')
plt.scatter(gamma_pfr,    p_pfr,    marker='.', s=5, color='tab:green', label='PFR')

# Highlight the first point from the island data
plt.scatter(gamma_island[0], p_island[0], marker='o', s=10, color='black', label='Separatrix')
plt.scatter(gamma_island[-1], p_island[-1], marker='o', s=10, color='tab:brown', label='O-point boudary')
plt.scatter(gamma_core[0], p_core[0], marker='o', s=10, color='tab:purple', label='Core boudary')
plt.scatter(gamma_pfr[-1], p_pfr[-1], marker='o', s=10, color='tab:orange', label='Core boudary')


# Local font size adjustments
plt.xlabel(r'$\gamma$', fontsize=16)        # Axis label font size
plt.ylabel(r'$p_\gamma$', fontsize=16)      # Axis label font size
plt.title(r'Phase space: Core, Island, and PFR grid nodes', fontsize=18)  # Title font size
#plt.legend(fontsize=14)                      # Legend font size
plt.xticks(fontsize=14)                      # X-axis tick font size
plt.yticks(fontsize=14)                      # Y-axis tick font size

# Grid and layout
plt.grid(True)
plt.tight_layout()
plt.show()