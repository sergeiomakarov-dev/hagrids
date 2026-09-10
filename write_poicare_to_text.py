import numpy as np
import input_params as ip
import pickle
from run_time_indepndent_magnetic_hamiltonian import ( 
R_lines, Z_lines, theta_lines)

if ip.flag_output == 0:

    for i_poincare in range(ip.n_poincare):
           
        # File name
        filename = f"poincare_map_{i_poincare}.dat"
        
        with open(filename, 'w') as file:
            header_top = """# TYPE poincare_map
    """        
            file.write(header_top)
        
        for i_lines in range(ip.n_lines_tot):
            header = f"""# P0{R_lines[i_poincare][i_lines][0]: .10E}{Z_lines[i_poincare][i_lines][0]: .10E}{theta_lines[i_poincare][i_lines][0]: .10E}
    # DIRECTION -1
    # PHIX {i_poincare: .10E} deg
    # NSYMMETRY {ip.nfp}
    # POINTS {len(R_lines[0][0])}
    """
        
            data = [
                R_lines[i_poincare][i_lines][:],
                Z_lines[i_poincare][i_lines][:],
                theta_lines[i_poincare][i_lines][:],
                np.zeros(len(R_lines[0][0])),
            ]
            
    
            data_transp = [list(row) for row in zip(*data)]
            # Open the file in write mode
            with open(filename, 'a') as file:
                # Write the header
                file.write(header)
                
                # Write the data
                for row in data_transp:
                    file.write(f" {row[0]: .10E} {row[1]: .10E} {row[2]: .10E} {row[3]: .10E}\n")
            
        print(f"Data has been written to {filename}")

elif ip.flag_output == 1:
    with open("R_lines.pkl", "wb") as file:
        pickle.dump(R_lines, file)
    with open("Z_lines.pkl", "wb") as file:
        pickle.dump(Z_lines, file)