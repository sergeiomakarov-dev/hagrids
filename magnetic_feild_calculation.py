import numpy as np

def Br_perturbation(r, vartheta, varphi, m, A, iota_res, R, varepsilon, psi_0, B0, k=1.0, flag_type_A_fun=0, flag_toroildal=1):
    """
    Calculate the radial component of the magnetic field perturbation (Br_pert).

    Parameters:
        r (float): Radial coordinate.
        vartheta (float): Poloidal angular coordinate (in radians).
        varphi (float): Toroidal angular coordinate (in radians).
        m (int): Mode number.
        A (float): Amplitude of perturbation.
        iota_res (float): Resonant iota parameter.
        R (float): Effective major radius.
        varepsilon (float): Inverse aspect ratio (r / R0).
        psi_0 (float): Resonant flux value.
        B0 (float): Base magnetic field.
        k (float): Scaling factor for perturbation.
        flag_type_A_fun (int): Type of perturbation function (0 or 1).
        flag_toroildal (int): 0 (uniform toroidal field) or 1 (geometry-modified field).

    Returns:
        float: Radial perturbation component (Br_pert).
    """
    base_factor = m * A / (2 * np.pi * r * R)
    if flag_toroildal == 0:
        phase = m * (vartheta - iota_res * varphi)
    
        if flag_type_A_fun == 0:
            Br_pert = base_factor * np.sin(phase)
        elif flag_type_A_fun == 1:
            floor_arg = 1 - k + k * (np.pi * r**2 * B0) / psi_0
            extended_factor = base_factor * floor_arg * np.heaviside(floor_arg, 1.0)
            Br_pert = extended_factor * np.sin(phase)
        else:
            raise ValueError("Invalid value for flag_type_A_fun. Must be 0 or 1.")
    elif flag_toroildal == 1:
        phase = m * (vartheta - varepsilon * np.sin(vartheta) - iota_res * varphi)
    
        if flag_type_A_fun == 0:
            Br_pert = base_factor * np.sin(phase) - varepsilon * np.cos(vartheta) * base_factor * np.sin(phase)
        elif flag_type_A_fun == 1:
            floor_arg = 1 - k + k * (np.pi * r**2 * B0) / psi_0
            extended_factor = base_factor * floor_arg * np.heaviside(floor_arg, 1.0)
            Br_pert = extended_factor * np.sin(phase) - varepsilon * np.cos(vartheta) * extended_factor * np.sin(phase)
        else:
            raise ValueError("Invalid value for flag_type_A_fun. Must be 0 or 1.")
    else:
        raise ValueError("Invalid value for flag_toroildal. Must be 0 or 1.")
    return Br_pert


def Bvartheta_perturbation(r, vartheta, varphi, m, A, iota_res, R, varepsilon, psi_0, B0, k=1.0, flag_type_A_fun=0, flag_toroildal=1):
    """
    Calculate the poloidal component of the magnetic field perturbation (Bvartheta_pert).

    Parameters:
        Same as `Br_perturbation`.

    Returns:
        float: Poloidal perturbation component (Bvartheta_pert).
    """
    base_factor = m * A / (2 * np.pi * r * R)
    if flag_toroildal == 0:
        phase = m * (vartheta - iota_res * varphi)
    
        if flag_type_A_fun == 0:
            Bvartheta_pert = 0
        elif flag_type_A_fun == 1:
            floor_arg = 1 - k + k * (np.pi * r**2 * B0) / psi_0
            step = np.heaviside(floor_arg, 1.0)
            additional_term = (A / (2 * np.pi * r * R)) * (k * (2 * np.pi * r * B0) / psi_0) * np.cos(phase) * step
            Bvartheta_pert = additional_term
        else:
            raise ValueError("Invalid value for flag_type_A_fun. Must be 0 or 1.")
    elif flag_toroildal == 1:
        phase = m * (vartheta - varepsilon * np.sin(vartheta) - iota_res * varphi)
    
        if flag_type_A_fun == 0:
            Bvartheta_pert = varepsilon * np.sin(vartheta) * base_factor * np.sin(phase)
        elif flag_type_A_fun == 1:
            floor_arg = 1 - k + k * (np.pi * r**2 * B0) / psi_0
            step = np.heaviside(floor_arg, 1.0)
            extended_factor = base_factor * floor_arg * step
            additional_term = (A / (2 * np.pi * r * R)) * (k * (2 * np.pi * r * B0) / psi_0) * np.cos(phase) * step
            Bvartheta_pert = varepsilon * np.sin(vartheta) * extended_factor * np.sin(phase) + additional_term
        else:
            raise ValueError("Invalid value for flag_type_A_fun. Must be 0 or 1.")
    else:
        raise ValueError("Invalid value for flag_toroildal. Must be 0 or 1.")
    return Bvartheta_pert


def magnetic_field_components(r, vartheta, varphi, B0, R0,
                              iota_a=0.3, iota_b=0.7,
                              m=5, m1=5, m2=5, m3=5,
                              iota_res=1, iota_res1=1, iota_res2=1, iota_res3=1,
                              A=0.0, A1=0.0, A2=0.0, A3=0.0, k=1.0,
                              A_arr=(), m_arr=(), iota_res_arr=(),
                              flag_form=0, flag_toroildal=1,
                              flag_type_A_fun=0):
    """
    Compute the magnetic field components (Br, Bvartheta, Bvarphi).

    Parameters:
        r (float): Radial coordinate.
        vartheta (float): Poloidal angular coordinate (in radians).
        varphi (float): Toroidal angular coordinate (in radians).
        B0 (float): Base magnetic field.
        R0 (float): Reference major radius.
        iota_a, iota_b (float): Coefficients for iota profile.
        m, m1, m2, m3 (int): Mode numbers for perturbations.
        iota_res, iota_res1, iota_res2, iota_res3 (float): Resonant iota parameters.
        A, A1, A2, A3 (float): Amplitudes of perturbations.
        k (float): Scaling factor for perturbations.
        A_arr, m_arr, iota_res_arr (sequence): Amplitudes, mode numbers and resonant iota of
            the additional perturbations (A_list, m_list, n_list of the input), used if flag_form = 2.
        flag_form (int): Perturbation set: 0 the single perturbation (m, A); 1 the perturbations
            (m1, A1), (m2, A2), (m3, A3); 2 the perturbation (m1, A1) plus the additional
            perturbations A_arr, m_arr, iota_res_arr.
        flag_toroildal (int): 0 (uniform toroidal field) or 1 (geometry-modified field).
        flag_type_A_fun (int): Type of perturbation function (0 or 1).

    Returns:
        tuple: Magnetic field components (Br, Bvartheta, Bvarphi).
    """
    if flag_form not in [0, 1, 2]:
        raise ValueError("Invalid value for flag_form. Must be 0, 1 or 2.")
    if flag_form == 2 and not (len(A_arr) == len(m_arr) == len(iota_res_arr)):
        raise ValueError("A_arr, m_arr and iota_res_arr must have the same length.")
    if flag_toroildal not in [0, 1]:
        raise ValueError("Invalid value for flag_toroildal. Must be 0 or 1.")

    varepsilon = r / R0
    R = R0 * (1 + varepsilon * np.cos(vartheta))
    iota_profile_r = iota_a * B0 * np.pi * r**2 + iota_b

    psi_0 = (iota_res - iota_b) / iota_a
    psi_01 = (iota_res1 - iota_b) / iota_a
    psi_02 = (iota_res2 - iota_b) / iota_a
    psi_03 = (iota_res3 - iota_b) / iota_a

    Br = 0
    Bvartheta = varepsilon * B0 / (1 + varepsilon * np.cos(vartheta)) * iota_profile_r
    if flag_toroildal == 0:
        Bvarphi = B0
    elif flag_toroildal == 1:    
        Bvarphi = B0 * (1 - varepsilon * np.cos(vartheta))  
    else:
        raise ValueError("Invalid value for flag_toroildal. Must be 0 or 1.")

    # Add perturbations based on flag_form
    if flag_form == 0:
        # Add the 1st perturbation
        Br += Br_perturbation(r, vartheta, varphi, m, A, iota_res, R, varepsilon, psi_0, B0, k=k, flag_type_A_fun=flag_type_A_fun, flag_toroildal=flag_toroildal)
        Bvartheta += Bvartheta_perturbation(r, vartheta, varphi, m, A, iota_res, R, varepsilon, psi_0, B0, k=k, flag_type_A_fun=flag_type_A_fun, flag_toroildal=flag_toroildal)
    elif flag_form == 1:
        # Add the 1st perturbation
        Br += Br_perturbation(r, vartheta, varphi, m1, A1, iota_res1, R, varepsilon, psi_01, B0, k=k, flag_type_A_fun=flag_type_A_fun, flag_toroildal=flag_toroildal)
        # Add the 2nd perturbation
        Br += Br_perturbation(r, vartheta, varphi, m2, A2, iota_res2, R, varepsilon, psi_02, B0, k=k, flag_type_A_fun=flag_type_A_fun, flag_toroildal=flag_toroildal)
        # Add the 3nd perturbation
        Br += Br_perturbation(r, vartheta, varphi, m3, A3, iota_res3, R, varepsilon, psi_03, B0, k=k, flag_type_A_fun=flag_type_A_fun, flag_toroildal=flag_toroildal)
        
        # Add the 1st perturbation
        Bvartheta += Bvartheta_perturbation(r, vartheta, varphi, m1, A1, iota_res1, R, varepsilon, psi_01, B0, k=k, flag_type_A_fun=flag_type_A_fun, flag_toroildal=flag_toroildal)
        # Add the 2nd perturbation
        Bvartheta += Bvartheta_perturbation(r, vartheta, varphi, m2, A2, iota_res2, R, varepsilon, psi_02, B0, k=k, flag_type_A_fun=flag_type_A_fun, flag_toroildal=flag_toroildal)
        # Add the 3nd perturbation
        Bvartheta += Bvartheta_perturbation(r, vartheta, varphi, m3, A3, iota_res3, R, varepsilon, psi_03, B0, k=k, flag_type_A_fun=flag_type_A_fun, flag_toroildal=flag_toroildal)
    elif flag_form == 2:
        # Add the 1st perturbation
        Br += Br_perturbation(r, vartheta, varphi, m1, A1, iota_res1, R, varepsilon, psi_01, B0, k=k, flag_type_A_fun=flag_type_A_fun, flag_toroildal=flag_toroildal)
        Bvartheta += Bvartheta_perturbation(r, vartheta, varphi, m1, A1, iota_res1, R, varepsilon, psi_01, B0, k=k, flag_type_A_fun=flag_type_A_fun, flag_toroildal=flag_toroildal)
        # Add the additional perturbations
        for m_i, A_i, iota_res_i in zip(m_arr, A_arr, iota_res_arr):
            psi_0i = (iota_res_i - iota_b) / iota_a
            Br += Br_perturbation(r, vartheta, varphi, m_i, A_i, iota_res_i, R, varepsilon, psi_0i, B0, k=k, flag_type_A_fun=flag_type_A_fun, flag_toroildal=flag_toroildal)
            Bvartheta += Bvartheta_perturbation(r, vartheta, varphi, m_i, A_i, iota_res_i, R, varepsilon, psi_0i, B0, k=k, flag_type_A_fun=flag_type_A_fun, flag_toroildal=flag_toroildal)
    else:
        raise ValueError("Invalid value for flag_form. Must be 0, 1 or 2.")

    return Br, Bvartheta, Bvarphi
