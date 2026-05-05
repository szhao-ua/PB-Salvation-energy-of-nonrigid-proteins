import numpy as np
import os

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def readin1(fname, fradius=None, dcel=None, base_path=None):
    """
    Reads atom data from an .xyzqr file and defines the computational grid.
    
    Args:
        fname (str): The filename (without extension).
        fradius (float): The padding radius (edge).
        dcel (float): Grid spacing.
        base_path (str): Optional directory containing the data files. If None, it will parse usrdata.in
        
    Returns:
        tuple: (atm_pos, atm_chr, atm_rad, grid_params)
    """

    base_path = 'data-TIP3P/'

    # 1. Construct File Path
    full_path = os.path.join(base_path, fname + ".xyzqr")
    print(f"Reading file: {full_path}")
    
    if not os.path.exists(full_path):
        raise FileNotFoundError(f"File not found: {full_path}")

    # 2. Read Data
    # Assumes columns are: X, Y, Z, Charge, Radius (whitespace separated)
    try:
        data = np.loadtxt(full_path)
    except ValueError:
        raise ValueError("File format error. Expected 5 numeric columns (X Y Z Q R).")

    # Extract columns
    atm_pos = data[:, 0:3].T   # Columns 1-3 -> X, Y, Z
    atm_chr = data[:, 3]       # Column 4 -> Charge
    atm_rad = data[:, 4]       # Column 5 -> Radius
    
    natm = atm_pos.shape[1]
    print(f"number of atoms = {natm}")

    grid_params = None

    return atm_pos, atm_chr, atm_rad, grid_params

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def readin1_vac(fname, base_path=None):
    """
    Reads explicit vacuum atomic coordinates from an .xyzqr file.
    
    Args:
        fname (str): The filename (without extension).
        base_path (str): Optional directory containing the vacuum data files. If None, it will parse usrdata.in
        
    Returns:
        np.array: Vacuum positions (vatmpos) with shape (3, n_atm).
    """

    base_path = 'data-Xtal/'

    # 1. Construct File Path
    full_path = os.path.join(base_path, fname + ".xyzqr")
    print(f"Reading vacuum file: {full_path}")
    
    if not os.path.exists(full_path):
        raise FileNotFoundError(f"Vacuum file not found: {full_path}")

    # 2. Read Data
    # Assumes columns are: X, Y, Z, Charge, Radius
    try:
        data = np.loadtxt(full_path)
    except ValueError:
        raise ValueError("File format error. Expected numeric columns.")

    # 3. Extract Coordinates
    vatmpos = data[:, 0:3].T
    
    natm = vatmpos.shape[1]
    print(f"number of atoms (vac) = {natm}")
    
    return vatmpos

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++