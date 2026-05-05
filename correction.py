r"""
correction.py
===========
Python translation of nonpolar.f90.

Implements the three nonrigid electrostatic energy contributions for the
Delphi nonpolar solvation workflow:

    CALC_E1            - vacuum Coulomb sum  E1
    CALC_E2            - in-solvent Coulomb sum  E2

Mathematical background
-----------------------
The pairwise Coulomb sums are:

    E = \Sum_{i < j}  q_i · q_j / r_ij          (guard: r_ij > 1e-15)

The nonrigid energy components are:

	\Delta E_{12} = (E2 - E1) · ec2 / epsp
"""

import numpy as np


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

def _pairwise_coulomb_sum(positions: np.ndarray, charges: np.ndarray) -> float:
    r"""
    Compute the upper-triangle Coulomb sum: \Sum_{i < j} q_i · q_j / r_ij.

    Parameters
    ----------
    positions : np.ndarray, shape (3, N)
        Atomic coordinates [Å].
    charges : np.ndarray, shape (N,)
        Atomic partial charges [e].

    Returns
    -------
    float
        Coulomb sum (un-scaled; multiply by ec2/eps outside).
    """
    n = positions.shape[1]
    if n < 2:
        return 0.0

    # Upper-triangle index pairs
    i_idx, j_idx = np.triu_indices(n, k=1)

    # Pairwise displacement vectors and distances
    diff = positions[:, i_idx] - positions[:, j_idx]    # shape (3, N_pairs)
    r    = np.sqrt(np.sum(diff ** 2, axis=0))           # shape (N_pairs,)

    # Guard:
    mask = r > 1.0e-15

    energy = np.sum(charges[i_idx[mask]] * charges[j_idx[mask]] / r[mask])
    return float(energy)


# ---------------------------------------------------------------------------
# Compute E1 and E2
# ---------------------------------------------------------------------------

def CALC_E1(vatmpos: np.ndarray, atmchr: np.ndarray) -> float:
    """
    Compute the vacuum electrostatic energy E1.

    Parameters
    ----------
    vatmpos : np.ndarray, shape (3, N)
        Vacuum atomic positions [Å]
        atmchr : np.ndarray, shape (N,)
        Atomic partial charges [e]

    Returns
    -------
    float
        E1 — raw Coulomb sum in vacuum (unscaled by ec2/epsp).
    """
    E1 = _pairwise_coulomb_sum(vatmpos, atmchr)
    print(f" CALCULATED E1 = {E1}")
    return E1


def CALC_E2(atmpos: np.ndarray, atmchr: np.ndarray) -> float:
    """
    Compute the in-solvent Coulomb sum E2.

    Parameters
    ----------
    atmpos : np.ndarray, shape (3, N)
        Solvated atomic positions [Å]
    atmchr : np.ndarray, shape (N,)
        Atomic partial charges [e]

    Returns
    -------
    float
        E2 — raw Coulomb sum in solvent geometry (unscaled).
    """
    E2 = _pairwise_coulomb_sum(atmpos, atmchr)
    print(f"CALCULATED E2 = {E2}")
    return E2


def nonrigid_energy(
    atmpos:       np.ndarray,
    vatmpos:      np.ndarray,
    atmchr:       np.ndarray,
    ec2:          float,
    epsp:         float,
) -> float:
    """
    Assemble the total nonrigid solvation energy and write a report file.

    Parameters
    ----------
    atmpos : np.ndarray, shape (3, N)
        Solvated atomic positions [Å].
    vatmpos : np.ndarray, shape (3, N)
        Vacuum atomic positions [Å].
    atmchr : np.ndarray, shape (N,)
        Atomic partial charges [e].
    ec2 : float
        Conversion constant (e.g. 332.06364 kcal·Å/mol/e²).
    epsp : float
        Protein interior dielectric constant.
    Returns
    -------
    float
        Total nonrigid energy [kcal/mol].
    """

    ec2 = 332.06364
    epsp = 1.0
    # --- Step 1: compute the raw Coulomb sums ---
    E1     = CALC_E1(vatmpos, atmchr)          # vacuum
    E2     = CALC_E2(atmpos, atmchr)           # solvent geometry

    # --- Step 2: energy correction ---
    delta_E12 = (E2 - E1) * ec2 / epsp

    return delta_E12

if __name__ == "__main__":
    import sys
    import os
    import re
    from readin import readin1, readin1_vac

    if len(sys.argv) < 2:
        print("Usage: python correction.py <pdb_id>")
        sys.exit(1)

    pdb_id = sys.argv[1]

    epsp = 1.0
    ec2 = 332.06364

    try:
        atm_pos, atm_chr, atm_rad, grid_params = readin1(pdb_id)
    except FileNotFoundError:
        atm_pos, atm_chr, atm_rad, grid_params = readin1(pdb_id, base_path='./')

    try:
        vatmpos = readin1_vac(pdb_id)
    except FileNotFoundError:
        vatmpos = readin1_vac(pdb_id, base_path='./')

    delta_E12 = nonrigid_energy(atm_pos, vatmpos, atm_chr, ec2, epsp)

    out_file = f"{pdb_id}_delta_E12.txt"
    with open(out_file, "w") as f:
        f.write(str(delta_E12))
    
    print(f"Saved delta_E12 = {delta_E12} to {out_file}")
