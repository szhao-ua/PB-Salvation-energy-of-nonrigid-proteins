# Protein Solvation Energy Calculations under Nonrigid Conditions

This repository contains scripts, configuration files, and documentation for calculating polar solvation energy of nonrigid proteins in the Poisson-Boltzmann theory. The primary goal of this repository is to compute the correction term in the calculation of the polar solvation energy.

## Requirements

To run the tools and scripts in this repository, ensure you have the following installed:
- **Python 3.6+**
- **NumPy** (`pip install numpy`)
<!-- - **GROMACS v5.0.5** (or compatible version for structure preparation)
- **DelPhi v8.5.0** (for solving the Poisson-Boltzmann Equation via traditional method)
- **APBS** (Adaptive Poisson-Boltzmann Solver) -->

## Repository Structure

```text
.
├── APBS.in               # APBS template configuration file
├── DelPhi.prm            # DelPhi template configuration file
├── correction.py         # Main script to compute nonrigid energy contributions (ΔE_12)
├── readin.py             # Utility script to parse atomic coordinates (.xyzqr)
├── 70set.txt             # List of the 70 PDB IDs analyzed
├── data-TIP3P/           # Directory for solvent-equilibrated atomic data (.xyzqr files)
├── data-Xtal/            # Directory for vacuum-equilibrated atomic data (.xyzqr files)
└── readme.md             # This documentation file
```

---

## Usage

### 1. Structure Preparation

All protein structures analyzed in this study were sourced from the Protein Data Bank (PDB) [1]. Structure preparation and molecular dynamics (MD) minimizations were conducted in **GROMACS v5.0.5** [2] using the **AMBER99SB** force field [3], with all titratable residues maintained in their standard charged states.

To evaluate the methods under nonrigid conditions, two distinct conformational states were utilized for each protein:

- **Solvent-equilibrated structure (Explicitly solvated state)**: The proteins were solvated using TIP3P water molecules, with ions added as necessary for neutralization. Further details on this structure preparation process can be found in Chakravorty et al. [4].
- **Vacuum-equilibrated structure (Reference state)**: Crystals were protonated and subjected to energy minimization with heavy atomic restraints (using a force constant of $1 \times 10^6\,\text{kJ}\,\text{mol}^{-1}\,\text{nm}^{-2}$) to preserve the native backbone topology. Further details on this structure preparation process can be found in Chakravorty et al. [4].

### 2. Calculating Solvation Energy

#### Using DelPhi

To run tests on the 70-protein set using DelPhi, we use the following input template. Detailed descriptions of these parameters can be found in the [DelPhi User Manual (v8.5.0)](https://compbio.clemson.edu/lab/media/download/DelPhiUserManual_v8_5_0.pdf).

```text
scale=2.0
perfil=70.0
in(pdb,file="1TG0.pdb")
in(siz,file="1TG0.siz")
in(crg,file="1TG0.crg")
indi=1.0
exdi=80.0
prbrad=1.4
salt=0.15
dencut=0.759
bndcon=2
maxc=0.0001
linit=1000
nonit=800
gaussian=0
energy(s)
out(energy,file="1TG0_TIP3P_TRAD_1.dat")
```

When calculating the energy, DelPhi provides the **Corrected reaction field energy**, which is the polar solvation energy in units of kT. This is converted to kcal/mol using the conversion factor:
$$1\,\text{kT} = 0.5922\,\text{kcal/mol}$$

#### Using APBS

Similarly, to run tests using APBS, we use the following input setup. Detailed descriptions of these parameters can be found in the [APBS User Manual](https://apbs.readthedocs.io/en/latest/using/examples/solvation-energies.html#polar-solvation).

```text
# Protein Solvation Energy Calculation

read
mol pqr 1TG0.pqr
end

# --- SOLVATED STATE (in water) ---
elec name solvated
mg-manual                # Specify the mode for APBS to run
dime 97 97 97            # The grid dimensions
grid 0.33 0.33 0.33      # Grid spacing
gcent mol 1              # Center the grid on molecule 1
mol 1                    # Perform the calculation on molecule 1
npbe                     # Solve the nonlinear Poisson-Boltzmann equation
bcfl sdh                 # single Debye-Hückel boundary conditions
pdie 1.0                 # Solute dielectric
sdie 80                  # Solvent dielectric
chgm spl0                # Spline-based discretization of the delta functions
srfm mol                 # Molecular surface definition
srad 1.4                 # Solvent probe radius (for molecular surface)
swin 0.3                 # Solvent surface spline window (not used here)
sdens 10.0               # Sphere density for accessibility object
temp 298.15              # Temperature
calcenergy total         # Calculate energies
calcforce no             # Do not calculate forces
end

# --- REFERENCE STATE (in vacuum, sdie=1) ---
elec name reference
mg-manual
dime 97 97 97
grid 0.33 0.33 0.33
gcent mol 1
mol 1
npbe
bcfl sdh
pdie 1.0
sdie 1.0
chgm spl0
srfm mol
srad 1.4
swin 0.3
sdens 10.0
temp 298.15
calcenergy total
calcforce no
end

# Print solvation energy = solvated - reference
print elecEnergy solvated - reference end

quit
```

Running APBS using this input file provides the **Global net ELEC energy** value in units of kJ/mol. This is converted to kcal/mol using the conversion factor: $$1\,\text{kcal/mol} = 4.184\,\text{kJ/mol}$$

### 3. Correcting for Nonrigid Contributions ($\Delta E_{12}$)

To evaluate conditions under nonrigid topologies, we must correct the rigid solvation energy by adding the $\Delta E_{12}$ nonrigid energy contribution. 

#### Data Preparation
The atomic coordinate data (`.xyzqr` files) must be placed in the respective directories:
- **Solvated data**: Place inside `data-TIP3P/`
- **Vacuum data**: Place inside `data-Xtal/`

*Note: The script `readin.py` is responsible for parsing this data. If different directory names are used, update the `base_path` variable inside `readin.py`.*

#### Running the Correction Script
The `correction.py` script computes the raw vacuum energy ($E_1$), the raw in-solvent energy ($E_2$), and the final nonrigid contribution ($\Delta E_{12}$). 

To run the calculation, use the following bash command:

```bash
python correction.py <PDB_ID>
```
*(e.g., `python correction.py 1TG0`)*

The output ($\Delta E_{12}$) is automatically saved to a text file in the working directory (e.g., `1TG0_delta_E12.txt`). 

#### Final Adjusted Solvation Energy
For both DelPhi and APBS, the final corrected polar solvation energy is simply the sum of the rigid polar solvation energy (converted to kcal/mol) and the nonrigid energy contribution ($\Delta E_{12}$):
$$\Delta E = E_{\text{rigid}} (\text{in kcal/mol}) + \Delta E_{12}$$

---
## References
[1] Berman HM, et al. *The Protein Data Bank*. Nucleic Acids Res. 2000.  
[2] Van Der Spoel D, et al. *GROMACS: fast, flexible, and free*. J Comput Chem. 2005.  
[3] Ponder JW, Case DA. *Force fields for protein simulations*. Adv Protein Chem. 2003.  
[4] Chakravorty A, et al. *Reproducing the ensemble average polar solvation energy of a protein from a single structure: Gaussian-based smooth dielectric function for macromolecular modeling*. J Chem Theory Comput. 2018.
