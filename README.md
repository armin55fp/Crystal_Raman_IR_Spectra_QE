# Raman and IR Spectra Generation of Crystal Structures using Quantum ESPRESSO

This repository contains Python codes and input files to generate Raman and IR spectra of crystal structures using the open-source Density Functional Theory (DFT) package, **Quantum ESPRESSO**. The scripts automate the workflow from setting DFT parameters to calculating and plotting spectra, applicable to any crystal structure (e.g., Quartz).

## Repository Contents

### Files and Scripts
1. **`parameters.py`**  
   - Contains the DFT parameters required for **SCF** (Self-Consistent Field) calculations.  
   - Includes parameters like `ibrav`, `nat`, `ecutwfc`, and supercell size.  
   - Ensure you update the molecule name (`Molecule`) and the number of atoms (`natoms`) before use.

2. **`relaxed_atomic_position.txt`**  
   - Stores the relaxed atomic positions of the crystal structure after geometry optimization.

3. **`IR_Raman_Step1_crystal.py`**  
   - Reads atomic positions from `relaxed_atomic_position.txt` and DFT parameters from `parameters.py`.  
   - Generates input files for SCF (`scf.in`) and phonon calculations (`ph.in`).  
   - Creates batch submission files (specific to your HPC cluster) and organizes all files into corresponding directories.

4. **`IR_Raman_Step2_crystal.py`**  
   - Processes outputs from phonon calculations (`ph.out`).  
   - Generates the `BORN` file required for subsequent IR and Raman calculations.

5. **`IR_Raman_Step3_crystal.py`**  
   - Aggregates outputs from all supercell DFT calculations.   
   - Calculates the IR spectra in a user-defined frequency range (`Nu_min` to `Nu_max`).  
   - Ensure you update the frequency range in the plot (`Nu_min`, `Nu_max`).

6. **`IR_Raman_Step4_crystal.py`**  
   - Gathers all outputs from supercell DFT calculations.  
   - Uses the [Phonopy-Spectroscopy](https://github.com/skelton-group/Phonopy-Spectroscopy) package to generate the final IR and Raman spectra plots and results.  
   - For more details about the `Phonopy-Spectroscopy` package and its functionality, refer to the [official GitHub repository](https://github.com/skelton-group/Phonopy-Spectroscopy).

7. **`qe2outcar.py`**  
   - Converts Quantum ESPRESSO outputs to **VASP-compatible** `OUTCAR` format.  
   - Utilized in `IR_Raman_Step3_crystal.py` and `IR_Raman_Step4_crystal.py`.

8. **`poscar2qe.py`**  
   - Converts VASP `POSCAR` files to Quantum ESPRESSO-compatible input formats.  
   - Used within `IR_Raman_Step3_crystal.py` and `IR_Raman_Step4_crystal.py`.

## How to Use

1. **Prepare Input Files**  
   - Update `parameters.py` with the required SCF and phonon parameters.  
   - Provide the optimized atomic structure in `relaxed_atomic_position.txt`.

2. **Run Step-by-Step Scripts**  
   - Execute `IR_Raman_Step1_crystal.py` to generate SCF and phonon input files and batch submission scripts.  
   - Run Quantum ESPRESSO calculations for SCF and phonon properties on your HPC cluster.  
   - Use `IR_Raman_Step2_crystal.py` to generate the `BORN` file.  
   - Execute `IR_Raman_Step3_crystal.py` and `IR_Raman_Step4_crystal.py` for spectra calculations and plotting.  

3. **Visualization**  
   - The `Phonopy-Spectroscopy` package is used in `IR_Raman_Step4_crystal.py` to create final plots of IR and Raman spectra.  
   - The resulting spectra are saved in the `Result` folder.

## Requirements

- **Software**:  
  - Quantum ESPRESSO  
  - [Phonopy-Spectroscopy](https://github.com/skelton-group/Phonopy-Spectroscopy)  
  - Python 3.x with the following libraries:  
    - `numpy`  
    - `matplotlib`  
    - `os`  

- **Cluster Environment**:  
  - Ensure batch submission scripts are configured for your specific cluster environment.

## Key Notes

- Update molecule-specific parameters like the name (`Molecule`), number of atoms (`natoms`), and frequency range (`Nu_min`, `Nu_max`) in relevant scripts.  
- Ensure the VASP and Quantum ESPRESSO conversion tools (`qe2outcar.py` and `poscar2qe.py`) are properly integrated into your workflow if needed.


## Contact

For any questions or support, please feel free to open an issue or contact the repository maintainer.
