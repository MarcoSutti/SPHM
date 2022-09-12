# SPHM: A MATLAB package for SPH simulations
 
This is a collection of MATLAB files for SPH simulations.

For documentation, please refer to the document:

- [SPHM: a MATLAB package for Smoothed Particle Hydrodynamics simulations](https://www.marcosutti.net/files/SPHM_a_MATLAB_package_for_SPH_simulations.pdf),
M. Sutti, Tech. report, September 2022.

The code is mainly based on the corresponding FORTRAN code contained in:

- Liu, G. R. and Liu, M. B. Smoothed particle hydrodynamics: a meshfree particle method. World Scientific, Singapore, 2003.

If there are any problems or bugs, feel free to email me at msutti (at) ncts.tw


## I) Version History

- Ver 2, 12 September 2022.
- Ver 1, 11 July 2021: initial release.


## II) Contents

- eos: contains the equations of state for ideal gas and artificial water.
- examples: contains parameters and data for the 1D shock tube and the
            2D shear cavity problems.
- plots: contains the plots generated from the simulation data.
- postprocessing: contains utilities for plotting and visualizing data.
- results: contains the matfiles generated by *Driver_SPHM*.
- utilities: contains all the core functions for the SPH simulations.
- videos: contains videos of the simulations.


## III) Installation and Usage

No installation is required.
- Use the *Driver_SPHM* to perform an SPH simulation.
The results of the simulations will be saved into matfiles
in the results folder. 
- The script *Driver_Animation* can be used to generate a video of 
the simulation. The script *Driver_shocktube_profiles* generates the 
profile plots for the shock tube example.
- The script *Driver_shear_cavity_steady_state* generates the plots for
the shear-driven cavity problem at the steady state.
- The script *Driver_Sod_shocktube* can be used to visualize the profiles
of the exact solution to the classical Sod shock tube problem.


## IV) License

Code written by me is GPL licensed.
