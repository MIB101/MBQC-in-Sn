# MBQC-in-Sn
Mathematica and AMBiT script used in the completion of my Honours in Advanced Science majoring in Physics at UNSW in 2020.
The thesis is included as "z5061445 Matt Brand Thesis"

All code and processed data has been sorted by language/program and has the following structure and dependancies.
(Note, all file directories for importing/exporting/saving files within the code must be re-entered to match the naming conventions and file structure used on your computer.)

## AMBiT
* CIMatrix.input
  * Sample input file for AMBiT to produce the CI matrices
* Spectra.input
  * Sample input file for AMBiT to produce the E1 transitions

The leading configurations, valence configuration, number of electrons and ConfigurationAverageEnergyRange should be altered depending on the ion being considered. The CIMatrix.input files are produced in bulk by the python code. 

## Data
* CSF_Average_AMBiT
  *  Includes AMBIT output files that list the average energy and number of states in the non-relativisitc configurations up to a large energy. This is done for Sn07+ - Sn14+.
* Excel_All_MBQC
  * Excel spreadsheets that collect and detail all MBQC parameters for the various ions and energies as obtained from Mathematica/Batch_MBQC.nb

## Katana
* SpectraCI.pbs
  * Sample PBS submission files to use AMBiT to evaluate Spectra.input on Katana using normal AMBiT.
* SpectraNoCI.pbs
  * Sample PBS submission files to use AMBiT to evaluate Spectra.input on Katana using a version of AMBiT where all off-diagonal elements of the Hamiltonian are zeroed.

## Mathematica
* Packages
  * AMBiT.m
    * Provided by A/Prof Julian Berengut. Reads AMBiT CI Hamiltonians.
  * MBQCSpectra.m
    * Provides all functions to calculate the line stengths and intensity from AMBiT output.
  * ReadAMBiTOutput.m
    * Provides functions to parse various sections of the AMBiT output file.
* Batch_MBQC.nb
  * Calculates the MBQC properties for all CI matrices.
* Heuristic.nb
  * Calculates the line strength and intensity spectra based on a trial and error, empirical heuristic.
* Level_Growth_Parameter.nb
  * Calculates the level growth parameter for the various ions.
* Plasma_Spectra.nb
  * Calculates the intensity spectra of a Plasma at a given temperature and ionic composition.
* Spectra_Calculation
  * Calculates the line strength and intensity spectra using no-CI, CI and MBQC theories.
  
## Matlab
* Matlab_Input
  * This folder includes Data/Excel_All_MBQC data in matlab format
* MBQC_Energy.m
  * Plots the J^\pi averaged MBQC parameters as a function of energy for a given ion
* MBQC_JPI.m
  * Plots the MBQC parameters as a function of J^\pi  for a given energy and ion

## Python
Includes Jupyter notebooks to create input files and batch script for AMBiT to calculate CI matrices for Sn7+ - Sn14+
Bash files to batch run the CImatrix.input files are also produced.
