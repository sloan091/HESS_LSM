# HESS_LSM
This repository contains the MATLAB codes for the dual-source, two-big-leaf land surface model used Sloan et al. (2021). The provided codes allow users to recreate the simulation results for the US-Me2 ponderosa pine site using the five separate transpiration downregulation schemes.

Sloan, B. P., Thompson, S. E., and Feng, X.: Plant Hydraulic Transport Controls Transpiration Response to Soil Water Stress, Hydrol. Earth Syst. Sci., https://doi.org/10.5194/hess-2020-671, 2021.

# Instructions
The LSM can be run either serially (LSM_Serial.m) or in parallel (LSM_Parallel.m) mode using MATLAB.  For the simulation in Sloan et al. (2021), the serial version takes approximately 30 minutes to run, whereas the parallel mode takes approximately 3 minutes using 24 cores. To run the LSM, simply run either the LSM_Serial.m or LSM_Parallel.m script in MATLAB and ensure all the folders are included in your MATLAB search path.  I have included three separate parameter files in the Parameter folder that can be used to create all the simulation in the paper. Additionally, I have included the LSM simulation results used in the paper in the Results folder. Each results .mat file contains 3 MATLAB tables of the selected states and fluxes from the LSM run assuming: 1) the well-watered conditions (_WW), 2) beta transpiration downregulation (_Beta), and 3) PHM transpiration downregulation (_PHM).  If the flags for either transpiration downregulation scheme are turned off (see below), the tables contain zeros. Please see the Utilities/storeOutputs.m for descriptions and units of the LSM results.    

There are three separate parameter files: Parameters_HESS_PHM_and_beta_s.mat, Parameters_HESS_PHM_and_beta_2L.mat, and Parameters_HESS_PHM_and_beta_dyn.mat. Each runs the same well-watered and PHM downregulation simulation, but uses a different beta downregulation scheme: 1) a single beta applied to both sunlit and shaded big leaves (beta_s), 2) a beta applied separately to sunlit and shaded big leaves (beta_2L), and 3) the 'dynamic' beta that varies with atmospheric moisture demand applied to separately to sunlit and shaded big leaves (beta_dyn). The user can update these parameter files to run only certain simulations by editing the Flag structure.  By default, any LSM simulation will calculate the well-watered simulation as it is used for the downregulation schemes.  Setting the Betaflag or PHMflag values can toggle these downregulation simulations on or off.  Please see the Parameter_File_Creator_HESS.m for the appropriate flag values.  

The other folders contain all the shells, solvers and sub-modules used by the LSM and sorted by with an appropriate descriptor of which portion of the soil-plant-atmosphere continuum they belong.  Each script and function contains a brief description and reference to the Supplement sections of Sloan et al. (2021) as well as relevant literature used to create it.  Furtheremore, I have tried to keep consistent naming and descriptions of model inputs and outputs so the dependence between functions are easily followed.  Additionally, I have included significant annotation in these codes to make units clear and have also attempted to name the variables similarly to those included in the Supplement of Sloan et al. (2021).  For those interested in diggin through the code, I would recommend starting with Solvers/runLSMParallel.m or Solvers/runLSMSerial.m as these are the shell functions that cooridnate all the submodules.  Please refer to Sect. S2-S3 Sloan et al. (2021) for the detailed model description.
