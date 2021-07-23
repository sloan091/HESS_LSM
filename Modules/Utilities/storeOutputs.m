% =========================================================================
% Name   : storeOutputs.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% Description   : This function populates the results table with the
% model solution at each timestep for the LSM used in Sloan et al. (2021).
% I have tried to make all definitions similar to the detailed LSM
% description in Sect. S2-S3 of Sloan et al. (2021). Note GA and LA refer
% to Ground Area and Leaf Area in the definitions below.
%
% INPUTS
%   Rad         - Outputs from the radiative transfer model
%   StepResults - Results of SEB_Solver at the current timestep
%   exitflag    - Flag from the Outer Solver at the current timestep
%
% OUTPUTS
%   LSM_Results   - A MATLAB of results for the current time step.
% =========================================================================

function LSM_Results = storeOutputs(Rad,StepResults,exitflag)


% Energy Balance Residuals [W/m^2]
%'''''''''''''''''''''''''''''''''
% Sunlit leaf
LSM_Results.EB_sl = StepResults.EB_sl;
% Shaded Leaf
LSM_Results.EB_sh = StepResults.EB_sh;
% Ground
LSM_Results.EB_g = StepResults.EB_g;
% MATLAB's nonlinear least squares solver exit flag for outer solver
LSM_Results.Opt_Flag = exitflag;
	
% Sensible heat fluxes [W/m^2]
%'''''''''''''''''''''''''''''
% Canopy airspace to atmosphere
LSM_Results.H = StepResults.H;
% Sunlit leaf to canopy airspace
LSM_Results.H_sl = StepResults.H_sl;
% Shaded leaf to canopy airspace
LSM_Results.H_sh = StepResults.H_sh;
% Ground to canopy airspace
LSM_Results.H_g = StepResults.H_g;
% Quick continuity check of H - H_sl - H_sh - H_g
LSM_Results.H_cont = StepResults.H_cont;

% Latent heat fluxes [W/m^2]
%'''''''''''''''''''''''''''''
% Canopy airspace to atmosphere
LSM_Results.E = StepResults.E;
% Sunlit leaf to canopy airspace
LSM_Results.E_sl = StepResults.E_sl;
% Shaded leaf to canopy airspace
LSM_Results.E_sh = StepResults.E_sh;
% Ground to canopy airspace
LSM_Results.E_g = StepResults.E_g;
% Quick continuity check of E - E_sl - E_sh - E_g
LSM_Results.E_cont = StepResults.E_cont;
% Well-watered sunlit leaf to canopy airspace
LSM_Results.E_sl_ww = StepResults.E_sl_ww;
% Well-watered shaded leaf to canopy airspace
LSM_Results.E_sh_ww = StepResults.E_sh_ww;

% Radiative fluxes [W/m^2]
%'''''''''''''''''''''''''
% Reflected PAR
LSM_Results.S_r_par = Rad.S_r_par;
% Reflected NIR
LSM_Results.S_r_nir = Rad.S_r_nir;
% Outgoing longwave radiation
LSM_Results.L_out = Rad.L_out;
% Overall net radiation of the soil-plant system
LSM_Results.R_n = Rad.Rn_all;

% Sunlit and Shaded Leaf Area Index [m^2 LA/m^2 GA]
%''''''''''''''''''''''''''''''''''''''''''''''''''
LSM_Results.LAI_sl = Rad.LAI_sl; 
LSM_Results.LAI_sh = Rad.LAI_sh;

	
% Obukhov Length indicating atmospheric stability [m]
%''''''''''''''''''''''''''''''''''''''''''''''''''''
LSM_Results.L_obkv = StepResults.L_obkv;
	
% Conductances [units vary]
%''''''''''''''''''''''''''
% Sunlit and shaded stomatal conductance [moles air/m^2 LA/s]
LSM_Results.g_s_sl = StepResults.g_s_sl;
LSM_Results.g_s_sh = StepResults.g_s_sh;
% Sunlit and shaded conductance from leaf to canopy air space 
% [m^3 air/m^2 GA/s or m/s]
LSM_Results.g_sl_c = StepResults.g_sl_c;
LSM_Results.g_sh_c = StepResults.g_sh_c;
% Heat and water vapor conductance from canopy air space to atmosphere
% [m^3 air/m^2 GA/s or m/s]
LSM_Results.g_ah = StepResults.g_ah;
LSM_Results.g_av = StepResults.g_av;
% Heat and water vapor conductance from ground to canopy air space
% [m^3 air/m^2 GA/s or m/s]
LSM_Results.g_ah_g = StepResults.g_ah_g;
LSM_Results.g_av_g = StepResults.g_av_g;
% Heat and water vapor leaf laminar boundary layer conductance
% [m^3 air/m^2 LA/s or m/s]
LSM_Results.g_bh = StepResults.g_bh;
LSM_Results.g_bv = StepResults.g_bv;
% Water vapor leaf laminar boundary layer conductance [moles air/m^2 LA/s]
LSM_Results.g_bvl = StepResults.g_bvl;

% Temperatures [degrees C]
%'''''''''''''''''''''''''
% Sunlit leaf
LSM_Results.T_l_sl = StepResults.T_l_sl;
% Shaded leaf
LSM_Results.T_l_sh = StepResults.T_l_sh;
% Ground
LSM_Results.T_g = StepResults.T_g;
% Canopy air space
LSM_Results.T_ca = StepResults.T_ca;
	
% Vapor Pressures [Pa]
%'''''''''''''''''''''
% Canopy air space
LSM_Results.e_ca = StepResults.e_ca;
% Surface of sunlit leaf
LSM_Results.e_s_sl = StepResults.e_s_sl;
% Surface of shaded leaf
LSM_Results.e_s_sh = StepResults.e_s_sh;
% Inside sunlit leaf
LSM_Results.e_i_sl = StepResults.e_i_sl;
% Inside shaded leaf
LSM_Results.e_i_sh = StepResults.e_i_sh;
% Ground
LSM_Results.e_g = StepResults.e_g;
	
% CO2 assimilation [units vary]
%''''''''''''''''''''''''''''''
% Sunlit and shaded net CO2 assimilation [micromoles CO2/m^2 LA/s]
LSM_Results.A_n_sl = StepResults.A_n_sl;
LSM_Results.A_n_sh = StepResults.A_n_sh;
% Sunlit and shaded CO2 assimilation [micromoles CO2/m^2 LA/s]
LSM_Results.A_sl = StepResults.A_sl;
LSM_Results.A_sh = StepResults.A_sh;
% Flag indicating Rubisco-, light-, or product-limited CO2 assimilation 
% (values of 1,2, or 3, respectively) 
LSM_Results.ContA_sl = StepResults.ContA_sl;
LSM_Results.ContA_sh = StepResults.ContA_sh;
% Net ecosystem exchange [micromoles CO2/m^2 GA/s]
LSM_Results.NEE = (StepResults.A_n_sl*Rad.LAI_sl + StepResults.A_n_sh*Rad.LAI_sh);
% Gross Primary Productivity [micromoles CO2/m^2 GA/s]
LSM_Results.GPP = (StepResults.A_sl*Rad.LAI_sl + StepResults.A_sh*Rad.LAI_sh);
% Sunlit and shaded dark respiration rate [micromoles CO2/m^2 LA/s]
LSM_Results.R_d_sl = StepResults.R_d_sl;
LSM_Results.R_d_sh = StepResults.R_d_sh;
% CO2 partial pressure at the sunlit and shaded leaf surface [Pa]
LSM_Results.c_s_sl = StepResults.c_s_sl;
LSM_Results.c_s_sh = StepResults.c_s_sh;
% CO2 partial pressure inside the sunlit and shaded leaf [Pa]
LSM_Results.c_i_sl = StepResults.c_i_sl;
LSM_Results.c_i_sh = StepResults.c_i_sh;

% Convert results to MATLAB table
LSM_Results = struct2table(LSM_Results);

end