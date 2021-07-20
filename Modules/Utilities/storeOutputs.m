% =========================================================================
% Name   : storeOutputs.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% Description   : This function populates the results table with the
% model solution at each timestep for the LSM used in Sloan et al. (2021).
%
% INPUTS
%   Rad         - Outputs from the radiative transfer model
%   StepResults - Results of SEB_Solver at the current timestep
%   exitflag    - Flag from the Outer Solver at the current timestep
%
% OUTPUTS
%   Results   - A MATLAB of results for the current time step.
% =========================================================================

function LSM_Results = storeOutputs(Rad,StepResults,exitflag)


% Energy Balance
LSM_Results.EB_sl = StepResults.EB_sl;
LSM_Results.EB_sh = StepResults.EB_sh;
LSM_Results.EB_g = StepResults.EB_g;
LSM_Results.Opt_Flag = exitflag;
	
% Sensible heat
LSM_Results.H = StepResults.H;
LSM_Results.H_sl = StepResults.H_sl;
LSM_Results.H_sh = StepResults.H_sh;
LSM_Results.H_g = StepResults.H_g;
LSM_Results.H_cont = StepResults.H_cont;

% Latent heat
LSM_Results.E = StepResults.E;
LSM_Results.E_sl = StepResults.E_sl;
LSM_Results.E_sh = StepResults.E_sh;
LSM_Results.E_g = StepResults.E_g;
LSM_Results.E_cont = StepResults.E_cont;
LSM_Results.E_sl_ww = StepResults.E_sl_ww;
LSM_Results.E_sh_ww = StepResults.E_sh_ww;

% Radiation ouputs
LSM_Results.S_r_par = Rad.S_r_par;
LSM_Results.S_r_nir = Rad.S_r_nir;
LSM_Results.L_out = Rad.L_out;
LSM_Results.LAI_sl = Rad.LAI_sl; 
LSM_Results.LAI_sh = Rad.LAI_sh;
LSM_Results.R_n = Rad.Rn_all;
	
% Stability
LSM_Results.L_obkv = StepResults.L_obkv;
	
% Conductances
LSM_Results.g_s_sl = StepResults.g_s_sl;
LSM_Results.g_s_sh = StepResults.g_s_sh;
LSM_Results.g_sl_c = StepResults.g_sl_c;
LSM_Results.g_sh_c = StepResults.g_sh_c;
LSM_Results.g_ah = StepResults.g_ah;
LSM_Results.g_av = StepResults.g_av;
LSM_Results.g_ah_g = StepResults.g_ah_g;
LSM_Results.g_av_g = StepResults.g_av_g;
LSM_Results.g_bvl = StepResults.g_bvl;
LSM_Results.g_bh = StepResults.g_bh;
LSM_Results.g_bv = StepResults.g_bv;
	
% Temperatures
LSM_Results.T_l_sl = StepResults.T_l_sl;
LSM_Results.T_l_sh = StepResults.T_l_sh;
LSM_Results.T_g = StepResults.T_g;
LSM_Results.T_ca = StepResults.T_ca;
	
% Vapor Pressures
LSM_Results.e_ca = StepResults.e_ca;
LSM_Results.e_s_sl = StepResults.e_s_sl;
LSM_Results.e_s_sh = StepResults.e_s_sh;
LSM_Results.e_i_sl = StepResults.e_i_sl;
LSM_Results.e_i_sh = StepResults.e_i_sh;
LSM_Results.e_g = StepResults.e_g;
	
% CO2 assimilation
LSM_Results.A_n_sl = StepResults.A_n_sl;
LSM_Results.A_n_sh = StepResults.A_n_sh;
LSM_Results.NEE = (StepResults.A_n_sl*Rad.LAI_sl + StepResults.A_n_sh*Rad.LAI_sh);
LSM_Results.ContA_sl = StepResults.ContA_sl;
LSM_Results.ContA_sh = StepResults.ContA_sh;
LSM_Results.A_sl = StepResults.A_sl;
LSM_Results.A_sh = StepResults.A_sh;
LSM_Results.GPP = (StepResults.A_sl*Rad.LAI_sl + StepResults.A_sh*Rad.LAI_sh);
LSM_Results.R_d_sl = StepResults.R_d_sl;
LSM_Results.R_d_sh = StepResults.R_d_sh;
LSM_Results.c_s_sl = StepResults.c_s_sl;
LSM_Results.c_s_sh = StepResults.c_s_sh;
LSM_Results.c_i_sl = StepResults.c_i_sl;
LSM_Results.c_i_sh = StepResults.c_i_sh;

% Convert to table
LSM_Results = struct2table(LSM_Results);

end