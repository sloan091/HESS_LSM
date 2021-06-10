% =========================================================================
% Script Name   : Store_Outputs.m
% Author        : Brandon Sloan
% Start Date    : Aug 13, 2019
% Last Updated  : Aug 13, 2019
%
% Description   : This function populates the LSM_Results table with the
% model solution at each timestep.
%
%   INPUTS:
%   Timestep    - The current timestep in the simulation
%   Rad         - Outputs from the radiation module at the current timestep
%   StepResults - Results of SEB_Solver at the current timestep
%   exitflag    - Flag from the SEB_Solver at the current timestep
%
%   OUTPUTS:
%   Results   - A MATLAB table of zeros for the specified model output
%               variables
% =========================================================================

function LSM_Results = Store_Outputs(Rad,StepResults,exitflag)
%global LSM_Results
Timestep = 1;
% Energy Balance
LSM_Results.EB_sl(Timestep) = StepResults.EB_sl;
LSM_Results.EB_sh(Timestep) = StepResults.EB_sh;
LSM_Results.EB_g(Timestep) = StepResults.EB_g;
LSM_Results.Opt_Flag(Timestep) = exitflag;
	
% Sensible heat
LSM_Results.H(Timestep) = StepResults.H;
LSM_Results.H_sl(Timestep) = StepResults.H_sl;
LSM_Results.H_sh(Timestep) = StepResults.H_sh;
LSM_Results.H_g(Timestep) = StepResults.H_g;
LSM_Results.H_cont(Timestep) = StepResults.H_cont;

% Latent heat
LSM_Results.E(Timestep) = StepResults.E;
LSM_Results.E_sl(Timestep) = StepResults.E_sl;
LSM_Results.E_sh(Timestep) = StepResults.E_sh;
LSM_Results.E_g(Timestep) = StepResults.E_g;
LSM_Results.E_cont(Timestep) = StepResults.E_cont;
LSM_Results.E_sl_ww(Timestep) = StepResults.E_sl_ww;
LSM_Results.E_sh_ww(Timestep) = StepResults.E_sh_ww;

% Radiation ouputs
LSM_Results.S_r_par(Timestep) = Rad.S_r_par;
LSM_Results.S_r_nir(Timestep) = Rad.S_r_nir;
LSM_Results.L_out(Timestep) = Rad.L_out;
LSM_Results.LAI_sl(Timestep) = Rad.LAI_sl; 
LSM_Results.LAI_sh(Timestep) = Rad.LAI_sh;
LSM_Results.R_n(Timestep) = Rad.Rn_all;
	
% Stability
LSM_Results.L_obkv(Timestep) = StepResults.L_obkv;
	
% Conductances
LSM_Results.g_s_sl(Timestep) = StepResults.g_s_sl;
LSM_Results.g_s_sh(Timestep) = StepResults.g_s_sh;
LSM_Results.g_sl_c(Timestep) = StepResults.g_sl_c;
LSM_Results.g_sh_c(Timestep) = StepResults.g_sh_c;
LSM_Results.g_ah(Timestep) = StepResults.g_ah;
LSM_Results.g_av(Timestep) = StepResults.g_av;
LSM_Results.g_ah_g(Timestep) = StepResults.g_ah_g;
LSM_Results.g_av_g(Timestep) = StepResults.g_av_g;
LSM_Results.g_bvl(Timestep) = StepResults.g_bvl;
LSM_Results.g_bh(Timestep) = StepResults.g_bh;
LSM_Results.g_bv(Timestep) = StepResults.g_bv;
	
% Temperatures
LSM_Results.T_l_sl(Timestep) = StepResults.T_l_sl;
LSM_Results.T_l_sh(Timestep) = StepResults.T_l_sh;
LSM_Results.T_g(Timestep) = StepResults.T_g;
LSM_Results.T_ca(Timestep) = StepResults.T_ca;
	
% Vapor Pressures
LSM_Results.e_ca(Timestep) = StepResults.e_ca;
LSM_Results.e_s_sl(Timestep) = StepResults.e_s_sl;
LSM_Results.e_s_sh(Timestep) = StepResults.e_s_sh;
LSM_Results.e_i_sl(Timestep) = StepResults.e_i_sl;
LSM_Results.e_i_sh(Timestep) = StepResults.e_i_sh;
LSM_Results.e_g(Timestep) = StepResults.e_g;
	
% CO2 assimilation
LSM_Results.A_n_sl(Timestep) = StepResults.A_n_sl;
LSM_Results.A_n_sh(Timestep) = StepResults.A_n_sh;
LSM_Results.NEE(Timestep) = (StepResults.A_n_sl*Rad.LAI_sl + StepResults.A_n_sh*Rad.LAI_sh);
LSM_Results.ContA_sl(Timestep) = StepResults.ContA_sl;
LSM_Results.ContA_sh(Timestep) = StepResults.ContA_sh;
LSM_Results.A_sl(Timestep) = StepResults.A_sl;
LSM_Results.A_sh(Timestep) = StepResults.A_sh;
LSM_Results.GPP(Timestep) = (StepResults.A_sl*Rad.LAI_sl + StepResults.A_sh*Rad.LAI_sh);
LSM_Results.R_d_sl(Timestep) = StepResults.R_d_sl;
LSM_Results.R_d_sh(Timestep) = StepResults.R_d_sh;
LSM_Results.c_s_sl(Timestep) = StepResults.c_s_sl;
LSM_Results.c_s_sh(Timestep) = StepResults.c_s_sh;
LSM_Results.c_i_sl(Timestep) = StepResults.c_i_sl;
LSM_Results.c_i_sh(Timestep) = StepResults.c_i_sh;

% Convert to table
LSM_Results = struct2table(LSM_Results);
end