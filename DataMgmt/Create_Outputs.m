% =========================================================================
% Script Name   : Create_Outputs.m
% Author        : Brandon Sloan
% Start Date    : Aug 12, 2019
% Last Updated  : Aug 12, 2019
%
% Description   : This function creates the MATLAB table that will be used
% to store the model outputs.
%
%   INPUTS:
%   Timesteps - The number of timesteps to be run in the simulation
%
%   OUTPUTS:
%   Results   - A MATLAB table of zeros for the specified model output
%               variables
% =========================================================================

function Results = Create_Outputs(Timesteps)

steps = zeros(Timesteps,1);

% Energy balance
EB_sl = steps;
EB_sh = steps;
EB_g = steps;
Opt_Flag = steps;

% Sensible heat
H = steps;
H_sl = steps;
H_sh = steps;
H_g = steps;
H_cont = steps;

% Latent heat
E = steps;
E_sl = steps;
E_sh = steps;
E_g = steps;
E_cont = steps;
E_sl_ww = steps;
E_sh_ww = steps;

% Radiation Ouputs
S_r_par = steps;
S_r_nir = steps;
L_out = steps;
LAI_sl = steps; 
LAI_sh = steps;
R_n = steps;

% Stability
L_obukhov = steps;

% Conductances
g_s_sl = steps;
g_s_sh = steps;
g_sl_c = steps;
g_sh_c = steps;
g_ah = steps;
g_av = steps;
g_ah_g = steps;
g_av_g = steps;
g_bvl = steps;
g_bh = steps;
g_bv = steps;

% Temperature
T_l_sl = steps;
T_l_sh = steps;
T_g = steps;
T_ca = steps;

% Vapor Pressure
e_ca = steps;
e_s_sl = steps;
e_s_sh = steps;
e_i_sl = steps;
e_i_sh = steps;
e_g = steps;

% Assimilation
A_n_sl = steps;
A_n_sh = steps;
NEE = steps;
ContA_sl = steps;
ContA_sh = steps;
A_sl = steps;
A_sh = steps;
GPP = steps;
R_d_sl = steps;
R_d_sh = steps;
c_s_sl = steps;
c_s_sh = steps;
c_i_sl = steps;
c_i_sh = steps;

Results = table(EB_sl,EB_sh,EB_g,Opt_Flag,H,H_sl,H_sh,H_g,H_cont,E,E_sl,E_sh,E_g,E_cont,E_sl_ww,E_sh_ww,...
    S_r_par,S_r_nir,L_out,LAI_sl,LAI_sh,R_n,L_obukhov,g_s_sl,g_s_sh,g_sl_c,g_sh_c,g_ah,g_av,...
    g_ah_g,g_av_g,g_bvl,g_bh,g_bv,T_l_sl,T_l_sh,T_g,T_ca,e_ca,e_s_sl,e_s_sh,e_i_sl,e_i_sh,e_g,...
    A_n_sl,A_n_sh,NEE,ContA_sl,ContA_sh,A_sl,A_sh,GPP,R_d_sl,R_d_sh,c_s_sl,c_s_sh,c_i_sl,c_i_sh);
end