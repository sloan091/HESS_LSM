% =========================================================================
% Name   : Inner_Solver_Beta.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This function is the inner solver of the LSM with a Beta downregulation
% scheme used in Sloan et al. (2021). Unlike the well-watered LSM, this
% inner solver is not actually a nested nonlinear least squares problem
% within the outer solver.  Rather, this function simply calculates the
% fluxes and states of carbon, water and heat transport for the dual
% source, two-big leaf LSM using the outer solver's decision variables 
% (T_l_sl, T_l_sh, T_g, g_s_sl, g_s_sh). These fluxes are used to check the
% energy balance in the outer solver and that the transpiration matches the
% values given by the well-watered tranpiration rate corrected by the
% empirical beta function. For reference, this code covers from the fourth
% box on in Fig. S3 of Sloan et al. (2021).
%
% INPUTS
%   OS_DVs - Decision variables (T_l_sl, T_l_sh, T_g, g_s_sl, g_s_sh)
%                 for the LSM outer solver.
%   Const  - Physical constants
%   Rad    - Outputs from the radiative transfer model
%   BC     - Boundary conditions (i.e., environmental forcings)
%   Plant  - Plant-specific parameters
%   Soil   - Soil-specific parameters
%   Flag   - Simulation settings
%
% OUTPUTS
%   IS_Results - Resulting fluxes and states from the inner solver.
% =========================================================================

function IS_Results = Inner_Solver_Beta(OS_DVs,Const,Rad,BC,...
    Plant,Soil,Flag)

% PART 1: UPDATE CONDUCTANCES AND STATES BASED ON TEMPERATURE
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Total LAI needed for some calculations
LAI_total = Rad.LAI_sl + Rad.LAI_sh;

% Calculate laminar leaf boundary layer conductance [m/s]
[IS_Results.g_bh,IS_Results.g_bv] = Leaf_BL_Cond(BC.u_star,Plant.d_l);
% Convert from m/s to mol/m^2/s.
IS_Results.g_bvl       = MpsFlux2Mol(IS_Results.g_bv,BC.P_atm,0.5*(OS_DVs.T_l_sl...
    + OS_DVs.T_l_sh)); 

% Calculate soil conductance [m/s]
[IS_Results.g_ah_g,IS_Results.g_av_g] = Soil_Cond(LAI_total,Plant.SAI,BC.u_star,...
    Soil.theta_i,Soil.theta_sat,BC.theta_s,Soil.psi_sat,Soil.b,OS_DVs.T_g);
% Correction factor for soil moisture
IS_Results.alpha_g = exp(BC.psi_s*Const.g/(Const.R_v*10^3*OS_DVs.T_g)); 
% Saturated vapor pressure at the soil temperature [Pa]
IS_Results.e_sat_g = CC(OS_DVs.T_g);        
% Vapor pressure at the ground surface corrected for soil moisture [Pa]
IS_Results.e_g = IS_Results.alpha_g*IS_Results.e_sat_g;                                 

% Calculate atmospheric conductance via Monin-Obukhov Similarity Theory
% [m/s]
[~,IS_Results.g_ah,IS_Results.g_av] = Atmos_Cond(Plant.h_v,...
    BC.z,Plant.z_0m_c,LAI_total,Plant.SAI,BC.U);

% Calculate internal leaf vapor pressure using Clausius-Clapeyron [Pa]
IS_Results.e_i_sl = CC(OS_DVs.T_l_sl);
IS_Results.e_i_sh = CC(OS_DVs.T_l_sh);


% PART 2: RE-CALCULATE CO2 ASSIMILATION AND CONCENTRATION UNDER
% DOWNREGULATION
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% CO2 conductance from inside leaf to canopy airspace
g_sl_CO2 = 1/(1.4/IS_Results.g_bvl + 1.6/(OS_DVs.g_s_sl));
g_sh_CO2 = 1/(1.4/IS_Results.g_bvl + 1.6/(OS_DVs.g_s_sh));

% Iteratively solve for CO2 assimilation and concentration given you know
% stomatal conductances
[c_i_sl,c_i_sh,~,A_n_sl,A_sl,R_d_sl,ContA_sl,~,A_n_sh,A_sh,R_d_sh,...
    ContA_sh] = An_ci_solver(OS_DVs,g_sl_CO2,g_sh_CO2,Const,...
    Rad,BC,Plant,Flag);

% Back calculate surface CO2 concentration [Pa]
c_s_sl = BC.c_a - (BC.P_atm*A_n_sl*1.4)/(IS_Results.g_bvl*Plant.ns*10^6);
c_s_sh = BC.c_a - (BC.P_atm*A_n_sh*1.4)/(IS_Results.g_bvl*Plant.ns*10^6);

% PART 3: UPDATE REMAINING STATES ALTERED BY TRANSPIRATION DOWNREGULATION
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Calculate the new leaf ot canopy airspace conductances [m/s]
g_sl_beta = Plant.ns*MolFlux2mps(Rad.LAI_sl*(IS_Results.g_bvl...
    *OS_DVs.g_s_sl)/(IS_Results.g_bvl +...
    OS_DVs.g_s_sl),BC.P_atm,OS_DVs.T_l_sl); 
g_sh_beta = Plant.ns*MolFlux2mps(Rad.LAI_sh*(IS_Results.g_bvl...
    *OS_DVs.g_s_sh)/(IS_Results.g_bvl +...
    OS_DVs.g_s_sh),BC.P_atm,OS_DVs.T_l_sh); 

% Calculate canopy airspace temperature [degrees C]
T_ca = (IS_Results.g_ah*BC.T_a + 2*Rad.LAI_sl*IS_Results.g_bh*OS_DVs.T_l_sl +...
    2*Rad.LAI_sh*IS_Results.g_bh*OS_DVs.T_l_sh + ...
    IS_Results.g_ah_g*OS_DVs.T_g)/(IS_Results.g_ah + 2*Rad.LAI_sl*IS_Results.g_bh...
    + 2*Rad.LAI_sh*IS_Results.g_bh + IS_Results.g_ah_g);

% Calculate canopy airspace vapor pressure [Pa]
e_ca = (IS_Results.g_av*BC.e_a + g_sl_beta*IS_Results.e_i_sl + ...
    g_sh_beta*IS_Results.e_i_sh + IS_Results.g_av_g*IS_Results.e_g)...
    /(IS_Results.g_av + g_sl_beta + g_sh_beta + IS_Results.g_av_g);

% Calculate leaf surface vapor pressure [Pa]
e_s_sl = (OS_DVs.g_s_sl*IS_Results.e_i_sl ...
    + IS_Results.g_bvl*e_ca)/(OS_DVs.g_s_sl...
    + IS_Results.g_bvl);
e_s_sh = (OS_DVs.g_s_sh*IS_Results.e_i_sh...
    + IS_Results.g_bvl*e_ca)/(OS_DVs.g_s_sh...
    + IS_Results.g_bvl);

% Re-calculate Medlyn g_s [moles H2O/m^2/s}. NOTE: Only used if the CLM
% downregulation method is used. Not used in Sloan et al. (2021)
IS_Results.g_s_sl_ww = Medlyn_gs(Plant.g_o,Plant.g_1,IS_Results.e_i_sl - e_s_sl,...
    A_n_sl,BC.c_a,BC.P_atm);
IS_Results.g_s_sh_ww = Medlyn_gs(Plant.g_o,Plant.g_1,IS_Results.e_i_sh - e_s_sh,...
    A_n_sh,BC.c_a,BC.P_atm);

% PART 4: CALCULATE FLUXES AND STORE WITH STATES
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% PART 4a: CALCULATE SENSIBLE HEAT FLUXES (H)
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Canopy airspace to atmosphere [W/m^2]
IS_Results.H = Const.rho_a*Const.c_p*IS_Results.g_ah*(T_ca - BC.T_a);
% Sunlit leaf to canopy airspace [W/m^2]
IS_Results.H_sl = 2*Rad.LAI_sl*Const.rho_a*Const.c_p*IS_Results.g_bh*...
    (OS_DVs.T_l_sl - T_ca);
% Shaded leaf to canopy airspace [W/m^2]
IS_Results.H_sh = 2*Rad.LAI_sh*Const.rho_a*Const.c_p*IS_Results.g_bh*...
    (OS_DVs.T_l_sh - T_ca);
% Ground to canopy airspace [W/m^2]
IS_Results.H_g  = Const.rho_a*Const.c_p*IS_Results.g_ah_g*...
    (OS_DVs.T_g - T_ca);
% continuity check
IS_Results.H_cont = IS_Results.H - (IS_Results.H_sl + ...
    IS_Results.H_sh + IS_Results.H_g);

% PART 4b: CALCULATE LATENT HEAT FLUXES (E)
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Canopy airspace to atmosphere [W/m^2]
IS_Results.E = Const.rho_a*IS_Results.g_av*(Const.M_w/Const.M_a/BC.P_atm)...
    *(e_ca - BC.e_a)*Const.L_vap;
% Sunlit leaf to canopy airspace [W/m^2]
IS_Results.E_sl = Const.rho_a*g_sl_beta*(Const.M_w/Const.M_a/BC.P_atm)*...
    (IS_Results.e_i_sl - e_ca)*Const.L_vap;
% Shaded leaf to canopy airspace [W/m^2]
IS_Results.E_sh = Const.rho_a*g_sh_beta*(Const.M_w/Const.M_a/BC.P_atm)*...
    (IS_Results.e_i_sh - e_ca)*Const.L_vap;
% Ground to canopy airspace [W/m^2]
IS_Results.E_g  = Const.rho_a*IS_Results.g_av_g*(Const.M_w/Const.M_a/...
    BC.P_atm)*(IS_Results.e_g - e_ca)*Const.L_vap;
% Continuity check
IS_Results.E_cont = IS_Results.E - (IS_Results.E_sl + ...
    IS_Results.E_sh + IS_Results.E_g);

% Calculate the new ETww values. NOTE: Only used if the CLM
% downregulation method is used. Not used in Sloan et al. (2021)
g_sl_ww = Plant.ns*MolFlux2mps(Rad.LAI_sl*(IS_Results.g_bvl...
    *IS_Results.g_s_sl_ww)/(IS_Results.g_bvl +...
    IS_Results.g_s_sl_ww),BC.P_atm,OS_DVs.T_l_sl); 
g_sh_ww = Plant.ns*MolFlux2mps(Rad.LAI_sh*(IS_Results.g_bvl...
    *IS_Results.g_s_sh_ww)/(IS_Results.g_bvl +...
    IS_Results.g_s_sh_ww),BC.P_atm,OS_DVs.T_l_sh); 
IS_Results.E_sl_ww = Const.rho_a*g_sl_ww*(Const.M_w/Const.M_a/BC.P_atm)...
    *(IS_Results.e_i_sl - e_ca)*Const.L_vap;
IS_Results.E_sh_ww = Const.rho_a*g_sh_ww*(Const.M_w/Const.M_a/BC.P_atm)...
    *(IS_Results.e_i_sh - e_ca)*Const.L_vap;

% Calculate Obukhov length just for reference to see how atmospheric
% instabilities may affect results
IS_Results.L_obkv = (-(BC.u_star^3)*Const.rho_a)/(Const.k*Const.g*...
    ((IS_Results.H/((BC.T_a + 273.15)*Const.c_p)) + ...
    0.61*(IS_Results.E/Const.L_vap)));

% PART 4c: CALCULATE GROUND AND LEAF SURFACE ENERGY BUDGET RESIDUALS
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

IS_Results.EB_sl = Rad.Rn_sl - IS_Results.E_sl - IS_Results.H_sl;
IS_Results.EB_sh = Rad.Rn_sh - IS_Results.E_sh - IS_Results.H_sh;
IS_Results.EB_g  = Rad.Rn_g  - IS_Results.E_g - IS_Results.H_g - BC.G;

% PART 4d: STORE REMAINING CONDUCTANCES
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

IS_Results.g_s_sl = OS_DVs.g_s_sl;
IS_Results.g_s_sh = OS_DVs.g_s_sh;
IS_Results.g_sl_c = g_sl_beta;
IS_Results.g_sh_c = g_sh_beta;

% PART 4E: STORE STATES AND CO2 ASSIMILATION
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Temperature
IS_Results.T_ca   = T_ca;
IS_Results.T_l_sl = OS_DVs.T_l_sl;
IS_Results.T_l_sh = OS_DVs.T_l_sh;
IS_Results.T_g    = OS_DVs.T_g;
    
% Vapor Pressure
IS_Results.e_ca   = e_ca;
IS_Results.e_s_sl = e_s_sl;
IS_Results.e_s_sh = e_s_sh;

% Assimilation
IS_Results.A_n_sl   = A_n_sl;
IS_Results.A_n_sh   = A_n_sh;
IS_Results.ContA_sl = ContA_sl;
IS_Results.ContA_sh = ContA_sh;
IS_Results.A_sl     = A_sl;
IS_Results.A_sh     = A_sh;
IS_Results.R_d_sl   = R_d_sl;
IS_Results.R_d_sh   = R_d_sh;
IS_Results.c_s_sl   = c_s_sl;
IS_Results.c_s_sh   = c_s_sh;
IS_Results.c_i_sl   = c_i_sl;
IS_Results.c_i_sh   = c_i_sh;

end
