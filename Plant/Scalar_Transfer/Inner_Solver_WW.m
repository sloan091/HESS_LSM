% =========================================================================
% Name   : Inner_Solver_WW.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This function is the well-watered inner solver of the LSM used in Sloan
% et al. (2021). The function solves for the internal carbon concentration
% of the sunlit (c_i_sl) and shaded (c_i_sh) big leaves as well as canopy
% airspace vapor pressure (e_ca) that balance leaf to atmosphere carbon,
% vapor and heat transport under well-watered conditions (i.e., no
% transpiration downregulation). The details of this solver are covered in 
% Sect. S2.6.1 of Sloan et al. (2021).
%
% INPUTS
%   IS_DVs_temp - Decision variables [c_i_sl, c_i_sh, e_ca] for the LSM
%                  inner solver. The '_temp' label is because the solver 
%                   requires a vector format, whereas my functions require
%                   a structure array.
%   OS_DVs      - Decision variables (T_l_sl, T_l_sh, T_g) for the LSM
%                 outer solver at the current iteration.
%   Const       - Physical constants
%   Rad         - Outputs from the radiative transfer model
%   BC          - Boundary conditions (i.e., environmental forcings)
%   Plant       - Plant-specific parameters
%   Soil        - Soil-specific parameters
%   opts        - Solver criteria for nonlinear least squares
%   Flag        - Simulation settings
%   RFlag       - Results flag determines whether to store inner solver (1)
%                 results or skip storage (0) to speed up solution.
%
% OUTPUTS
%   IS_Resid    - Inner solver residuals for nonlinear least squares
%   IS_Results  - Resulting fluxes and states from the inner solver.
% =========================================================================

function [IS_Resid,IS_Results] = Inner_Solver_WW(IS_DVs_temp,OS_DVs,...
    Const,Rad,BC,Plant,Soil,Flag,Rflag)


% PART 1: UNPACK DECISION VARIABLES
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Internal leaf CO2 concentration for sunlit and shaded leaves [Pa] and
% canopy water vapor pressure [Pa].
IS_DVs.c_i_sl = IS_DVs_temp(1); 
IS_DVs.c_i_sh = IS_DVs_temp(2); 
IS_DVs.e_ca   = IS_DVs_temp(3);

% Total LAI needed for some calculations
LAI_total      = Rad.LAI_sl + Rad.LAI_sh;

% PART 2: CALCULATE LEAF SCALE CO2 ASSIMILATION
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Scale photosynthesis parameters to account for nonlinear nitrogen and
% light profiles
[V_cmax_sl,V_cmax_sh,J_cmax_sl,J_cmax_sh] = Scale_V_J(Plant.V_max25,...
    Plant.J_max25,Plant.K_n,Rad.K_b,Rad.K_d,LAI_total);

% Convert parameters from unit ground area (GA) to unit leaf area (LA)

% Absorbed PAR [W/m^2 LA]
Q_sl     = Rad.S_sl_par/Rad.LAI_sl;
Q_sh     = Rad.S_sh_par/Rad.LAI_sh;
% Max Rubisco-limited carboxylation rate [micromoles CO2/m^2 LA/s]
V_max_sl = V_cmax_sl/Rad.LAI_sl;     
V_max_sh = V_cmax_sh/Rad.LAI_sh;
% Max RuBP-limited (light) carboxylation rate [micromoles CO2/m^2 LA/s]
J_max_sl = J_cmax_sl/Rad.LAI_sl;     
J_max_sh = J_cmax_sh/Rad.LAI_sh;   
  

% Calculate net CO2 assimilation (A_n) using the selected form of the
% Farquhar, von Caemmerer, and Berry photosynthetis model for sunlit and
% shaded big leaves
if isequal(Flag.Tflag,1) && isequal(Flag.Cflag,1)
    % C3 photosynthesis with temperature dependence
    [A_n_sl,A_sl,R_d_sl,ContA_sl] = Farquhar_C3(V_max_sl,J_max_sl,...
        Plant.Sig_psii,Q_sl,Plant.Theta_psii,BC.P_atm,IS_DVs.c_i_sl,1,...
        OS_DVs.T_l_sl,Const.R_g);
    [A_n_sh,A_sh,R_d_sh,ContA_sh] = Farquhar_C3(V_max_sh,J_max_sh,...
        Plant.Sig_psii,Q_sh,Plant.Theta_psii,BC.P_atm,IS_DVs.c_i_sh,...
        1,OS_DVs.T_l_sh,Const.R_g);
elseif isequal(Flag.Tflag,1) && isequal(Flag.Cflag,2)
    % C4 photosynthesis with temperature dependence
    [A_n_sl,A_sl,R_d_sl,ContA_sl] = Farquhar_C4(V_max_sl,Q_sl,BC.P_atm,...
        IS_DVs.c_i_sl,1,OS_DVs.T_l_sl);
    [A_n_sh,A_sh,R_d_sh,ContA_sh] = Farquhar_C4(V_max_sh,Q_sh,BC.P_atm,...
        IS_DVs.c_i_sh,1,OS_DVs.T_l_sh);
elseif isequal(Flag.Tflag,0) && isequal(Flag.Cflag,2)
    % C4 photosynthesis without temperature dependence
    [A_n_sl,A_sl,R_d_sl,ContA_sl] = Farquhar_C4_NoTemp(V_max_sl,Q_sl,...
        BC.P_atm,IS_DVs.c_i_sl,1);
    [A_n_sh,A_sh,R_d_sh,ContA_sh] = Farquhar_C4_NoTemp(V_max_sh,Q_sh,...
        BC.P_atm,IS_DVs.c_i_sh,1);
else
    % C3 photosynthesis without temperature dependence
    [A_n_sl,A_sl,R_d_sl,ContA_sl] = Farquhar_C3_NoTemp(V_max_sl,...
        J_max_sl,Plant.Sig_psii,Q_sl,Plant.Theta_psii,BC.P_atm,...
        IS_DVs.c_i_sl,1);
    [A_n_sh,A_sh,R_d_sh,ContA_sh] = Farquhar_C3_NoTemp(V_max_sh,...
        J_max_sh,Plant.Sig_psii,Q_sh,Plant.Theta_psii,BC.P_atm,...
        IS_DVs.c_i_sh,1);
end

% PART 3: CALCULATE LEAF STOMATAL CONDUCTANCE AND STATES FROM DIFFUSION
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Calculate laminar leaf boundary layer conductance [m/s]
[g_bh,g_bv] = Leaf_BL_Cond(BC.u_star,Plant.d_l);
% Convert from m/s to mol/m^2/s.
g_bvl       = MpsFlux2Mol(g_bv,BC.P_atm,0.5*(OS_DVs.T_l_sl...
    + OS_DVs.T_l_sh)); 

% Back-calculate leaf surface CO2 pressure for Medlyn Eqn. [Pa]
c_s_sl = BC.c_a - (BC.P_atm*A_n_sl*1.4)/(g_bvl*Plant.ns*10^6); 
c_s_sh = BC.c_a - (BC.P_atm*A_n_sh*1.4)/(g_bvl*Plant.ns*10^6); 

% Calculate internal leaf vapor pressure using Clausius-Clapeyron [Pa]
e_i_sl = CC(OS_DVs.T_l_sl);
e_i_sh = CC(OS_DVs.T_l_sh);

% Define some constants to perform quadratic estimation of g_s
C_1_sl = 1.6*A_n_sl/((c_s_sl/BC.P_atm)*10^6); 
C_1_sh = 1.6*A_n_sh/((c_s_sh/BC.P_atm)*10^6);
C_2_sl = (e_i_sl - IS_DVs.e_ca)/1000; 
C_2_sh = (e_i_sh - IS_DVs.e_ca)/1000;

% Take larger root of quadratic equation to solver for g_s as done in CLM
% v5. This trick substitutes the diffusion equations for vapor and CO2 and
% H2O into Medlyn's stomatal response equation for an simpler solution.
g_s_sl = max(roots([1,-(2*Plant.g_o + 2*C_1_sl + (C_1_sl^2*Plant.g_1^2/...
    (g_bvl*C_2_sl))),Plant.g_o^2 + 2*C_1_sl*Plant.g_o + C_1_sl^2*...
    (1 - Plant.g_1^2/C_2_sl)]));
g_s_sh = max(roots([1,-(2*Plant.g_o + 2*C_1_sh + (C_1_sh^2*Plant.g_1^2/...
    (g_bvl*C_2_sh))),Plant.g_o^2 + 2*C_1_sh*Plant.g_o + C_1_sh^2*...
    (1 - Plant.g_1^2/C_2_sh)]));

% Handle numerical issues due to small g_s values
if ~isreal(g_s_sl) || g_s_sl<0
    g_s_sl = 1e-8;
else
end

if ~isreal(g_s_sh) || g_s_sh<0
    g_s_sh = 1e-8;
else
end

% Convert from moles to m/s
g_s_slc = MolFlux2mps(g_s_sl,BC.P_atm,OS_DVs.T_l_sl);
g_s_shc = MolFlux2mps(g_s_sh,BC.P_atm,OS_DVs.T_l_sh);

% Calculate leaf surface vapor pressure [Pa]
e_s_sl = (g_s_sl*e_i_sl + g_bvl*IS_DVs.e_ca)/(g_s_sl + g_bvl);
e_s_sh = (g_s_sh*e_i_sh + g_bvl*IS_DVs.e_ca)/(g_s_sh + g_bvl);

% Back-calculate internal leaf CO2
IS_DVs.c_i_slch = BC.c_a - (1.4/g_bvl + 1.6/g_s_sl)*BC.P_atm*A_n_sl/10^6;
IS_DVs.c_i_shch = BC.c_a - (1.4/g_bvl + 1.6/g_s_sh)*BC.P_atm*A_n_sh/10^6;

% First two residuals for the inner solver
IS_Resid(1) = IS_DVs.c_i_sl - IS_DVs.c_i_slch;
IS_Resid(2) = IS_DVs.c_i_sh - IS_DVs.c_i_shch;

% PART 4: CALCULATE CANOPY AIR SPACE VAPOR PRESSURE
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Calculate soil conductance [m/s]
[g_ah_g,g_av_g] = Soil_Cond(LAI_total,Plant.SAI,BC.u_star,Soil.theta_i...
    ,Soil.theta_sat,BC.theta_s,Soil.psi_sat,Soil.b,OS_DVs.T_g);
% Correction factor for soil moisture
alpha_g = exp(BC.psi_s*Const.g/(Const.R_v*10^3*OS_DVs.T_g)); 
% Saturated vapor pressure at the soil temperature [Pa]
e_sat_g = CC(OS_DVs.T_g);        
% Vapor pressure at the ground surface corrected for soil moisture [Pa]
e_g     = alpha_g*e_sat_g;                                 

% Calculate atmospheric conductance via Monin-Obukhov Similarity Theory
% [m/s]
[~,g_ah,g_av] = Atmos_Cond(Plant.h_v,BC.z,Plant.z_0m_c,LAI_total,...
    Plant.SAI,BC.U);

% Calculate canopy airspace temperature [degrees C]
T_ca = (g_ah*BC.T_a + 2*Rad.LAI_sl*g_bh*OS_DVs.T_l_sl + ...
    2*Rad.LAI_sh*g_bh*OS_DVs.T_l_sh + g_ah_g*OS_DVs.T_g)...
    /(g_ah + 2*Rad.LAI_sl*g_bh + 2*Rad.LAI_sh*g_bh + g_ah_g);

% Collect calculate leaf to canopy airspace conductance [m/s]
g_sl = Plant.ns*Rad.LAI_sl*(g_bv*g_s_slc)/(g_bv + g_s_slc);
g_sh = Plant.ns*Rad.LAI_sh*(g_bv*g_s_shc)/(g_bv + g_s_shc);

% Calculate canopy airspace vapor pressure [Pa]
IS_DVs.e_ca_ch = (g_av*BC.e_a + g_sl*e_i_sl + g_sh*e_i_sh + ...
    g_av_g*e_g)/(g_av + g_sl + g_sh + g_av_g);

% Final residual for the inner solver
IS_Resid(3) = IS_DVs.e_ca - IS_DVs.e_ca_ch;

% PART 4: CALCULATE FLUXES AND STORE WITH STATES
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

if isequal(Rflag,0)
    % If initial solver run do not store results to speed up solution
    IS_Results = [];
else
    
    % PART 4a: CALCULATE SENSIBLE HEAT FLUXES (H)
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    % Canopy airspace to atmosphere [W/m^2]
    IS_Results.H    = Const.rho_a*Const.c_p*g_ah*(T_ca - BC.T_a);
    % Sunlit leaf to canopy airspace [W/m^2]
    IS_Results.H_sl = 2*Rad.LAI_sl*Const.rho_a*Const.c_p*g_bh...
        *(OS_DVs.T_l_sl - T_ca);
    % Shaded leaf to canopy airspace [W/m^2]
    IS_Results.H_sh = 2*Rad.LAI_sh*Const.rho_a*Const.c_p*g_bh...
        *(OS_DVs.T_l_sh - T_ca);
    % Ground to canopy airspace [W/m^2]
    IS_Results.H_g  = Const.rho_a*Const.c_p*g_ah_g*(OS_DVs.T_g - T_ca);
    % Continuity check
    IS_Results.H_cont = IS_Results.H - (IS_Results.H_sl +...
        IS_Results.H_sh + IS_Results.H_g);
    
    % PART 4b: CALCULATE LATENT HEAT FLUXES (E)
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    % Canopy airspace to atmosphere [W/m^2]
    IS_Results.E    = Const.rho_a*g_av*(Const.M_w/Const.M_a/BC.P_atm)...
        *(IS_DVs.e_ca - BC.e_a)*Const.L_vap;
    % Sunlit leaf to canopy airspace [W/m^2]
    IS_Results.E_sl = Const.rho_a*g_sl*(Const.M_w/Const.M_a/BC.P_atm)...
        *(e_i_sl - IS_DVs.e_ca)*Const.L_vap;
    % Shaded leaf to canopy airspace [W/m^2]
    IS_Results.E_sh = Const.rho_a*g_sh*(Const.M_w/Const.M_a/BC.P_atm)...
        *(e_i_sh - IS_DVs.e_ca)*Const.L_vap;
    % Ground to canopy airspace [W/m^2]
    IS_Results.E_g  = Const.rho_a*g_av_g*(Const.M_w/Const.M_a/BC.P_atm)...
        *(e_g - IS_DVs.e_ca)*Const.L_vap;
    % Quick continuity check
    IS_Results.E_cont = IS_Results.E - (IS_Results.E_sl +...
        IS_Results.E_sh + IS_Results.E_g);
    % Store well-watered latent heat fluxes separately in case
    % transpiration downregulation is enabled.
    IS_Results.E_sl_ww = IS_Results.E_sl;
    IS_Results.E_sh_ww = IS_Results.E_sh;
    
    % Calculate Obukhov length just for reference to see how atmospheric
    % instabilities may affect results
    IS_Results.L_obkv = (-(BC.u_star^3)*Const.rho_a)/(Const.k*...
        Const.g*((IS_Results.H/((BC.T_a + 273.15)*Const.c_p)) +...
        0.61*(IS_Results.E/Const.L_vap)));
    
    % PART 4c: CALCULATE GROUND AND LEAF SURFACE ENERGY BUDGET RESIDUALS
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    % Sunlit big leaf [W/m^2]
    IS_Results.EB_sl = Rad.Rn_sl - IS_Results.E_sl - IS_Results.H_sl;
    % Shaded big leaf [W/m^2]
    IS_Results.EB_sh = Rad.Rn_sh - IS_Results.E_sh - IS_Results.H_sh;
    % Ground [W/m^2]
    IS_Results.EB_g  = Rad.Rn_g  - IS_Results.E_g - IS_Results.H_g - BC.G;
   
    % PART 4d: STORE CONDUCTANCES
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    IS_Results.g_s_sl = g_s_sl;
    IS_Results.g_s_sh = g_s_sh;
    IS_Results.g_sl_c = g_sl;
    IS_Results.g_sh_c = g_sh;
    IS_Results.g_ah   = g_ah;
    IS_Results.g_av   = g_av;
    IS_Results.g_ah_g = g_ah_g;
    IS_Results.g_av_g = g_av_g;
    IS_Results.g_bvl  = g_bvl;
    IS_Results.g_bh   = g_bh;
    IS_Results.g_bv   = g_bv;
    
    % PART 4E: STORE STATES AND CO2 ASSIMILATION
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    % Temperature
    IS_Results.T_ca   = T_ca;
    IS_Results.T_l_sl = OS_DVs.T_l_sl;
    IS_Results.T_l_sh = OS_DVs.T_l_sh;
    IS_Results.T_g    = OS_DVs.T_g;
    
    % Vapor Pressure
    IS_Results.e_ca   = IS_DVs.e_ca;
    IS_Results.e_s_sl = e_s_sl;
    IS_Results.e_s_sh = e_s_sh;
    IS_Results.e_i_sl = e_i_sl;
    IS_Results.e_i_sh = e_i_sh;
    IS_Results.e_g    = e_g;
    
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
    IS_Results.c_i_sl   = IS_DVs.c_i_sl;
    IS_Results.c_i_sh   = IS_DVs.c_i_sh;
    
end

end
