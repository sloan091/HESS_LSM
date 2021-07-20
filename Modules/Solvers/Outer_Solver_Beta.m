% =========================================================================
% Name   : Outer_Solver_PHM.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This function is the outer solver of the LSM with a PHM downregulation
% scheme used in Sloan et al. (2021). The function solves for the
% temperatures of the sunlit big leaf (T_l_sl), shaded big leaf (T_l_sh),
% and ground (T_g) as well as stomatal conductance of sunlit (g_s_sl) and
% shaded big leaves (g_s_sh) that balance the surface energy budget under
% conditions soil and atmospheric water stress. The details of this solver
% are covered in Sect. S2.6.2 of Sloan et al. (2021).
%
% INPUTS
%   OS_DVs_temp - Decision variables (T_l_sl, T_l_sh, T_g, g_s_sl, g_s_sh)
%                 for the LSM outer solver. The '_temp' label is because 
%                 the solver requires a vector format, whereas my functions
%                 require a structure array.
%   WW_Results  - Fluxes and states from well-watered LSM simulation
%   Const       - Physical constants
%   BC          - Boundary conditions (i.e., environmental forcings)
%   Plant       - Plant-specific parameters
%   Soil        - Soil-specific parameters
%   Flag        - Simulation settings
%
%   OUTPUTS:
%   OS_Resid    - Outer solver residuals for nonlinear least squares
%   Rad         - Outputs from the radiative transfer model
%   IS_Results  - Resulting fluxes and states from the inner solver.
% =========================================================================

function [OS_Resid,Rad,IS_Results] = Outer_Solver_Beta(OS_DVs_temp,...
    WW_Results,Const,BC,Plant,Soil,Flag)

 % PART 1: INTIALIZE DECISION VARIABLES AND RUN RADIATIVE TRANSFER MODEL
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Repackage decision variables into structure array
OS_DVs.T_l_sl = OS_DVs_temp(1);
OS_DVs.T_l_sh = OS_DVs_temp(2);
OS_DVs.T_g    = OS_DVs_temp(3);
OS_DVs.g_s_sl = OS_DVs_temp(4);
OS_DVs.g_s_sh = OS_DVs_temp(5);

% Run radiative transfer model
Rad = Rad_Transfer(BC,Plant,Soil,OS_DVs);

% PART 2: SETUP AND RUN THE INNER SOLVER
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Run inner solver 
IS_Results = Inner_Solver_PHM(OS_DVs,Const,Rad,BC,Plant,Soil,Flag);

% PART 3: STORE RESIDUALS FOR OUTER SOLVER
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Energy Balance Residuals
OS_Resid(1) = IS_Results.EB_sl;  
OS_Resid(2) = IS_Results.EB_sh;  
OS_Resid(3) = IS_Results.EB_g;

% Specify our downregulation method versus the CLM downregulation
% method. The difference is whether well-watered transpiration rate is
% fixed or not during downregulation. See third paragraph of Sect.
% S2.6.2 of Sloan et al. (2021) for discussion. For either method, the
% residual is the difference in the transpiration rate based on the outer
% solver decision variables and the well-watered transpiration rate
% mulitplied by the empirical beta factor.
if isequal(Flag.DownRegMethod,1)
    
    % Fixing well-watered stomatal conductance (and transpiration) as
    % in Sloan et al. (2021). 
    OS_Resid(4) = IS_Results.E_sl - WW_Results.E_sl*BC.beta_sl;
    OS_Resid(5) = IS_Results.E_sh - WW_Results.E_sh*BC.beta_sh;
    
else
    
    % Letting well-watered stomatal conductance change with plant
    % microclimate as in CLM v5.
    OS_Resid(4) = IS_Results.E_sl - IS_Results.E_sl_ww*BC.beta_sl;
    OS_Resid(5) = IS_Results.E_sh - IS_Results.E_sh_ww*BC.beta_sh;
    
end
    

end


