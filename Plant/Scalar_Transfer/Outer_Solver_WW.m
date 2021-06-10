% =========================================================================
% Name   : Outer_Solver_WW.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This function is the well-watered outer solver of the LSM used in Sloan
% et al. (2021). The function solves for the temperatures of the sunlit big
% leaf (T_l_sl), shaded big leaf (T_l_sh), and ground (T_g) that balances
% the surface energy budget under well-watered conditions (i.e., no
% transpiration downregulation). The details of this solver are covered in 
% Sect. S2.6.1 of Sloan et al. (2021).
%
% INPUTS
%   OS_DVs_temp - Decision variables (T_l_sl, T_l_sh, T_g) for the LSM
%                 outer solver. The '_temp' label is because the solver 
%                 requires a vector format, whereas my functions require
%                 a structure array.
%   Const       - Physical constants
%   BC          - Boundary conditions (i.e., environmental forcings)
%   Plant       - Plant-specific parameters
%   Soil        - Soil-specific parameters
%   opts        - Solver criteria for nonlinear least squares
%   Flag        - Simulation settings
%
% OUTPUTS
%   OS_Resid    - Outer solver residuals for nonlinear least squares
%   Rad         - Outputs from the radiative transfer model
%   IS_Results  - Resulting fluxes and states from the inner solver.
% =========================================================================

function [OS_Resid,Rad,IS_Results] = Outer_Solver_WW(OS_DVs_temp,...
    Const,BC,Plant,Soil,opts,Flag)

% PART 1: INTIALIZE DECISION VARIABLES AND RUN RADIATIVE TRANSFER MODEL
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Sets inner solver inputs as a persistent variable. Since there are two
% nested least squares problems (outer and inner solver), this allows each
% iteration of inner solver to begin with the solution of the previous
% outer solver's iteration. This vastly speeds up the solution process.
persistent iter_count IS_DVs_i

if isempty(iter_count)
    
    iter_count = 0;
    % ICs:        c_i_sl [Pa] | c_i_sh [Pa] |  e_ca [Pa]
    IS_DVs_i = [0.8*BC.c_a,   0.8*BC.c_a,   0.5*BC.e_a];
    
else
    iter_count = iter_count+1;
end

% Repackage decision variables into structure array
OS_DVs.T_l_sl = OS_DVs_temp(1);
OS_DVs.T_l_sh = OS_DVs_temp(2);
OS_DVs.T_g    = OS_DVs_temp(3);

% Run radiative transfer model
Rad = Rad_Transfer(BC,Plant,Soil,OS_DVs,Flag);

% PART 2: SETUP AND RUN THE INNER SOLVER
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Bounds for the inner solver
lb = [0,0,0];
ub = [0.99*BC.c_a,0.99*BC.c_a,5000];

% Run inner solver
IS_Results = IS_shell(IS_DVs_i,OS_DVs,Const,Rad,BC,Plant,Soil,...
    Flag,lb,ub,opts);

% PART 3: STORE RESIDUALS FOR OUTER SOLVER
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

OS_Resid(1) = IS_Results.EB_sl;  
OS_Resid(2) = IS_Results.EB_sh;  
OS_Resid(3) = IS_Results.EB_g;

% Set inner solver initial conditions for the next outer solver iteration.
IS_DVs_i    = [IS_Results.c_i_sl,IS_Results.c_i_sh,IS_Results.e_ca];

end

% PART 4: SHELL FUNCTIONS TO FACILITATE SOLVERS
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

function IS_Results = IS_shell(IS_DVs_i,OS_DVs,Const,Rad,BC,Plant,...
    Soil,Flag,lb,ub,opts)

% Run and store inner solver without atmospheric stability effects
Resid           = @(IS_DVs)Inner_Solver_WW(IS_DVs,OS_DVs,Const,Rad,BC,...
    Plant,Soil,Flag,0);
IS_DVs_opt     = lsqnonlin(Resid,IS_DVs_i,lb,ub,opts.lsqnl);
[~,IS_Results] = Inner_Solver_WW(IS_DVs_opt,OS_DVs,Const,Rad,BC,...
    Plant,Soil,Flag,1);

end