% =========================================================================
% Name   : An_ci_solver.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This function back-calculates the sunlit and shaded big leaf net CO2
% assimilation and internal leaf CO2 concentration given stomatal
% conductance values.  Specifically, we use this function for the PHM and
% beta downregulation shchemes in Sloan et al. (2021).  
%
% INPUTS
%   OS_DVs   - Decision variables of the outer solver.  Only need T_l_sl
%              and T_l_sh for this function.
%   g_sl_CO2 - CO2 conductance from in sunlit leaf to canopy airspace
%              [moles air/m^2 LA/s]
%   g_sh_CO2 - CO2 conductance from in shaded leaf to canopy airspace
%              [moles air/m^2 LA/s]
%   Const    - Physical constants
%   Rad      - Outputs from the radiative transfer model
%   BC       - Boundary conditions (i.e., environmental forcings)
%   Plant    - Plant-specific parameters
%   Flag     - Simulation settings
%
% OUTPUTS
% Note: _sl and _sh indicate sunlit and shaded big leaf, respectively.
%   c_i_sl - CO2 partial pressure inside the leaf [Pa]
%   Resid  - Leaf residual between CO2 assimilation from
%            photosynthesis and CO2 diffusive transport into the leaf
%   A_n    - Net CO2 assimilation rate [micromoles CO2/m^2 LA/s]
%   A      - CO2 assimilation rate [micromoles CO2/m^2 LA/s]
%   R_d    - Dark respiration rate [micromoles CO2/m^2 LA/s]
%   ContA  - Which process is limiting: 1 Rubisco, 2 light, and 3 production
% =========================================================================

function [c_i_sl,c_i_sh,Resid_sl,A_n_sl,A_sl,R_d_sl,ContA_sl,Resid_sh,...
    A_n_sh,A_sh,R_d_sh,ContA_sh] = An_ci_solver(OS_DVs,g_sl_CO2,g_sh_CO2,...
    Const,Rad,BC,Plant,Flag)

% PART 1: UNPACK PARAMETERS AND STATES
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Atmospheric pressure [Pa]
P_atm = BC.P_atm;

% Selects whether photosynthesis depends on temperature or not
Tflag = Flag.Tflag;

% Selects C3 or C4 photosynthesis model
Cflag = Flag.Cflag;

% Universal gas constant [J/K/kmol]
R_g = Const.R_g;

% Maximum Rubisco carboxylation rate [micromoles CO2 /m^2 LA/s]
V_max25 = Plant.V_max25;

% Maximum electron transport rate [micromoles electrons/m^2 LA/s]
J_max25 = Plant.J_max25;

% Quantum efficiency of photosystem II [micromoles electrons/ micromoles photons]
Sig_psii = Plant.Sig_psii;

% Curvature parameter for light-limited assimilation [-]
Theta_psii = Plant.Theta_psii;

% Extinction coefficient of leaf nitrogen allocation [-]
K_n = Plant.K_n;

% Sunlit and shaded leaf absorbed PAR [W/m^2]
S_sl_par = Rad.S_sl_par;
S_sh_par = Rad.S_sh_par;

% Sunlit and shaded leaf area index [m^2 LA/m^2 GA]
LAI_sl = Rad.LAI_sl;
LAI_sh = Rad.LAI_sh;

% Total leaf area index [m^2 LA/m^2 GA]
LAI = LAI_sl + LAI_sh;

% Direct beam and diffuse light extinction coefficients [-]
K_b = Rad.K_b;
K_d = Rad.K_d;

% Sunlit and shaded leaf temperature [degrees C]
T_l_sl = OS_DVs.T_l_sl;
T_l_sh = OS_DVs.T_l_sh;

% PART 2: SCALE PHOTOSYNTHETIC PARAMETERS
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Scale V_max and J_max to account for nonlinear nitrogen and light in the
% plant canopy. All outputs are per unit ground area (GA).
[V_cmax_sl,V_cmax_sh,J_cmax_sl,J_cmax_sh] = ...
    Scale_V_J(V_max25,J_max25,K_n,K_b,K_d,LAI);

% Convert absorbed PAR from units GA to LA [W/m^2 LA]
Q_sl = S_sl_par/LAI_sl;      
Q_sh = S_sh_par/LAI_sh; 

% Convert scaled V_max and J_max from units GA to LA [W/m^2 LA]
V_max_sl = V_cmax_sl/LAI_sl;
V_max_sh = V_cmax_sh/LAI_sh;
J_max_sl = J_cmax_sl/LAI_sl;
J_max_sh = J_cmax_sh/LAI_sh;


% PART 3: SOLVE FOR CO2 ASSIMILATION AND INTERNAL CO2 PRESSURE
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Select optimization and give initial interval for internal CO2
options = optimset('Display','none','TolX',1e-5);
c_int = [0,BC.c_a]; 

% Setup, run and store sunlit solver
fun_sl = @(c_i_sl) Solver_sl(c_i_sl);
c_i_sl = fzero(fun_sl,c_int,options);
[Resid_sl,A_n_sl] = Solver_sl(c_i_sl);

% Setup, run and store shaded solver
fun_sh = @(c_i_sh) Solver_sh(c_i_sh);
c_i_sh = fzero(fun_sh,c_int,options);
[Resid_sh,A_n_sh] = Solver_sh(c_i_sh);

% PART 4: SOLVER FUNCTIONS FOR THE SUNLIT AND SHADED BIG LEAVES
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

    function [Resid,A_n_sl] = Solver_sl(c_i_sl)
        
        % Select appropriate photosynthesis model based on user settings
        if isequal(Tflag,1) && isequal(Cflag,1)
            [A_n_sl,A_sl,R_d_sl,ContA_sl] = FvCB_C3wTemp(V_max_sl,...
                J_max_sl,Sig_psii,Q_sl,Theta_psii,P_atm,c_i_sl,T_l_sl,R_g);
        elseif isequal(Tflag,1) && isequal(Cflag,2)
            [A_n_sl,A_sl,R_d_sl,ContA_sl] = Collatz_C4wTemp(V_max_sl,...
                Q_sl,P_atm,c_i_sl,T_l_sl);
        elseif isequal(Tflag,0) && isequal(Cflag,2)
            [A_n_sl,A_sl,R_d_sl,ContA_sl] = Collatz_C4(V_max_sl,...
                Q_sl,P_atm,c_i_sl);
        else
            [A_n_sl,A_sl,R_d_sl,ContA_sl] = FvCB_C3(V_max_sl,...
                J_max_sl,Sig_psii,Q_sl,Theta_psii,P_atm,c_i_sl);
        end
        
        % CO2 transport into the leaf via diffusion [micromoles CO2/m^2/s]
        A_n_diff = g_sl_CO2*(BC.c_a - c_i_sl)/P_atm*10^6;
        
        % Difference between assimilated CO2 and CO2 transported into the
        % leaf
        Resid = A_n_diff - A_n_sl;
        
        
    end

    function [Resid,A_n_sh] = Solver_sh(c_i_sh)
        
        % Select appropriate photosynthesis model based on user settings
        if isequal(Tflag,1) && isequal(Cflag,1)
            [A_n_sh,A_sh,R_d_sh,ContA_sh] = FvCB_C3wTemp(V_max_sh,...
                J_max_sh,Sig_psii,Q_sh,Theta_psii,P_atm,c_i_sh,T_l_sh,R_g);
        elseif isequal(Tflag,1) && isequal(Cflag,2)
            [A_n_sh,A_sh,R_d_sh,ContA_sh] = Collatz_C4wTemp(V_max_sh,...
                Q_sh,P_atm,c_i_sh,T_l_sh);
        elseif isequal(Tflag,0) && isequal(Cflag,2)
            [A_n_sh,A_sh,R_d_sh,ContA_sh] = Collatz_C4(V_max_sh,...
                Q_sh,P_atm,c_i_sh);
        else
            [A_n_sh,A_sh,R_d_sh,ContA_sh] = FvCB_C3(V_max_sh,...
                J_max_sh,Sig_psii,Q_sh,Theta_psii,P_atm,c_i_sh);
        end
        
        % CO2 transport into the leaf via diffusion [micromoles CO2/m^2/s]
        A_n_diff = g_sh_CO2*(BC.c_a - c_i_sh)/P_atm*10^6;
        
        % Difference between assimilated CO2 and CO2 transported into the
        % leaf
        Resid = A_n_diff - A_n_sh;
           
    end

end
