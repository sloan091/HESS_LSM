% =========================================================================
% Name   : calcPHM.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This function solves the complex Plant Hydraulic Model (PHM) in Sloan
% et al. (2021). The PHM is used to downregulate transpiration for the LSM
% used in Sloan et al. (2021). This code solves the PHM by finding the
% potentials for xylem (psi_x), sunlit (psi_l_sl) and shaded big leaves 
% (psi_l_sh) that balances water transport through each segment of the 
% soil-plant system using nonlinear least squares. See Sect. 2.2 and 
% Sect. S2.5.3 of Sloan et al. (2021) for more details.
%
% INPUTS
%   BC          - Boundary conditions (i.e., environmental forcings)
%   Plant       - Plant-specific parameters
%   Soil        - Soil-specific parameters
%   Rad         - Radiative transfer model outputs
%   Const       - Physical constants
%   IS_Results  - Resulting fluxes and states from the LSM inner solver.
%   opts        - Solver criteria for nonlinear least squares
%
% OUTPUTS
%   PHM_Results - PHM fluxes and states.
% =========================================================================

function PHM_Results = calcPHM(BC,Plant,Soil,Rad,Const,Results,opts)

% Sets PHM inputs as a persistent variable. Since the PHM is nested within
% another least squares problem (the outer solver), we chose to specify the
% initial condition of the PHM problem as the solution during the last
% iteration of the outer solver. This vastly speeds up the solution process.
persistent iter_count_PHM psi_vec

if isempty(iter_count_PHM)
    
    iter_count_PHM = 0;
    
    % ICs:        psi_l_sl [MPa] | psi_l_sh [MPa] |  psi_x [MPa]
    psi_vec = [BC.psi_s,BC.psi_s,BC.psi_s];
    
else
    iter_count_PHM = iter_count_PHM+1;
end

% Create anonymous function for solver
fx = @(x) calcPHMshell(x);

% Set solver algorithm
opts.lsqnl.Algorithm = 'levenberg-marquardt';

% Solve PHM
psi_star = lsqnonlin(fx,psi_vec,[],[],opts.lsqnl);

% Calculate and store PHM solution
[Resid,PHM_Results] = calcPHMshell(psi_star);

% Store solution as peristant variable to be used as intiial conditions for
% the next iteration of the outer solver if necessary.
psi_vec = [PHM_Results.psi_l_sl,PHM_Results.psi_l_sh,PHM_Results.psi_x];

% Hard set low ET supply solutions to 0 since many of the numerical issues
% come when dealing with near 0 transpiration.
if PHM_Results.E_sx < 1
    
    Resid(1) = 0; Resid(2) = 0; Resid(3) = 0;
    PHM_Results.E_sx = 0;
    PHM_Results.E_xl_sl = 0;
    PHM_Results.E_xl_sh = 0;
    PHM_Results.E_c_sl = 0;
    PHM_Results.E_c_sh = 0;
    PHM_Results.g_s_sl = 0;
    PHM_Results.g_s_sh = 0;
    
else
end

% Check if the solver found an acceptable closure in water transport
% between the PHM compartments.  If not, this line will stop the solver and
% throw an error for the time step.  
assert(isequal(sum(abs(Resid)<1e-1),3));

% PART 2: PHM SHELL FUNCTION USED BY SOLVER
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    function  [Resid,PHM_Results] = calcPHMshell(psi_vec)
        
        % PART 2A: DECISION VARIABLES
        %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        
        % Sunlit leaf water potential [MPa]
        psi_l_sl = psi_vec(1);
        
        % Shaded leaf water potential [MPa]
        psi_l_sh = psi_vec(2);
        
        % Xylem water potential [MPa]
        psi_x = psi_vec(3);
        
        % PART 2B: CONSTANTS
        %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        
        % Average density of water [kg/m^3]
        rho_w = Const.rho_w;
        
        % Average density of air [kg/m^3]
        rho_a = Const.rho_a;
        
        % Molecular weight of water [kg/mol]
        M_w = Const.M_w;
        
        % Molecular weight of AIR [kg/mol]
        M_a = Const.M_a;
        
        % Gravitational acceleration [m/s^2]
        g = Const.g;           
        
        % Latent heat of vaporization [J/kg]
        L_vap = Const.L_vap;
        
        % PART 2C: SOIL TO XYLEM PARAMETERS AND STATES
        %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        
        % Soil water potential near saturation [MPa]
        psi_sat = Soil.psi_sat;
        
        % Brooks-Corey soil water retention exponent [-]
        b = Soil.b;
        
        % Unsaturated hydraulic conductivity exponent
        c = 2*b+3;
        
        % Correction facotr for root extension in drying soil [-]
        d = Soil.d;
        
        % Soil water potential of current time step [MPa]
        psi_s = BC.psi_s;          
        
        % Saturated hydraulic conductivity [m/s]
        ksat = Soil.ksat/86400; 
        
        % Rooting area index [m^2 root surface/m^2 ground]
        RAI = Plant.RAI;
        
        % Fine root diameter [m]
        d_r = Plant.d_r;
        
        % Effective rooting depth [m]
        Z_r = Plant.Z_r;

        % Maximum soil to xylem conductance [m/s/MPa]
        g_sx_max = ksat/(rho_w*g)*sqrt(RAI/(d_r*Z_r))*10^6;
        
        % PART 2D: XYLEM TO LEAF PARAMETERS AND STATES
        %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        
        % Sapwood-area specific hydraulic conductivity [kg/m/s/MPA]
        k_sap = Plant.k_sap;
        
        % Sapwood area index [m^2 sapwood/m^2 ground]
        SapAI = Plant.SapAI;
        
        % Vegetation height [m]
        h_v = Plant.h_v; 
        
        % Maximum xylem to leaf conductance [m/s/MPa]
        g_xl_max = k_sap*SapAI/(rho_w*h_v); 
        
        % Xylem vulnerability curve parameter [-]
        a = Plant.a;                  
        
        % Xylem water potential at 50% loss of xylem conductivity [MPa]
        psi_x_50 = Plant.psi_x_50;
        
        % PART 2E: LEAF TO ATMOSPHERE PARAMETERS AND STATES
        %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        
        % Sunlit and shaded big leaf temperature [degrees C]
        T_l_sl = Results.T_l_sl;
        T_l_sh = Results.T_l_sh;
        
        % Atmospheric pressure [Pa]
        P_atm = BC.P_atm;  
        
        % Sunlit and shaded well-watered stomatal conductance
        % for H2O [mol air/m^2 LA/s]
        g_s_sl = Results.g_s_sl_ww;    
        g_s_sh = Results.g_s_sh_ww; 
        
        % Convert g_s to units of m/s
        g_s_slc = MolFlux2mps(g_s_sl,P_atm,T_l_sl); 
        g_s_shc = MolFlux2mps(g_s_sh,P_atm,T_l_sh);
        
        % Sunlit and shaded leaf area index [m^2 LA/m^2 GA]
        LAI_sl = Rad.LAI_sl;         
        LAI_sh = Rad.LAI_sh;
        
        % Leaf water potential at 50% stomatal closure [MPa]
        psi_l_50 = Plant.psi_l_50; 
        
        % Leaf vulnerability curve shape parameter [-]
        b_l = Plant.b_l;   
        
        % Sunlit and shaded leaf internal vapor pressure [Pa]
        e_i_sl = Results.e_i_sl;
        e_i_sh = Results.e_i_sh;
        
        % Sunlit and shaded leaf surface vapor pressure [Pa]
        e_s_sl = Results.e_s_sl; 
        e_s_sh = Results.e_s_sh;
        
        
        % PART 2F: CALCULATE FLUXES FOR EACH PHM SEGMENT
        %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        
        % Soil to xylem matric flux potential [m/s]
        Phi_sx = @(psi) (b*g_sx_max*psi)/(b-c+d)*...
            (psi_sat/psi)^((c-d)/b);
        
        % Xylem to leaf matric flux potential [m/s]
        Phi_xl = @(psi) g_xl_max*(log(exp(-psi*a) +  exp(-psi_x_50*a))/...
            a + psi);
        
        % Actual stomatal conductance (i.e., downregulated) for H2O [m/s]
        g_c = @(psi,g_s_ww,LAI) g_s_ww*LAI*2.^(-(psi/psi_l_50)^b_l);   
        
        % Soil to xylem liquid water flux [W/m^2]
        E_sx = (Phi_sx(psi_s) - Phi_sx(psi_x))*rho_w*L_vap;
        
        % The matric flux potential equations for xylem to leaf can cause
        % numerical issues, so I have added a step to deal with infinite
        % matric flux potentials.  These values should be 0, but are not
        % due to dividing by very small values.
        
        % Matric flux potential at xylem segment endpoint [m/s]
        ET_xl_x = checkInf(Phi_xl(psi_x));
        
        % Matric flux potential at sunlit and shaded leaf segment 
        % endpoints [m/s]
        ET_xl_sl = checkInf(Phi_xl(psi_l_sl));
        ET_xl_sh = checkInf(Phi_xl(psi_l_sh));

        % Xylem to sunlit leaf liquid water flux [W/m^2]
        E_xl_sl = (ET_xl_x - ET_xl_sl)*rho_w*L_vap;
        
        % Xylem to sunlit leaf liquid water flux [W/m^2]
        E_xl_sh = (ET_xl_x - ET_xl_sh)*rho_w*L_vap;
        
        % Sunlit and shaded leaf to atmosphere transpiration [W/m^2]
        E_c_sl = g_c(psi_l_sl,g_s_slc,LAI_sl)*(e_i_sl - e_s_sl)*...
            (M_w/M_a/P_atm)*rho_a*L_vap;
        E_c_sh = g_c(psi_l_sh,g_s_shc,LAI_sh)*(e_i_sh - e_s_sh)*...
            (M_w/M_a/P_atm)*rho_a*L_vap;
        
        
        % PART 2F: CALCULATE RESIDUALS AND STORE RESULTS
        %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        
        % Calculate residuals showing flow balance between PHM segments
        Resid(1) = (E_xl_sl - E_c_sl); 
        Resid(2) = (E_xl_sh - E_c_sh); 
        Resid(3) = E_xl_sl + E_xl_sh - E_sx;
        
        % Store fluxes and states
        PHM_Results.E_sx = E_sx;
        PHM_Results.E_xl_sl = E_xl_sl;
        PHM_Results.E_xl_sh = E_xl_sh;
        PHM_Results.E_c_sl = E_c_sl;
        PHM_Results.E_c_sh = E_c_sh;
        PHM_Results.psi_l_sl = psi_l_sl;
        PHM_Results.psi_l_sh = psi_l_sh;
        PHM_Results.psi_x = psi_x;
        
        % Recalculate and store downregulated stomatal conductance in molar
        % units [mol air/m^2 LA/s]
        PHM_Results.g_s_sl = g_s_sl*2.^(-(psi_l_sl/psi_l_50)^b_l); 
        PHM_Results.g_s_sh = g_s_sh*2.^(-(psi_l_sh/psi_l_50)^b_l); 
  
    end

% This small helper function is used to curb numerical issues.
    function cx = checkInf(x)
        if isinf(x)
            cx = 0;
        else
            cx = x;
        end
    end
end