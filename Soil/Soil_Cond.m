% =========================================================================
% Name   : Soil_Cond.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% Description   : This function calculates the heat and vapor conductance
% between the soil and canopy airspace as laid out in Oleson et al. (2018).
% See Sect. S2.3.1 of Sloan et al. (2021) for full details.
%
% INPUTS
%   LAI       - Leaf area index [m^2 LA/m^2 GA]
%   SAI       - Stem area index [m^2 SA/m^2 GA]
%   u_star    - Friction velocity [m/s]
%   theta_i   - Water content at which soil restricts conductance 
%               [m^3 water/m^3 soil]
%   theta_sat - Saturated soil water content or porosity
%               [m^3 water/m^3 soil]
%   theta_s   - Soil water content [m^3 water/m^3 soil]
%   psi_sat   - Saturated soil water content [MPa]
%   b         - Brooks-Corey soil water retention exponent [-]
%   T_g       - Soil temperature [degress C]
%
%   OUTPUTS:
%   g_ah_g - Soil to canopy air space heat conductance [m/s]
%   g_ah_g - Soil to canopy air space vapor conductance  [m/s]
% =========================================================================

function [g_ah_g,g_av_g] = Soil_Cond(LAI,SAI,u_star,theta_i,theta_sat,...
    theta_s,psi_sat,b,T_g)

% CONSTANTS
k = 0.4;                   % von Karmen constant
z_0m_g = 0.01;             % Momentum roughness length for bare soil [m]
C_s_dense = 0.004;         % Turbulent transfer coefficient for dense canopy (Eqn 5.120)
D_max = 0.015;             % Max thickness of the dry layer [m]
psi_air = (-10^4)*9.8*10^(-3);   % Potential of air dry soil [MPa]
nu = 1.5*10^-5;            % Kinematic viscosity of water [m^2/s]
g = 9.8;                   % gravitational acceleration (m/s^2)

% PART 1: CALCULATE HEAT TRANSFER CONDUCTANCE FROM SOIL TO CANOPY AIR
W = exp(-LAI - SAI);                            % Weighting factor (Eqn 5.118)
C_s_bare  = k/0.13*(z_0m_g*u_star/nu)^(-0.45);  % Turbulent transfer coefficient for bare soil (Eqn 5.121)
C_s = C_s_bare*W + C_s_dense*(1 - W);           % Weighted turbulent transfer coefficient (Eqn 5.118)
g_ah_g = C_s * u_star;                          % 5.118 from Oleson (2018) [m/s]

% PART 2: CALCULATE VAPOR TRANSFER CONDUCTANCE FROM SOIL PORES TO SURFACE
% AND ADD IN SERIES TO CONDUCTANCE FROM SOIL SURFACE TO CANOPY AIR

theta_air = psi2theta(theta_sat,psi_air,psi_sat,b);      % Air dry water content [-] (Eqn 5.78)

% Dry surface layer thickness [m] (Eqn 5.77)
if theta_s < theta_i
    DSL = D_max*(theta_i - theta_s)/(theta_i - theta_air);
    Phi_air = theta_sat - theta_air;                       % Air filled pore space [-] (Eqn 5.80)
    tau = Phi_air^2*(Phi_air/theta_sat)^(3/b);             % Soil tortuosity [-] (Eqn 5.79)
    D_v = (2.12*10^(-5))*((T_g + 273.15)/273.15)^1.75;     % Vapor Diffusivity as a function of T [m^2/s] (Eqn 5.81)
    g_soil = D_v*tau/DSL;                                  % The conductance from soil (Eqn 5.76)
    g_av_g = (g_soil*g_ah_g)/(g_soil + g_ah_g);            % Total vapor conductance from soil [m/s]
else
    g_av_g = g_ah_g;    % No additional resistance due to soil porosity
end

end