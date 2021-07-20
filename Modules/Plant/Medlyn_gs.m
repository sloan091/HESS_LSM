% =========================================================================
% Name   : Medlyn_gs.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This function calculates stomatal conductance (g_s) using the optimality
% model derived by Medlyn et al. (2011).  The Medlyn model assumes plants
% minimize water lost for a certain carbon gain and is derived using the
% calculus of variations. This formulation only applies for well-watered
% condtions and is used in conjunction with a transpiration downregulation
% scheme to represent soil water stress.
%
% INPUTS
%   g_o   - Minimum stomatal conductance for H2O [moles air/m^2/s]
%   g_1   - Empirical Medlyn parameter [kPa^0.5]
%   VPD   - Vapor pressure difference between inside the leaf and at 
%           the leaf surface [Pa]
%   A_n   - Estimate of net CO2 assimilation [micromoles CO2/m^2/s]
%   c_s   - CO2 partial pressure at the leaf surface [Pa]
%   P_atm - Atmospheric pressure [Pa]   
%
% OUTPUTS
%   g_s - Stomatal conductance of H20 [moles air/m^2/s]
%
% REFERENCES
%   (1) Medlyn, B. E. et al. (2011). Reconciling the optimal and empirical 
%   approaches to modelling stomatal conductance. Global Change Biology, 
%   17(6), 2134–2144. https://doi.org/10.1111/j.1365-2486.2010.02375.x
% =========================================================================

function g_s = Medlyn_gs(g_o,g_1,VPD,A_n,c_s,P_atm)

% Converts from Pa to moles CO2/moles air with a factor to convert A_n to
% units of moles CO2/m^2/s
c_s_conv = c_s/P_atm*10^6;

% Stomatal conductance of H20 [moles air/m^2/s]
g_s = g_o + 1.6*(1 + g_1/sqrt(VPD/1000))*A_n/c_s_conv;  

end
        





