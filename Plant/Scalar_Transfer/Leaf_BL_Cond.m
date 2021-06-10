% =========================================================================
% Name   : Leaf_BL_cond.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% Description   : This function calculates the heat and vapor conductance
% for the leaf boundary layer using simplified empirical functions provided
% in Oleson et al. (2018).  Note the original equation from Goudriaan (1977)
% and laid out in Montieth and Unsworth (2013) already accounts for
% transport from both sides of the leaf.  Here I have added a factor of 0.5
% to make the conductance value account for one side of the leaf. My
% sensible and latent heat formulation explicitly account for one or two
% sides of the leaf. Overall this value is empirical and has minor
% influence on the results.
%
% INPUTS
%   d_l    - Characteristic leaf length [m]
%   u_star - Friction velocity [m/s]
%
% OUTPUTS
%   g_bh - Heat conductance for laminar leaf boundary layer [m/s]
%   g_bh - Vapor conductance for laminar leaf boundary layer [m/s]
%
% REFERENCES
%   (1) Oleson, K. W. et al. (2018). Technical Description of the version 5
%   of the Community Land Model (CLM). 
%
%   (2) Goudriaan, J. (1977). Crop micrometeorology: a simulation study.
%   Wageningen.
%
%   (3) Monteith, J., & Unsworth, M. (2013). Principles of Environmental
%   Physics: Plants, Animals, and the Atmosphere: Fourth Edition.
%   https://doi.org/10.1016/C2010-0-66393-0
% =========================================================================

function [g_bh,g_bv] = Leaf_BL_Cond(u_star,d_l)

% Turbulent transfer coefficient for velocity on leaves [m/s] - (1)
C_s = 0.01;                

% Heat conductance for laminar leaf boundary layer [m/s] - (1)
g_bh = 0.5*C_s*(u_star/d_l);  

%Vapor conductance for laminar leaf boundary layer [m/s] - (1)
g_bv = g_bh;

end