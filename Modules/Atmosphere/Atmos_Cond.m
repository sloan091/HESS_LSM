% =========================================================================
% Script Name   : Atmos_Cond.m
% Author        : Brandon Sloan
% Start Date    : Mar 19, 2019
% Last Updated  : Mar 19, 2019
%
% Description   : This function calculates the conductances from the canopy air
% to the measurement height using Monin Obukhov Simlarity Theory (MOST), following the
% common equations used shown by Brutsaert (1982)
%
%   INPUTS:
%   h_v    - Vegetation height [m]
%   L_g    - Obukhov length initial guess; - value indicates unstable, + indicates stable [m]
%   z      - Measurement height [m]
%   z_0m_c - Momentum roughness length for corn canopy from Moneith (2013) [m];
%   LAI    - Leaf area index [m^2 leaf/m^2 ground]
%   SAI    - Stem area index [m^2 stem/m^2 ground] 
%   U      - Streamwise velocity at measurement height [m/s]
%
%   OUTPUTS:
%   g_am - Momentum conductance from the canopy air space to measurement height [m/s]
%   g_ah - Heat conductance from the canopy air space to measurement height [m/s]
%   g_av - Vapor conductance from the canopy air space to measurement height [m/s]
% =========================================================================

function [g_am,g_ah,g_av] = Atmos_Cond(h_v,z,z_0m_c,LAI,SAI,U)
% CONSTANTS
k = 0.4;                   % von Karmen constant
R_d0 = 2/3;                % Ratio of zero plane displacement height to canopy height.
z_0m_g = 0.01;             % Momentum roughness length for bare soil [m]

% PART 1: CALCULATE THE ROUGHNESS LENGTHS, DISPLACEMENT HEIGHT AND
% STABILITY PARAMETER
%theta_a = T_a - 0.0098*z;  % This is the potential temperature at the measurment height
V = (1 - exp(-min([LAI + SAI,2])))/(1 - exp(-2));  % Fractional weight accounting for veg density [-] (Eqn 5.127)
R_z0m = z_0m_c/h_v;                                % Ratio of momentum roughness length at canopy top to canopy height.
d_0 = R_d0*h_v*V;                                  % Zero-plane displacement height [m]
z_0m = exp(V*log(h_v*R_z0m) + (1-V)*log(z_0m_g));    % Momentum roughness length for the vegetation [m] (Eqn. 5.125)
z_0h = 0.1*z_0m;                                       % Heat and vapor roughness length is assumed equal to momentum according to Oleson (2018).  DO NOT LIKE THIS        

% PART 2: CALCULATE STABILITY CORRECTED CONDUCTANCES USING MOST
% Resistance formulations from Thom 1972
r_am = (log((z-d_0)/z_0m)).^2./(U*k^2);
r_ah = (log((z-d_0)/z_0m)).*(log((z-d_0)/z_0h))./(U*k^2);

% Convert to conductance form.
g_am = 1/r_am;
g_ah = 1/r_ah;
g_av = g_ah;

end