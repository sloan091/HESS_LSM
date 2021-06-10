function [K_b,K_d] = Rad_Calc_K(x_l,LAI,zen)
% =========================================================================
% Description   : This function  calculates the diffuse radiation
% extinction coefficient by integrating the direct beam coefficient, Kb,
% over the sky hemisphere. This methodology is laid out in Bonan (2019)
%
%   INPUTS:
%   x_l - Ross parameter which shows deviation from spherical leaf
%   LAI - Leaf area index [m^2 leaf area/m^2 ground area]
%   zen - Zenith of the sun used to calculate the current Kb
%
%   OUTPUTS:
%   K_d - Diffuse radiation extinction coefficient
% =========================================================================

% Ross-Goudrian Function
dpsi = 0.01;
zenith = 0:dpsi:pi()/2;
phi1 = 0.5 - 0.633*x_l - 0.33*x_l^2;
phi2 = 0.877*(1 - 2*phi1);
G_z = phi1 + phi2*cos(zenith);
K_b = min(G_z./cos(zenith),20); % Direct beam extinction coefficient without scattering.

% Direct beam transmission coefficient for a black leaf at a specific
% cumulative LAI (Eq. 15.1)
tau_b = exp(-K_b.*LAI);
% Contribution of the sky region
contZ = 2.*sin(zenith).*cos(zenith);

% Incremental extinction coefficient for dpsi (Part of Eq. 15.5)
tau_d = tau_b*contZ'*dpsi;
K_d = -log(tau_d)/LAI;

% Calculate Kb for the current step in the model
phi1 = 0.5 - 0.633*x_l - 0.33*x_l^2;
phi2 = 0.877*(1 - 2*phi1);
G_z = phi1 + phi2*cos(zen);
K_b = min(G_z./cos(zen),20); % Direct beam extinction coefficient without scattering.

end