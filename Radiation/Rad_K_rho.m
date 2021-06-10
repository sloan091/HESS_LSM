function [K_b,K_bp_par,K_bp_nir,K_d,K_dp_par,K_dp_nir,rho_cb_par,rho_cd_par,rho_cb_nir,rho_cd_nir] = Rad_K_rho(Plant,Soil,BC,Flag)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1: GRAB CORRECT LAI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isequal(Flag.LAIflag,1)
    LAI = Plant.LAI(BC.Step);
else
    LAI = Plant.LAI;
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 2: CALCULATE DIFFUSE EXTINCTION COEFFICIENT, K_d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ross-Goudrian Function
dpsi = 0.01;
zen_range = 0:dpsi:pi()/2;
phi1 = 0.5 - 0.633*Plant.x_lad - 0.33*Plant.x_lad^2;
phi2 = 0.877*(1 - 2*phi1);
G_z = phi1 + phi2*cos(zen_range);
K_b = min(G_z./cos(zen_range),20); % Direct beam extinction coefficient without scattering.

% Direct beam transmission coefficient for a black leaf at a specific cumulative LAI
tau_b = exp(-K_b.*LAI);
% Contribution of the sky region
contZ = 2.*sin(zen_range).*cos(zen_range);

% Incremental extinction coefficient for dpsi (Part of Eq. 15.5)
tau_d = tau_b*contZ'*dpsi;
K_d = -log(tau_d)/LAI;
K_dp_par = sqrt(Plant.a_l_par)*K_d; % Diffuse extinction coefficient adjusted for scattering (Kd' in DPF(1997))
K_dp_nir = sqrt(Plant.a_l_nir)*K_d; % Diffuse extinction coefficient adjusted for scattering (Kd' in DPF(1997))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 3: CALCULATE DIFFUSE CANOPY REFLECTANCE FOR PAR AND NIR WITHOUT SOIL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate canopy reflectances without accounting for soil
rho_h_par = (1 - sqrt(Plant.a_l_par))/(1 + sqrt(Plant.a_l_par));
rho_b_par = 2*K_b./(K_b + K_d)*rho_h_par;
rho_d_par = rho_b_par*contZ'*dpsi;
rho_h_nir = (1 - sqrt(Plant.a_l_nir))/(1 + sqrt(Plant.a_l_nir));
rho_b_nir = 2*K_b./(K_b + K_d)*rho_h_nir;
rho_d_nir = rho_b_nir*contZ'*dpsi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 3: CALCULATE Kb FOR CURRENT ZENITH AND BEAM CANOPY REFLECTANCE WITHOUT SOIL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Kb for the current step in the model
phi1 = 0.5 - 0.633*Plant.x_lad - 0.33*Plant.x_lad^2;
phi2 = 0.877*(1 - 2*phi1);
G_z = phi1 + phi2*cos(BC.zenith);
K_b = min(G_z./cos(BC.zenith),20); % Direct beam extinction coefficient without scattering.
K_bp_par = sqrt(Plant.a_l_par)*K_b; % Diffuse extinction coefficient adjusted for scattering (Kd' in DPF(1997))
K_bp_nir = sqrt(Plant.a_l_nir)*K_b; % Diffuse extinction coefficient adjusted for scattering (Kd' in DPF(1997))

% Direct beam canopy reflectance without soil
rho_b_par = 2*K_b./(K_b + K_d)*rho_h_par;
rho_b_nir = 2*K_b./(K_b + K_d)*rho_h_nir;

% Canopy reflectance including soil
rho_cb_par = rho_b_par + (Soil.rho_g_par - rho_b_par).*exp(-2*K_bp_par*LAI);
rho_cb_nir = rho_b_nir + (Soil.rho_g_nir - rho_b_nir).*exp(-2*K_bp_nir*LAI);
rho_cd_par = rho_d_par + (Soil.rho_g_par - rho_d_par).*exp(-2*K_dp_par*LAI);
rho_cd_nir = rho_d_nir + (Soil.rho_g_nir - rho_d_nir).*exp(-2*K_dp_nir*LAI);

end