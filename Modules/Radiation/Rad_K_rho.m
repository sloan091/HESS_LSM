% =========================================================================
% Name   : Rad_K_rho.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This function calculates the direct and diffuse beam extinction
% coefficients as well as reflectances used by the radiative transfer model
% in Sloan et al. (2021). The formulations used here follow closely
% Goudriaan (1977), Goudriaan and van Laar (1994) and Bonan (2019). See
% Sect. S2.2 of Sloan et al. (2021) for full details.
%
% INPUTS
%   Plant  - Plant-specific parameters
%   Soil   - Soil-specific parameters
%   BC     - Boundary conditions (i.e., environmental forcings)
%   
% OUTPUTS
%   K_b        - Direct beam extinction coefficient without scattering [-]
%   K_bp_par   - Direct beam PAR extinction coefficient with scattering [-]
%   K_bp_nir   - Direct beam NIR extinction coefficient with scattering [-]
%   K_d        - Diffuse beam extinction coefficient without scattering [-]
%   K_dp_par   - Diffuse beam PAR extinction coefficient with scattering [-]
%   K_dp_nir   - Diffuse beam NIR extinction coefficient with scattering [-]
%   rho_cb_par - Direct beam PAR soil-canopy reflectance [-]
%   rho_cd_par - Diffuse beam PAR soil-canopy reflectance [-]
%   rho_cb_nir - Direct beam NIR soil-canopy reflectance [-]
%   rho_cd_nir - Diffuse beam NIR soil-canopy reflectance [-]
%
% REFERENCES
%   (1) Goudriaan, J. (1977). Crop micrometeorology: a simulation 
%   study. Wageningen.
%
%   (1) Goudriaan, J., & Laar, H. H. van. (1994). Modelling potential crop 
%   growth processes : textbook with exercises. Current Issues in 
%   Production Ecology Volume 2 (First). Wageningen: Springer Science 
%   and Business Media Dordrecht. https://doi.org/10.1007/978-94-011-0750-1
%
%   (3) Bonan, G. (2019). Climate Change and Terrestrial Ecosystem 
%   Modeling. Cambridge University Press.
%   https://doi.org/10.1017/9781107339217
% =========================================================================
function [K_b,K_bp_par,K_bp_nir,K_d,K_dp_par,K_dp_nir,rho_cb_par,...
    rho_cd_par,rho_cb_nir,rho_cd_nir] = Rad_K_rho(Plant,Soil,BC)


% PART 1: CALCULATE DIFFUSE BEAM EXTINCTION COEFFICIENTS AND REFLECTANCES
%'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' 

% PART 1a: DIFFUSE BEAM EXTINCTION COEFFICIENTS
%''''''''''''''''''''''''''''''''''''''''''''''

% Leaf area index [m^2 LA/m^2 GA]
LAI = Plant.LAI;

% Zenith angle increment [rad]
dpsi = 0.01;

% Range of zenith angle to integrate extinction coefficient over [rad]
zen_range = 0:dpsi:pi()/2;

% Ross-Goudrian Function (next 3 lines)
phi1 = 0.5 - 0.633*Plant.x_lad - 0.33*Plant.x_lad^2;
phi2 = 0.877*(1 - 2*phi1);
G_z = phi1 + phi2*cos(zen_range);

% Direct beam extinction coefficient without scattering at each zenith
% angle [-]
K_b = min(G_z./cos(zen_range),20); 

% Direct beam transmission coefficient for a black leaf at the specified
% cumulative LAI [-]
tau_b = exp(-K_b.*LAI);

% Contribution of each sky region to the diffuse transmission coefficient
% [-]
contZ = 2.*sin(zen_range).*cos(zen_range);

% Diffuse beam transmission coefficient for the hemisphere [-]
tau_d = tau_b*contZ'*dpsi;

% Diffuse beam extinction coefficient without scattering [-]
K_d = -log(tau_d)/LAI;

% Diffuse beam PAR extinction coefficient with scattering [-]
K_dp_par = sqrt(Plant.a_l_par)*K_d; 

% Diffuse beam NIR extinction coefficient with scattering [-]
K_dp_nir = sqrt(Plant.a_l_nir)*K_d;


% PART 1B: DIFFUSE SOIL-PLANT REFLECTANCES
%'''''''''''''''''''''''''''''''''''''''''

% PAR reflectance of a thick canopy (i.e., light does not reach soil)
% with horizontal leaves [-]
rho_h_par = (1 - sqrt(Plant.a_l_par))/(1 + sqrt(Plant.a_l_par));

% Direct beam PAR reflectance for a thick canopy at every zenith angle [-]
rho_b_par = 2*K_b./(K_b + K_d)*rho_h_par;

% Diffuse beam PAR reflectance for a thick canopy [-]
rho_d_par = rho_b_par*contZ'*dpsi;

% Diffuse beam NIR soil-canopy reflectance [-]
rho_cd_par = rho_d_par + (Soil.rho_g_par - rho_d_par).*...
    exp(-2*K_dp_par*LAI);

% NIR reflectance of a thick canopy (i.e., light does not reach soil)
% with horizontal leaves [-]
rho_h_nir = (1 - sqrt(Plant.a_l_nir))/(1 + sqrt(Plant.a_l_nir));

% Direct beam NIR reflectance for a thick canopy at every zenith angle [-]
rho_b_nir = 2*K_b./(K_b + K_d)*rho_h_nir;

% Diffuse beam NIR reflectance for a thick canopy [-]
rho_d_nir = rho_b_nir*contZ'*dpsi;

% Diffuse beam NIR soil-canopy reflectance [-]
rho_cd_nir = rho_d_nir + (Soil.rho_g_nir - rho_d_nir).*...
    exp(-2*K_dp_nir*LAI);


% PART 2: CALCULATE DIRECT BEAM EXTINCTION COEFFICIENTS AND REFLECTANCES
%'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' 

% Ross-Goudrian Function for current time step's zenith angle (next 3 lines)
phi1 = 0.5 - 0.633*Plant.x_lad - 0.33*Plant.x_lad^2;
phi2 = 0.877*(1 - 2*phi1);
G_z = phi1 + phi2*cos(BC.zenith);

% Direct beam extinction coefficient without scattering [-]
K_b = min(G_z./cos(BC.zenith),20);

% Direct beam PAR extinction coefficient with scattering [-]
K_bp_par = sqrt(Plant.a_l_par)*K_b;

% Direct beam NIR extinction coefficient with scattering [-]
K_bp_nir = sqrt(Plant.a_l_nir)*K_b;

% Direct beam PAR reflectance of a thick canopy (i.e., light does not 
% reach soil)
rho_b_par = 2*K_b./(K_b + K_d)*rho_h_par;

% Direct beam PAR soil-canopy reflectance [-]
rho_cb_par = rho_b_par + (Soil.rho_g_par - rho_b_par).*exp(-2*K_bp_par*LAI);

% Direct beam NIR reflectance of a thick canopy (i.e., light does not 
% reach soil)
rho_b_nir = 2*K_b./(K_b + K_d)*rho_h_nir;

% Direct beam NIR soil-canopy reflectance [-]
rho_cb_nir = rho_b_nir + (Soil.rho_g_nir - rho_b_nir).*exp(-2*K_bp_nir*LAI);

end