% =========================================================================
% Name   : createBCs.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This function creats environmental forcings or boundary conditions for
% the LSM in Sloan et al. (2021). See Fig. S1 and Sect S5 of Sloan et al.
% (2021) for further details on the forcing data taken from the Ameriflux
% US-Me2 site.
%
% INPUTS
%   FluxData - Environmental forcing data
%   Soil     - Soil-specific parameters
%   Const    - Physical constants
%
% OUTPUTS
%   BC - Structure array of required boundary conditions
% =========================================================================

function BC = createBCs(FluxData,Soil,Const)

% PART 1: MEASUREMENT HEIGHT
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Flux tower measurement height [m]
BC.z = Const.z;  

% PART 2: SUBSURFACE BOUNDARY CONDITIONS
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Measured soil water content [-]
BC.theta_s = FluxData.SWC_1_5_1/100; 

% Calculated soil water potential [MPa]
BC.psi_s = theta2psi(Soil.psi_sat,Soil.b,BC.theta_s/Soil.theta_sat);  

% PART 2: ATMOSPHERIC BOUNDARY CONDITIONS
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Friction velocity [m/s]
BC.u_star = FluxData.USTAR;

% Streamwise velocity at measurement height [m/s]
BC.U = FluxData.WS_F;   

% Air temperature [C]
BC.T_a = FluxData.TA_F;  

% Vapor pressure [Pa]
BC.e_a = CC(BC.T_a).*FluxData.RH/100;  

% Atmospheric pressure [Pa]
BC.P_atm = FluxData.PA_F*1000;     

% CO2 partial pressure of air [Pa]
BC.c_a = FluxData.CO2_F_MDS/10^6.*BC.P_atm;   

% RADIATION BOUNDARY CONDITIONS
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Solar elevation angle [radians]. NOTE: The inputs of SolarAxEl()
% corresponds to the US-Me2 site.  YOU MUST CHANGE THE INPUTS IF USING LSM
% AT A DIFFERENT LOCATION OTHERWISE THE RADIATIVE TRANSFER MODEL WILL BE
% WRONG.
[~,El] = SolarAzEl(FluxData.TIMESTAMP_START + 8/24,44.4523,-121.5574,1253);

% Solar zenith angle [radians]
BC.zenith = deg2rad(90-El);

% Direct beam PAR on horizontal plane above canopy [W/m^2]
BC.S_ob_par = (FluxData.PPFD_IN - FluxData.PPFD_DIF)/4.6;

% Diffuse beam PAR on horizontal plane above canopy [W/m^2]
BC.S_od_par = FluxData.PPFD_DIF/4.6; 

% Fraction of PAR that is diffuse
frac_dif = min(max(BC.S_od_par./(FluxData.PPFD_IN/4.6),0),1);

% Calculated direct beam NIR on horizontal plane above canopy [W/m^2]
BC.S_ob_nir = (FluxData.SW_IN_F - FluxData.PPFD_IN/4.6).*(1 - frac_dif);  

% Calculated diffuse beam NIR on horizontal plane above canopy [W/m^2]
BC.S_od_nir = (FluxData.SW_IN_F - FluxData.PPFD_IN/4.6).*frac_dif;   

% Incoming longwave radiation at canopy top [W/m^2]
BC.L_in = FluxData.LW_IN_F; 

% Ground heat flux [W/m^2]
BC.G = FluxData.G_F_MDS;                                                 

end