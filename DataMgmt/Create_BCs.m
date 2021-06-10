% =========================================================================
% Script Name   : Create_Outputs.m
% Author        : Brandon Sloan
% Start Date    : Aug 12, 2019
% Last Updated  : Aug 12, 2019
%
% Description   : This function creates the MATLAB table that will be used
% to store the model outputs.
%
%   INPUTS:
%   FluxData - A MATLAB table of the selected flux tower data.
%   Const    - Structure input of physical and site specific constants
%   Soil     - Structure input of the soil-specific parameters
%   Timestep - Number of the current timestep
%
%   OUTPUTS:
%   BC       - Structure output of the boundary conditions
% =========================================================================

function BC = Create_BCs(FluxData,Plant,Soil,Const,Flag)

% MISCELLANEOUS
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% % Current timestep
% BC.Step = Timestep;

% Flux tower measurement height [m]
BC.z = Const.z;  

% SOIL
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% Measured soil water content [-]
BC.theta_s = FluxData.SWC_1_5_1/100; 

% Calculated soil water potential [MPa]
BC.psi_s = theta2psi(Soil.psi_sat,Soil.b,BC.theta_s/Soil.theta_sat);  

% % Soil moisture correction factor for sunlit and shaded leaf [-]
% BC.beta_sl = Plant.beta_sl(Timestep);
% BC.beta_sh = Plant.beta_sh(Timestep);

% Soil moisture correction factor [-]
q = Plant.q;
if isequal(Flag.SMflag,1)
    BC.betas = ((Plant.psi_c - BC.psi_s)/(Plant.psi_c - Plant.psi_o))^q;
    if ~isreal(BC.betas)
        BC.betas = 0;
    else
        BC.betas = min(BC.betas,1);
    end
else
    BC.betas = 1;
end

BC.beta_sl = BC.betas;
BC.beta_sh = BC.betas;
% ATMOSPHERE
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

% RADIATION
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% Solar elevation angle [radians]
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