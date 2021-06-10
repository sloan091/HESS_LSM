% =========================================================================
% Script Name   : Rad_Transfer.m
% Author        : Brandon Sloan
% Start Date    : Mar 25, 2019
% Last Updated  : Mar 25, 2019
%
% Description   : This function calculates the shortwave (solar) radiative forcing at the leaf
% level and also the fraction of sunlit and shaded leaves used to scale up the
% results to the canopy level. The code largely follows Depury and Farquhar
% (1997) as well as Goudriaan and Van Laar (1994) with some help from
% Norman (1998).  It also calculates the thermal radiation fluxes using the
% methodlogies laid out in Dai (2004)
%
%   INPUTS:
%   S_ob   - Observed direct beam radiation [W/m^2 or mole PAR/m^2/s]
%   S_od   - Observed diffuse radiation [W/m^2 or mole PAR/m^2/s]
%   LAI    - Leaf area index [m^2 leaf area/m^2 ground area]
%   a_l    - Leaf absorptivity at radiation waveband [-]
%   x      - Leaf angle distribution parameter [-]
%   BC.zenith - Solar BC.zenith angle [-]
%   BC.L_in   - Incoming longwave radiation at canopy top [W/m^2]
%   SEB_DVs.T_l_sl - Sunlit leaf temperature [degrees C]
%   SEB_DVs.T_l_sh - Shaded leaf temperature [degrees C]
%   SEB_DVs.T_g    - Soil temperature [degrees C]
%   
%
%   OUTPUTS:
%   S_sl - Shortwave radiation absorbed by sunlit leaves [W/m^2 or mole PAR/m^2/s]
%   S_sh - Shortwave radiation absorbed by shaded leaves [W/m^2 or mole PAR/m^2/s]
%   S_c  - Total canopy absorption of irradiance not separating sunlit and shaded [W/m^2 or mole PAR/m^2/s]
%   S_r  - Total shortwave radiation reflected [W/m^2 or mole PAR/m^2/s]
%   S_g  - Shortwave radiation absorbed by ground [W/m^2 or mole PAR/m^2/s]
%   F_sl - Fraction of LAI that is sunlit [-]
%   F_sl - Fraction of LAI that is shaded [-]
%   L_sl - Longwave radiation absorbed by sunlit leaves [W/m^2]
%   L_sh - Longwave radiation absorbed by shaded leaves [W/m^2]
%   L_g  - Longwave radiation absorbed by ground [W/m^2]
%   L_out- Longwave radiation leaving the canopy [W/m^2]
% =========================================================================

function Rad = Rad_Transfer(BC,Plant,Soil,SEB_DVs,Flag)

% Grab the correct LAI 
if isequal(Flag.LAIflag,1)
    LAI = Plant.LAI(BC.Step);
else
    LAI = Plant.LAI;
end                                                                   

[Rad.K_b,Rad.K_bp_par,Rad.K_bp_nir,Rad.K_d,Rad.K_dp_par,Rad.K_dp_nir,Rad.rho_cb_par,Rad.rho_cd_par,...
    Rad.rho_cb_nir,Rad.rho_cd_nir] = Rad_K_rho(Plant,Soil,BC,Flag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1: SHORTWAVE RADIATION BALANCE FOR PAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A) CALCULATE TOTAL CANOPY ABSORPTION OF IRRADIANCE FROM DEPURY (1997)
% Total beam absorbed by all leaves (scattered and unscattered)
Rad.S_bt_par = (1 - Rad.rho_cb_par).*BC.S_ob_par.*(1 - exp(-Rad.K_bp_par*LAI)); 
% Total diffues absorbed by all leaves
Rad.S_dt_par = (1 - Rad.rho_cd_par)*BC.S_od_par*(1 - exp(-Rad.K_dp_par*LAI));
% Total Canopy absorption of irradiance not separating sunlit and shaded (Eqn 13)
Rad.S_c_par = Rad.S_bt_par + Rad.S_dt_par;

% B) CALCULATE SUNLIT CANOPY ABSORPTION OF IRRADIANCE FROM DEPURY (1997)
% Direct beam absorption (Eqn 20b)
Rad.S_sl_b_par = BC.S_ob_par.*Plant.a_l_par.*(1 - exp(-Rad.K_b*LAI)); 
% Diffuse beam absorption (Eqn 20c)
Rad.S_sl_d_par = BC.S_od_par.*(1-Rad.rho_cd_par).*(1 - exp(-(Rad.K_dp_par + Rad.K_b)*LAI))...
    .*Rad.K_dp_par./(Rad.K_dp_par + Rad.K_b);
% Scattered direct beam absorption (Eqn 20d)
Rad.S_sl_bs_par = BC.S_ob_par.*((1-Rad.rho_cb_par).*(1 - exp(-(Rad.K_bp_par + Rad.K_b)*LAI))...
    .*Rad.K_bp_par./(Rad.K_bp_par + Rad.K_b) - Plant.a_l_par.*(1 - exp(-2*Rad.K_b*LAI))/2); 
% Total sunlit canopy absorption of irradiance (Eqn 20a)
Rad.S_sl_par = Rad.S_sl_b_par + Rad.S_sl_d_par + Rad.S_sl_bs_par;


% C) CALCULATE SHADED CANOPY ABSORPTION OF IRRADIANCE FROM DEPURY (1997)
% Total sunlit canopy absorption of irradiance (Eqn 21)
Rad.S_sh_par = Rad.S_c_par - Rad.S_sl_par;

% D) CALCULATE SHORTWAVE RADIATION ABSORBED BY SOIL
Rad.S_g_par = (BC.S_ob_par.*(1-Rad.rho_cb_par) + BC.S_od_par.*(1-Rad.rho_cd_par)) - Rad.S_c_par;

% E) CALCULATE TOTAL CANOPY-SOIL SYSTEM REFLECTANCE
Rad.S_r_par = BC.S_ob_par.*(Rad.rho_cb_par) + BC.S_od_par.*(Rad.rho_cd_par);

% F) CALCULATE SUNLIT AND SHADED FRACTIONS OF CANOPY
% Sunlit LAI and shaded LAI
Rad.LAI_sl = (1 - exp(-Rad.K_b*LAI))./Rad.K_b;  % Eqn 15.23 in Norman (1998) 
Rad.LAI_sh = LAI - Rad.LAI_sl;  

% Fraction of sunlit and shaded leaves
Rad.F_sl = Rad.LAI_sl./LAI;
Rad.F_sh = 1 - Rad.F_sl;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 2: SHORTWAVE RADIATION BALANCE FOR NIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A) CALCULATE TOTAL CANOPY ABSORPTION OF IRRADIANCE FROM DEPURY (1997)
% Total beam absorbed by all leaves (scattered and unscattered)
Rad.S_bt_nir = (1 - Rad.rho_cb_nir).*BC.S_ob_nir.*(1 - exp(-Rad.K_bp_nir.*LAI)); 
% Total diffues absorbed by all leaves
Rad.S_dt_nir = (1 - Rad.rho_cd_nir).*BC.S_od_nir.*(1 - exp(-Rad.K_dp_nir.*LAI));
% Total Canopy absorption of irradiance not separating sunlit and shaded (Eqn 13)
Rad.S_c_nir = Rad.S_bt_nir + Rad.S_dt_nir;

% B) CALCULATE SUNLIT CANOPY ABSORPTION OF IRRADIANCE FROM DEPURY (1997)
% Direct beam absorption (Eqn 20b)
Rad.S_sl_b_nir = BC.S_ob_nir.*Plant.a_l_nir.*(1 - exp(-Rad.K_b.*LAI)); 
% Diffuse beam absorption (Eqn 20c)
Rad.S_sl_d_nir = BC.S_od_nir.*(1-Rad.rho_cd_nir).*(1 - exp(-(Rad.K_dp_nir + Rad.K_b).*LAI))...
    .*Rad.K_dp_nir./(Rad.K_dp_nir + Rad.K_b);
% Scattered direct beam absorption (Eqn 20d)
Rad.S_sl_bs_nir = BC.S_ob_nir.*((1-Rad.rho_cb_nir).*(1 - exp(-(Rad.K_bp_nir + Rad.K_b).*LAI))...
    .*Rad.K_bp_nir./(Rad.K_bp_nir + Rad.K_b) - Plant.a_l_nir.*(1 - exp(-2.*Rad.K_b.*LAI))/2); 
% Total sunlit canopy absorption of irradiance (Eqn 20a)
Rad.S_sl_nir = Rad.S_sl_b_nir + Rad.S_sl_d_nir + Rad.S_sl_bs_nir;


% C) CALCULATE SHADED CANOPY ABSORPTION OF IRRADIANCE FROM DEPURY (1997)
% Total sunlit canopy absorption of irradiance (Eqn 21)
Rad.S_sh_nir = Rad.S_c_nir - Rad.S_sl_nir;

% D) CALCULATE SHORTWAVE RADIATION ABSORBED BY SOIL
Rad.S_g_nir = (BC.S_ob_nir.*(1-Rad.rho_cb_nir) + BC.S_od_nir.*(1-Rad.rho_cd_nir)) - Rad.S_c_nir;

% E) CALCULATE TOTAL CANOPY-SOIL SYSTEM REFLECTANCE
Rad.S_r_nir = BC.S_ob_nir.*(Rad.rho_cb_nir) + BC.S_od_nir.*(Rad.rho_cd_nir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 3: TERRESTRIAL (LONGWAVE) RADIATION PARTITIONING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate sunlit and shaded radiation partitioning using eqns 18a and 18b
% in Dai (2004).  Basically just splits the radiation using the sunlit and
% shaded fractions and the respective temperatures.  This also uses a bit
% of Sellers (1996) and 

sig = 5.67*10^(-8);    % Stefan-Boltzmann constant [W/m^2/K]
delt = 1 - exp(-LAI);  % This number gives the fraction of longwave radiation absorbed by the canopy (Eqn 18a in Dai (2004))
Rad.L_sl = (BC.L_in.*delt - 2.*sig.*(SEB_DVs.T_l_sl + 273.15)^4.*delt + sig.*(SEB_DVs.T_g + 273.15)^4.*delt).*Rad.F_sl;
Rad.L_sh = (BC.L_in.*delt - 2.*sig.*(SEB_DVs.T_l_sh + 273.15)^4.*delt + sig.*(SEB_DVs.T_g + 273.15)^4.*delt).*Rad.F_sh;
Rad.L_g  = (1 - delt).*BC.L_in + Rad.F_sl.*delt.*sig.*(SEB_DVs.T_l_sl + 273.15)^4 + Rad.F_sh.*delt.*sig.*(SEB_DVs.T_l_sh + 273.15)^4 - sig.*(SEB_DVs.T_g + 273.15)^4;
Rad.L_out = delt.*Rad.F_sl.*sig.*(SEB_DVs.T_l_sl + 273.15)^4 + delt.*Rad.F_sh.*sig.*(SEB_DVs.T_l_sh + 273.15)^4 + (1 - delt).*sig.*(SEB_DVs.T_g + 273.15)^4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 4: CALCULATE NET RADIATIONS AND CHECK ENERGY BALANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rad.Rn_sl = Rad.S_sl_par + Rad.S_sl_nir + Rad.L_sl;  % Big sunlit leaf
Rad.Rn_sh = Rad.S_sh_par + Rad.S_sh_nir + Rad.L_sh;  % Big shaded leaf
Rad.Rn_g  = Rad.S_g_par + Rad.S_g_nir + Rad.L_g;     % Ground (Soil)

% Conservation of Energy Check
Rad.S_in = BC.S_ob_par + BC.S_od_par + BC.S_ob_nir + BC.S_od_nir;
Rad.S_out = Rad.S_r_par + Rad.S_r_nir;
Rad.Rn_all = Rad.S_in - Rad.S_out + BC.L_in - Rad.L_out;
check = Rad.Rn_all - (Rad.Rn_sl + Rad.Rn_sh + Rad.Rn_g);
if abs(check) > 1e-10
    error('Energy Balance Failure')
else
end

end
