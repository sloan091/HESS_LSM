% =========================================================================
% Name   : Rad_Transfer.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This function calculates the shortwave (s) longwave (L) radiative forcing 
% at the sunlit and shaded big leaf and ground used used by the
% dual-source, two-big-leaf land surface model in Sloan et al. (2021). For 
% shortwave radiation, we use the Goudriaan and van Laar (1994) model 
% derived for sunlit and shaded big leaves as laid out by Depury and 
% Farquhar (1997) and Bonan (2019). The model performs calculations
% separately for photosynthetically active (PAR) and near infrared (NIR)
% radiation bands. For longwave radiation, we use the
% methodologies laid out in Dai (2004).  See Sect. S2.2 of Sloan et al.
% (2021) for full details.
%
% INPUTS
%   BC     - Boundary conditions (i.e., environmental forcings)
%   Plant  - Plant-specific parameters
%   Soil   - Soil-specific parameters
%   OS_DVs - Decision variables (T_l_sl, T_l_sh, T_g) for the LSM
%   
% OUTPUTS
%   Rad - Outputs from the radiative transfer model
%
% REFERENCES
%   (1) Goudriaan, J., & Laar, H. H. van. (1994). Modelling potential crop 
%   growth processes : textbook with exercises. Current Issues in 
%   Production Ecology Volume 2 (First). Wageningen: Springer Science 
%   and Business Media Dordrecht. https://doi.org/10.1007/978-94-011-0750-1
%
%   (2) De Pury, D. G. G., & Farquhar, G. D. (1997). Simple scaling of
%   photosynthesis from leaves to canopies without the errors of big-leaf
%   models. Plant, Cell and Environment, 20, 537–557. 
%   https://doi.org/10.1111/j.1365-3040.1997.00094.x
%
%   (3) Bonan, G. (2019). Climate Change and Terrestrial Ecosystem 
%   Modeling. Cambridge University Press.
%   https://doi.org/10.1017/9781107339217
%
%   (4) Dai, Y., et al. (2004). A Two-Big-Leaf Model for Canopy 
%   Temperature, Photosynthesis, and Stomatal Conductance. 
%   Journal of Climate, 17(12), 2281–2299. 
%   https://doi.org/10.1175/1520-0442(2004)017<2281:ATMFCT>2.0.CO;2
% =========================================================================

function Rad = Rad_Transfer(BC,Plant,Soil,OS_DVs)


% PART 1: EXTINCTION COEFFICIENTS, REFLECTANCES, AND LAI - (1-3)
%'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' 

% Calculate extinction coefficients and reflectances for PAR and NIR 
[Rad.K_b,Rad.K_bp_par,Rad.K_bp_nir,Rad.K_d,Rad.K_dp_par,Rad.K_dp_nir,...
    Rad.rho_cb_par,Rad.rho_cd_par,Rad.rho_cb_nir,Rad.rho_cd_nir] = ...
    Rad_K_rho(Plant,Soil,BC);

% Leaf area index [m^2 LA/m^2 GA]
LAI = Plant.LAI;

% Calculate sunlit and shaded LAI [m^2 LA/m^2 GA]
Rad.LAI_sl = (1 - exp(-Rad.K_b*LAI))./Rad.K_b;  
Rad.LAI_sh = LAI - Rad.LAI_sl;  

% Calculate fraction of sunlit and shaded leaves [-]
Rad.F_sl = Rad.LAI_sl./LAI;
Rad.F_sh = 1 - Rad.F_sl;


% PART 2: PAR SHORTWAVE RADIATION BALANCE - (1-3)
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''                                                               


% PART 2A: CALCULATE TOTAL CANOPY PAR ABSORPTION
%'''''''''''''''''''''''''''''''''''''''''''''''

% Total direct beam PAR absorbed by all leaves [W/m^2]
Rad.S_bt_par = (1 - Rad.rho_cb_par).*BC.S_ob_par.*...
    (1 - exp(-Rad.K_bp_par*LAI)); 

% Total diffuse PAR absorbed by all leaves [W/m^2]
Rad.S_dt_par = (1 - Rad.rho_cd_par)*BC.S_od_par*...
    (1 - exp(-Rad.K_dp_par*LAI));

% Total canopy PAR absorption [W/m^2]
Rad.S_c_par = Rad.S_bt_par + Rad.S_dt_par;


% PART 2B: CALCULATE SUNLIT BIG LEAF PAR ABSORPTION
%''''''''''''''''''''''''''''''''''''''''''''''''''

% Unscattered direct beam PAR absorbed by sunlit leaves [W/m^2]
Rad.S_sl_b_par = BC.S_ob_par.*Plant.a_l_par.*(1 - exp(-Rad.K_b*LAI)); 

% Scattered direct beam PAR absorbed by sunlit leaves [W/m^2]
Rad.S_sl_bs_par = BC.S_ob_par.*((1-Rad.rho_cb_par).*...
    (1 - exp(-(Rad.K_bp_par + Rad.K_b)*LAI)).*Rad.K_bp_par./...
    (Rad.K_bp_par + Rad.K_b) - Plant.a_l_par.*(1 - exp(-2*Rad.K_b*LAI))/2);

% Diffuse beam PAR absorbed by sunlit leaves [W/m^2]
Rad.S_sl_d_par = BC.S_od_par.*(1-Rad.rho_cd_par).*...
    (1 - exp(-(Rad.K_dp_par + Rad.K_b)*LAI)).*Rad.K_dp_par./...
    (Rad.K_dp_par + Rad.K_b);

% Total PAR absorbed by sunlit leaves [W/m^2]
Rad.S_sl_par = Rad.S_sl_b_par + Rad.S_sl_d_par + Rad.S_sl_bs_par;

% PART 1C: CALCULATE SHADED BIG LEAF PAR ABSORPTION
%''''''''''''''''''''''''''''''''''''''''''''''''''

% Total PAR absorbed by shaded leaves [W/m^2]
Rad.S_sh_par = Rad.S_c_par - Rad.S_sl_par;


% PART 2D: CALCULATE SOIL PAR ABSORPTION AND SUNLIT LEAF AREAS
%'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Total PAR absorbed by soil [W/m^2]
Rad.S_g_par = (BC.S_ob_par.*(1-Rad.rho_cb_par) + BC.S_od_par.*...
    (1-Rad.rho_cd_par)) - Rad.S_c_par;

% Total PAR reflected by canopy and soil [W/m^2]
Rad.S_r_par = BC.S_ob_par.*(Rad.rho_cb_par) + BC.S_od_par.*(Rad.rho_cd_par);


% PART 3: NIR SHORTWAVE RADIATION BALANCE - (1-3)
%'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' 

% PART 3A: CALCULATE TOTAL CANOPY NIR ABSORPTION
%'''''''''''''''''''''''''''''''''''''''''''''''

% Total direct beam NIR absorbed by all leaves [W/m^2]
Rad.S_bt_nir = (1 - Rad.rho_cb_nir).*BC.S_ob_nir.*...
    (1 - exp(-Rad.K_bp_nir.*LAI)); 

% Total diffuse NIR absorbed by all leaves [W/m^2]
Rad.S_dt_nir = (1 - Rad.rho_cd_nir).*BC.S_od_nir.*...
    (1 - exp(-Rad.K_dp_nir.*LAI));

% Total canopy NIR absorption [W/m^2]
Rad.S_c_nir = Rad.S_bt_nir + Rad.S_dt_nir;

% PART 3B: CALCULATE SUNLIT BIG LEAF NIR ABSORPTION
%''''''''''''''''''''''''''''''''''''''''''''''''''

% Unscattered direct beam NIR absorbed by sunlit leaves [W/m^2]
Rad.S_sl_b_nir = BC.S_ob_nir.*Plant.a_l_nir.*(1 - exp(-Rad.K_b.*LAI)); 

% Scattered direct beam NIR absorbed by sunlit leaves [W/m^2]
Rad.S_sl_bs_nir = BC.S_ob_nir.*((1-Rad.rho_cb_nir).*...
    (1 - exp(-(Rad.K_bp_nir + Rad.K_b).*LAI)).*Rad.K_bp_nir./...
    (Rad.K_bp_nir + Rad.K_b) - Plant.a_l_nir.*(1 - exp(-2.*Rad.K_b.*LAI))/2);

% Diffuse beam NIR absorbed by sunlit leaves [W/m^2]
Rad.S_sl_d_nir = BC.S_od_nir.*(1-Rad.rho_cd_nir).*...
    (1 - exp(-(Rad.K_dp_nir + Rad.K_b).*LAI)).*Rad.K_dp_nir./...
    (Rad.K_dp_nir + Rad.K_b);
 
% Total NIR absorbed by sunlit leaves [W/m^2]
Rad.S_sl_nir = Rad.S_sl_b_nir + Rad.S_sl_d_nir + Rad.S_sl_bs_nir;


% PART 3C: CALCULATE SHADED BIG LEAF NIR ABSORPTION
%''''''''''''''''''''''''''''''''''''''''''''''''''

% Total NIR absorbed by shaded leaves [W/m^2]
Rad.S_sh_nir = Rad.S_c_nir - Rad.S_sl_nir;

% PART 3D: CALCULATE SOIL NIR ABSORPTION AND REFLECTANCE
%'''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Total NIR absorbed by soil [W/m^2]
Rad.S_g_nir = (BC.S_ob_nir.*(1-Rad.rho_cb_nir) + BC.S_od_nir.*...
    (1-Rad.rho_cd_nir)) - Rad.S_c_nir;

% Total NIR reflected by canopy and soil [W/m^2]
Rad.S_r_nir = BC.S_ob_nir.*(Rad.rho_cb_nir) + BC.S_od_nir.*(Rad.rho_cd_nir);


% PART 4: LONGWAVE RADIATION BALANCE - (4)
%'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' 

% Stefan-Boltzmann constant [W/m^2/K]
sig = 5.67*10^(-8);    

% Fraction of longwave radiation absorbed by the canopy [-]
delt = 1 - exp(-LAI);  

% Total longwave radiation absorbed by sunlit leaves [W/m^2]
Rad.L_sl = (BC.L_in.*delt - 2.*sig.*(OS_DVs.T_l_sl + 273.15)^4.*delt +...
    sig.*(OS_DVs.T_g + 273.15)^4.*delt).*Rad.F_sl;

% Total longwave radiation absorbed by shaded leaves [W/m^2]
Rad.L_sh = (BC.L_in.*delt - 2.*sig.*(OS_DVs.T_l_sh + 273.15)^4.*delt +...
    sig.*(OS_DVs.T_g + 273.15)^4.*delt).*Rad.F_sh;

% Total longwave radiation absorbed by the ground [W/m^2]
Rad.L_g  = (1 - delt).*BC.L_in + Rad.F_sl.*delt.*sig.*(OS_DVs.T_l_sl +...
    273.15)^4 + Rad.F_sh.*delt.*sig.*(OS_DVs.T_l_sh + 273.15)^4 -...
    sig.*(OS_DVs.T_g + 273.15)^4;

% Total longwave radiation reflected by canopy and soil [W/m^2]
Rad.L_out = delt.*Rad.F_sl.*sig.*(OS_DVs.T_l_sl + 273.15)^4 +...
    delt.*Rad.F_sh.*sig.*(OS_DVs.T_l_sh + 273.15)^4 +...
    (1 - delt).*sig.*(OS_DVs.T_g + 273.15)^4;


% PART 5: NET RADIATION AND ENERGY BALANCE CHECK
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Net radiation absorbed of sunlit leaves (i.e., big leaf) [W/m^2]
Rad.Rn_sl = Rad.S_sl_par + Rad.S_sl_nir + Rad.L_sl;

% Net radiation absorbed of shaded leaves (i.e., big leaf) [W/m^2]
Rad.Rn_sh = Rad.S_sh_par + Rad.S_sh_nir + Rad.L_sh;

% Net radiation absorbed of ground (i.e., soil) [W/m^2]
Rad.Rn_g  = Rad.S_g_par + Rad.S_g_nir + Rad.L_g;

% Total measured incoming shortwave radiation [W/m^2]
Rad.S_in = BC.S_ob_par + BC.S_od_par + BC.S_ob_nir + BC.S_od_nir;

% Total calculated outgoing shortwave radiation [W/m^2]
Rad.S_out = Rad.S_r_par + Rad.S_r_nir;

% Total calculated net radiation of soil-canopy system [W/m^2]
Rad.Rn_all = Rad.S_in - Rad.S_out + BC.L_in - Rad.L_out;

% Check energy balance. This should only trip if there are erroneous
% parameter values or forcings that cause the underlying equations to blow
% up.
check = Rad.Rn_all - (Rad.Rn_sl + Rad.Rn_sh + Rad.Rn_g);
if abs(check) > 1e-10
    error('Energy Balance Failure')
else
end

end
