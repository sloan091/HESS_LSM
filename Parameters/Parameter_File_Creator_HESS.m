% =========================================================================
% Name   : Parameter_File_Creator_HESS.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This script creates the parameter file for running the LSM from Sloan et
% al. (2021) for the US-Me2 Ameriflux site.  See Sect. S3 and S4 of Sloan
% et al. (2021) for futher details on these parameter values.  Below I have
% included a list of references for these parameter values that are cited
% in the code. Overall, this code creates 4 structure arrays used by my LSM
% code.
%
% REFERENCES
%   (1) Ameriflux Base Measurement Height file named 
%   "BASE_MeasurementHeight_20191025.csv"
%
%   (2) Oleson, K. W. et al. (2018). Technical Description of the version 5
%   of the Community Land Model (CLM). 
%
%   (3) Ameriflux Site Data file named "AMF_US-Me2_BIF_LATEST.xlsx"
%
%   (4) Monteith, J., & Unsworth, M. (2013). Principles of Environmental
%   Physics: Plants, Animals, and the Atmosphere: Fourth Edition.
%   https://doi.org/10.1016/C2010-0-66393-0
%
%   (5) Schwarz, P. A. et al. (2004). Climatic versus biotic constraints on
%   carbon and water fluxes in seasonally drought-affected ponderosa 
%   pine ecosystems. Global Biogeochemical Cycles, 18(4), 1–17. 
%   https://doi.org/10.1029/2004GB002234
%
%   (6) Law, B. E. et al. (2001). Carbon storage and fluxes in ponderosa
%   pine forests at different developmental stages. Global Change Biology,
%   7(7), 755–777. https://doi.org/10.1046/j.1354-1013.2001.00439.x
%
%   (7) Clapp, R. B., & Hornberger, G. M. (1978). Empirical equations for 
%   some soil hydraulic properties. Water Resources Research, 14(4), 
%   601–604. https://doi.org/10.1029/WR014i004p00601
%
%   (8) Jackson, R. B., Mooney, H. A., & Schulze, E. D. (1997). A global
%   budget for fine root biomass, surface area, and nutrient contents. 
%   Proceedings of the National Academy of Sciences of the United States 
%   of America, 94(14), 7362–7366. https://doi.org/10.1073/pnas.94.14.7362
%
%   (9) Irvine, J. etal. (2004). Age-related changes in ecosystem...
%   structure and function and effects on water and carbon exchange
%   in ponderosa pine. Tree Physiology (Vol. 24). Oxford Academic. 
%   https://doi.org/10.1093/TREEPHYS/24.7.753
%
%   (10) Choat, B. et al. (2012). Global convergence in the vulnerability
%   of forests to drought. Nature, 491(7426), 752–755.
%   https://doi.org/10.1038/nature11688
%
%   (11) Maherali, H., et al. (2000). Xylem conductivity and vulnerability
%   to cavitation of ponderosa pine growing in contrasting climates. 
%   Tree Physiology, 20(13), 859–867. 
%   https://doi.org/10.1093/treephys/20.13.859
%
%   (12) DeLucia, E. H., & Heckathorn, S. A. (1989). The effect of soil 
%   drought on water?use efficiency in a contrasting Great Basin desert 
%   and Sierran montane species. Plant, Cell & Environment, 12(9), 935–940. 
%   https://doi.org/10.1111/j.1365-3040.1989.tb01973.x
%
%   (13) Feng, X., et al. (2018). The ecohydrological context of drought 
%   and classification of plant responses. Ecology Letters. 
%   https://doi.org/10.1111/ele.13139
%
%   (14) Ameriflux Me2 Half-Hourly flux data file named
%   "AMF_US-Me2_BASE_HH_9-5.csv".
%
%   (15) Dai, Y., et al. (2004). A Two-Big-Leaf Model for Canopy 
%   Temperature, Photosynthesis, and Stomatal Conductance. 
%   Journal of Climate, 17(12), 2281–2299. 
%   https://doi.org/10.1175/1520-0442(2004)017<2281:ATMFCT>2.0.CO;2
%
% =========================================================================

clc
clear

% PHYSICAL CONSTANTS
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Specific heat at constant pressure [J/kg/K]
Const.c_p = 1004;     

% Specific gas constant for water vapor [J/k/kg]
Const.R_v = 461;  

% von Karmen constant
Const.k = 0.4;   

% Gravitational acceleration (m/s^2)
Const.g = 9.8;    

% Average density of air (kg/m^3)
Const.rho_a = 1.2;

% Average density of water (kg/m^3)
Const.rho_w = 1000;  

% Kinematic viscosity of water [m^2/s]
Const.nu = 1.5*10^-5;   

% Molecular weight of water [kg/mol]
Const.M_w = 0.018; 

% Molecular weight of air [kg/mol]
Const.M_a = 0.029;        

% Latent heat of vaporization [J/kg]
Const.L_vap = 2.5*10^6;   

% Universal gas constant [J/K/kmol]
Const.R_g = 8314;   

% Measurement height of flux tower [m] - (1)
Const.z = 33;                    

% PLANT TRAITS
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% PHYSIOLOGY
% """"""""""

% Characteristic leaf length [m] - (2)
Plant.d_l = 0.04;    

% Leaf area index [m^2 LA/m^2 GA] - Calibrated
Plant.LAI = 3.232;    

% Stem area index [m^2 SA/m^2 GA] - Educated guess
Plant.SAI = 0.5;     

% Vegetation height [m] - (3)
Plant.h_v = 18;   

% Number of sides of leaf that have stomata. 1 for hypo- or hyperstomatous
% and 2 for amphistomatous. Always have set to 1 since most stomatal
% conductance measurements do not differentiate between one and two-sided
% stomatal coverage.
Plant.ns = 1; 

% Canopy momentum roughness length for coniferous forest [m] - (4)
Plant.z_0m_c = 1;   


% COUPLED STOMATAL CONDUCTANCE PHOTOSYNTHESIS MODEL
% """""""""""""""""""""""""""""""""""""""""""""""""

% Medlyn slope parameter [kPa^0.5] - Calibrated
Plant.g_1 = 0.8765;  

% Medlyn intercept parameter [moles/m^2/s] - (2)
Plant.g_o = 10^(-4); 

% Maximum rate of carboxylation [micromoles/m^2 LA/s] - Calibrated
Plant.V_max25 = 122.38;  

% Maximum electron transport rate [micromoles electrons/m^2 LA/s] -
% Assumption proposed by (15).
Plant.J_max25 = 122.38*2.1;  

% Quantum efficiency of photosystem II 
% [micromoles electrons/ micromoles photons] - (2)
Plant.Sig_psii = 0.85;  

% Curvature parameter for light-limited assimilation [-] - (2)
Plant.Theta_psii = 0.7; 

% Extinction coefficient of leaf nitrogen allocation [-] - (2)
Plant.K_n = 0.7;     


% RADIATIVE TRANSPORT
% """""""""""""""""""

% Leaf angle distribution parameter [-] - Calibrated
Plant.x_lad = 0.1059;     

 % Leaf absorptivity in PAR waveband [-] - Calibrated
Plant.a_l_par = 0.735;   

 % Leaf absorptivity in NIR waveband [-] - Calibrated
Plant.a_l_nir = 0.431;   



% EMPIRICAL STATIC BETA FOR TRANSPIRATION DOWNREGULATION
% """"""""""""""""""""""""""""""""""""""""""""""""""""""

% Soil water potential at 50% loss of stomatal closure [MPa] for sunlit big
% leaf - Calibrated
Plant.psi_s_50_sl = -0.695;   

% Soil water potential at 50% loss of stomatal closure [MPa] for shaded big
% leaf - Calibrated
Plant.psi_s_50_sh = -0.7791;   

% Stomatal sensitivity parameter for sunlit big leaf - Calibrated
Plant.b_s_sl = 2.74;   

% Stomatal sensitivity parameter for shaded big leaf - Calibrated
Plant.b_s_sh = 4.055;   


% EMPIRICAL DYNAMIC BETA FOR TRANSPIRATION DOWNREGULATION
% """"""""""""""""""""""""""""""""""""""""""""""""""""""

% Slope parameter for linear dependence of psi_s_50 on T_ww [MPa*day/mm] 
% for sunlit big leaf - Calibrated
Plant.psi_s_50_slope_sl = 0.06675;   

% Slope parameter for linear dependence of psi_s_50 on T_ww [MPa*day/mm] 
% for shaded big leaf - Calibrated
Plant.psi_s_50_slope_sh = 0.08104;   

% Intercept parameter for linear dependence of psi_s_50 on T_ww [MPa] 
% for sunlit big leaf - Calibrated
Plant.psi_s_50_int_sl = -0.9114;   

% Intercept parameter for linear dependence of psi_s_50 on T_ww [MPa*day/mm] 
% for shaded big leaf - Calibrated
Plant.psi_s_50_int_sh = -0.8912;   

% Slope parameter for linear dependence of b_s on T_ww [day/mm] 
% for sunlit big leaf - Calibrated
Plant.b_s_slope_sl = -0.4412;   

% Slope parameter for linear dependence of b_s on T_ww [day/mm] 
% for shaded big leaf - Calibrated
Plant.b_s_slope_sh = -0.8954;   

% Intercept parameter for linear dependence of b_s on T_ww [-] 
% for sunlit big leaf - Calibrated
Plant.b_s_int_sl = 4.379;   

% Intercept parameter for linear dependence of b_s on T_ww [-] 
% for shaded big leaf - Calibrated
Plant.b_s_int_sh = 5.444; 


% PLANT HYDRAULICS MODEL
% """"""""""""""""""""""

% Rooting area index [m^2 root surface/m^2 ground] Temperate Evergreen - (8)
Plant.RAI =  11;     

% Fine root diameter [m] - (8)
Plant.d_r = 0.5/1000;      

% Effective rooting depth [m] - (5)
Plant.Z_r = 1.1;                    

% Sapwood-area specific hydraulic conductivity [kg/m/s/MPA] - Calibrated
Plant.k_sap = 1.33; 

% Sapwood area index [m^2 SapA/m^2 GA] - (9)
Plant.SapAI   =  20e-4;  

% Water potential at 50% loss of conductivity [MPa] - (11)
Plant.psi_x_50 = -2.6;    

% Xylem vulnerability curve parameter [-] - Calibrated
Plant.a = 0.54;       

% Water potential at stomatal closure [MPa] - (12)
Plant.psi_l_50 = -1; 

% Stomatal sensitivity parameter [-] - Calibrated
Plant.b_l = 5;                

%% SOIL PROPERTIES
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% RETENTION CURVE
% """""""""""""""

% Saturated hydraulic conductivity [m/day] - Calibrated
Soil.ksat = 0.81;   

% Soil water potential near saturation [MPa] - Calibrated
Soil.psi_sat = -5.5*10^(-3); 

% Brooks-Corey soil water retention exponent [-] - Calibrated
Soil.b = 3.86;    

% Accounts for root extension in drying soil for PHM [-] - Originally set
% to 4 based on (13), but then set to 0 based on calibration.  See Sect. S4
% in Sloan et al. (2021) for details.
Soil.d = 0;                   

% Saturated water content or porosity [-]. There is debate on this value.
% Source (5) claims a value of 0.57, however, the data I am using (14)
% shows only a max value of 0.3. I have stuck with the reported value in
% (5).
Soil.theta_sat = 0.57; 

% Water content at which soil restricts conductance [-] - Calibrated
Soil.theta_i = 0.57;          


% RADIATIVE TRANSFER
% """"""""""""""""""

% Ground reflectance for near infrared [-] - (2)
Soil.rho_g_nir = 0.2;  

% Ground reflectance for photosynthetically active radiation [-] - (2)
Soil.rho_g_par = 0.1;         


%% SIMULATION SETTINGS
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% PHOTOSYNTHESIS
% """"""""""""""

% 1: Photosynthetis model with temperature dependence; 0: Photosynthesis 
% model independent of temperature
Flag.Tflag = 0;  

% 1: FcVB C3 photosynthesis model; 2: C4 Collatz photosynthesis model
Flag.Cflag = 1;  


% TRANSPIRATION DOWNREGULATION
% """"""""""""""""""""""""""""

% 1: Static beta function dependent only on soil water potential; 
% 2: Dynamic beta function with additoinal dependence on well-watered
% transpiration rate (T_ww); 0: Do not perform LSM sim with beta
% downregulation.
Flag.Betaflag = 2; 

% 1: Use Plant Hydraulics Model (PHM); 0: Do not perform LSM sim with PHM
% downregulation.
Flag.PHMflag = 0; 

% 1: My downregulation scheme with fixed T_ww; 2 CLM scheme with adaptive
% T_ww.
Flag.DownRegMethod = 1; 

save('Parameters_HESS_Final_beta_dyn_Weight.mat','Const','Plant','Soil','Flag')