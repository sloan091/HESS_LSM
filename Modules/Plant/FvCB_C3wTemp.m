% =========================================================================
% Name   : FvCB_C3wTemp.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% Description   : This function calculates the net CO2 assimilation using
% the photosynthesis model by Farquhar, von Caemmerer and Berry (1980) for
% C3 plants. Here, we follow the formulation and assumptions laid out by
% Oleson et al. (2013) including temperature dependence, which was ignored 
% for simplicity in Sloan et al. (2021). See Sect. S2.4.2 for full details.
% Area units are specified by LA for per unit leaf area and GA for per unit
% ground area
%
% INPUTS
%   V_max      - Maximum scaled Rubisco carboxylation rate 
%                [micromoles CO2 /m^2 LA/s]
%   J_max      - Maximum scaled electron transport rate 
%                [micromoles electron s /m^2 LA/s] 
%   Sig_psii   - Quantum efficiency of photosystem II 
%                [micromoles electrons/ micromoles photons]
%   Q          - Photosynthetically Active Radiation [W/m^2 LA]
%   Theta_psii - Curvature parameter for light-limited assimilation
%   P_atm      - Atmospheric pressure [Pa]
%   c_i        - CO2 partial pressure inside the leaf [Pa]
%   T_l        - Leaf temperature [C]
%   R_g        - Universal gas constant [J/K/kmol]
%
% OUTPUTS
%   A_n   - Net CO2 assimilation rate [micromoles CO2/m^2 LA/s]
%   A     - CO2 assimilation rate [micromoles CO2/m^2 LA/s]
%   R_d   - Dark respiration rate [micromoles CO2/m^2 LA/s]
%   ContA - Which process is limiting: 1 Rubisco, 2 light, and 3 production
%
% REFERENCES
%   (1) Farquhar, G. D., von Caemmerer, S., & Berry, J. A. (1980). 
%   A biochemical model of photosynthetic CO2 assimilation in leaves of C3 
%   species. Planta. https://doi.org/10.1007/BF00386231
%
%   (2) Oleson, K. W. et al. (2013). Technical Description of version 4.5 
%   of the Community Land Model (CLM). NCAR Technical Notes. 
%   https://doi.org/10.5065/D6RR1W7M
%
% =========================================================================

function [A_n,A,R_d,ContA] = FvCB_C3wTemp(V_max,J_max,Sig_psii,Q,Theta_psii,...
    P_atm,c_i,T_l,R_g)

% PART 1: DEFINE TEMPERATURE DEPENDENT PARAMETERS
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% V_max activation energy [J/mol]
delHa_V = 65330;  
% V_max deactivation energy [J/mol]
delHd_V = 149250; 
% V_max entropy term [J/mol/K]
delS_V = 485;

% J_max activation energy [J/mol]
delHa_J = 43540;
% J_max deactivation energy [J/mol]
delHd_J = 152040; 
% J_max entropy term [J/mol/K]
delS_J = 495;    

% K_c activation energy [J/mol]
delHa_Kc = 79430;  
% K_o activation energy [J/mol]
delHa_Ko = 36380;

% T_p activation energy [J/mol]
delHa_Tp = 65330;
% T_p deactivation energy [J/mol]
delHd_Tp = 149250; 
% T_p entropy term [J/mol/K]
delS_Tp = 485;

% R_d activation energy [J/mol]
delHa_Rd = 46930;
% R_d activation energy [J/mol]
delHd_Rd = 150650;
% R_d entropy term [J/mol/K]
delS_Rd = 490;     

% Gamma activation energy [J/mol]
delHa_G = 37830;  

% PART 2: CALCULATE LIGHT-LIMITED (RUBP-LIMITED) CO2 ASSIMILATION RATE
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Calculate J_max temperature correction factor
[f_T_J,f_TH_J] = f_T_corr(delHa_J,delHd_J,delS_J,R_g,T_l);

% Corrected max electron transport rate [micromoles electrons/m^2 LA/s)
J_max = J_max*f_T_J*f_TH_J;                                

% Electron transport rate from photosystem II [micromoles electrons/m^2 LA/s]
I_psii = 0.5*Sig_psii*4.6*Q; 

% Electron transport rate [micromoles/m^2 LA/s)
J = min(roots([Theta_psii,-(I_psii + J_max),I_psii*J_max]));   

% Calculate Gamma temperature correction factor
[f_T_G,~] = f_T_corr(delHa_G,0,0,R_g,T_l);                 

% Corrected CO2 compensation point [Pa]
Gamma = 42.75*10^(-6)*P_atm*f_T_G;

% Light-limited (RuBP-limited) CO2 assimilation rate [micromoles CO2/m^2 LA/s]
A_j = J*(c_i - Gamma)/(4*c_i + 8*Gamma);                   

% PART 3: CALCULATE RUBISCO-LIMITED CO2 ASSIMILATION RATE
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Calculate V_max temperature correction factor
[f_T_V,f_TH_V] = f_T_corr(delHa_V,delHd_V,delS_V,R_g,T_l);

% Corrected max Rubisco carboxylation rate [micromoles CO2/m^2 LA/s)
V_max = V_max*f_T_V*f_TH_V;    

% Calculate temperature correction factor for K_c
[f_T_Kc,~] = f_T_corr(delHa_Kc,0,0,R_g,T_l);

% Corrected Michaelis-Menten Constant for CO2 [Pa]
K_c = 404.9*10^(-6)*P_atm*f_T_Kc;  

% Calculate temperature correction for K_o
[f_T_Ko,~] = f_T_corr(delHa_Ko,0,0,R_g,T_l);

% Corrected Michaelis-Menten Constant for O2 [Pa]
K_o = 278.4*10^(-3)*P_atm*f_T_Ko;                                 

% Partial pressure of O2 [Pa]
O = 0.209*P_atm; 

% Rubisco limited CO2 assimilation rate [micromoles CO2/m^2 LA/s]
A_c = (V_max*(c_i - Gamma))/(c_i + K_c*(1+O/K_o));                

% PART 4: CALCULATE PRODUCT-LIMITED CO2 ASSIMILATION RATE
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Calculate temperature correction factor for T_p
[f_T_Tp,f_TH_Tp] = f_T_corr(delHa_Tp,delHd_Tp,delS_Tp,R_g,T_l); 

% Corrected Triose Phosphate utilization rate [micromoles CO2/m^2 LA/s]
T_p = 0.167*V_max*f_T_Tp*f_TH_Tp;  

% Product-limited assimilation rate [micromoles CO2/m^2 LA/s]
A_p = 3*T_p;                       

% PART 5: CALCULATE NET ASSIMILATION RATE USING COLIMITATION
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Empirical curvature parameters
Theta_cj = 0.98; 
Theta_ip = 0.95;

% Calculate co-limited CO2 assimilation rate [micromoles CO2/m^2 LA/s] - (2)
A_i = min(roots([Theta_cj,-(A_c + A_j),A_c*A_j]));
A = min(roots([Theta_ip,-(A_i + A_p),A_i*A_p]));     

% Calculate R_d temperature correction factor
[f_T_Rd,f_TH_Rd] = f_T_corr(delHa_Rd,delHd_Rd,delS_Rd,R_g,T_l);

% Corrected dark respiration rate [micromoles CO2/m^2 LA/s]
R_d = 0.015*V_max*f_T_Rd*f_TH_Rd;                          

% Net CO2 assimilation rate [micromoles CO2/m^2 LA/s]
A_n = A - R_d;
% Remove negative A_n values
A_n = max(A_n,0);
% Save which mechanism control assimilation
[~,ContA] = min([A_c,A_j,A_p]); 

end

% PART 6: TEMPERATURE CORRECTION FUNCTION
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

function [f_T,f_TH] = f_T_corr(delHa,delHd,delS,R_g,T_l)

% Temperature correction constants and functions based on Eq. 8.10-8.11 in
% Reference (2).
f_T = exp((delHa/(298.15*0.001*R_g))*(1-298.15/(T_l + 273.15))); 
f_TH = (1+exp((298.15*delS - delHd)/(298.15*0.001*R_g)))/...
    (1+exp((delS*(T_l + 273.15) - delHd)/(0.001*R_g*(T_l + 273.15))));  

end
