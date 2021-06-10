% =========================================================================
% Name   : FvCB_C3.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% Description   : This function calculates the net CO2 assimilation using
% the photosynthesis model by Farquhar, von Caemmerer and Berry (1980) for
% C3 plants. Here, we follow the formulation and assumptions laid out by
% Oleson et al. (2013) excluding temperature dependence as done in Sloan 
% et al. (2021). See Sect. S2.4.2 for full details. Area units are 
% specified by LA for per unit leaf area and GA for per unit ground area.
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

function [A_n,A,R_d,ContA] = FvCB_C3(V_max,J_max,Sig_psii,Q,Theta_psii,P_atm,c_i)


% PART 1: CALCULATE LIGHT-LIMITED (RUBP-LIMITED) CO2 ASSIMILATION RATE
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''                             

% Electron transport rate from photosystem II [micromoles electrons/m^2 LA/s]
I_psii = 0.5*Sig_psii*4.6*Q; 

% Electron transport rate [micromoles/m^2 LA/s)
J = min(roots([Theta_psii,-(I_psii + J_max),I_psii*J_max]));                 

% CO2 compensation point [Pa]
Gamma = 42.75*10^(-6)*P_atm;

% Light-limited (RuBP-limited) CO2 assimilation rate [micromoles CO2/m^2 LA/s]
A_j = J*(c_i - Gamma)/(4*c_i + 8*Gamma);  


% PART 2: CALCULATE RUBISCO-LIMITED CO2 ASSIMILATION RATE
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Michaelis-Menten Constant for CO2 [Pa]
K_c = 404.9*10^(-6)*P_atm;  

% Michaelis-Menten Constant for O2 [Pa]
K_o = 278.4*10^(-3)*P_atm;                                 

% Partial pressure of O2 [Pa]
O = 0.209*P_atm; 

% Rubisco limited CO2 assimilation rate [micromoles CO2/m^2 LA/s]
A_c = (V_max*(c_i - Gamma))/(c_i + K_c*(1+O/K_o));  


% PART 3: CALCULATE PRODUCT-LIMITED CO2 ASSIMILATION RATE
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Triose Phosphate utilization rate [micromoles CO2/m^2 LA/s]
T_p = 0.167*V_max;  

% Product-limited assimilation rate [micromoles CO2/m^2 LA/s]
A_p = 3*T_p;  


% PART 4: CALCULATE NET ASSIMILATION RATE USING COLIMITATION
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Empirical curvature parameters
Theta_cj = 0.98; 
Theta_ip = 0.95;

% Calculate co-limited CO2 assimilation rate [micromoles CO2/m^2 LA/s] - (2)
A_i = min(roots([Theta_cj,-(A_c + A_j),A_c*A_j]));
A = min(roots([Theta_ip,-(A_i + A_p),A_i*A_p]));     

% Dark respiration rate [micromoles CO2/m^2 LA/s]
R_d = 0.015*V_max;                          

% Net CO2 assimilation rate [micromoles CO2/m^2 LA/s]
A_n = A - R_d;
% Remove negative A_n values
A_n = max(A_n,0);
% Save which mechanism control assimilation
[~,ContA] = min([A_c,A_j,A_p]); 

end

