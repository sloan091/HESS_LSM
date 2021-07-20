% =========================================================================
% Name   : Collatz_C4wTemp.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% Description   : This function calculates the net CO2 assimilation using
% the photosynthesis model by Collatz et al. (1992) for
% C4 plants. Here, we follow the formulation and assumptions laid out by
% Oleson et al. (2013) including temperature dependence. Area units are 
% specified by LA for per unit leaf area and GA for per unit ground area.
%
% INPUTS
%   V_max      - Maximum scaled Rubisco carboxylation rate 
%                [micromoles CO2 /m^2 LA/s]
%   Q          - Photosynthetically Active Radiation [W/m^2 LA]
%   P_atm      - Atmospheric pressure [Pa]
%   c_i        - CO2 partial pressure inside the leaf [Pa]
%   T_l        - Leaf temperature [C]
%
% OUTPUTS
%   A_n   - Net CO2 assimilation rate [micromoles CO2/m^2 LA/s]
%   A     - CO2 assimilation rate [micromoles CO2/m^2 LA/s]
%   R_d   - Dark respiration rate [micromoles CO2/m^2 LA/s]
%   ContA - Which process is limiting: 1 Rubisco, 2 light, and 3 production
%
% REFERENCES
%   (1) Collatz, G., Ribas-Carbo, M., & Berry, J. (1992). Coupled 
%   Photosynthesis-Stomatal Conductance Model for Leaves of C4 Plants. 
%   Functional Plant Biology. https://doi.org/10.1071/pp9920519
%
%   (2) Oleson, K. W. et al. (2013). Technical Description of version 4.5 
%   of the Community Land Model (CLM). NCAR Technical Notes. 
%   https://doi.org/10.5065/D6RR1W7M
%
% =========================================================================

function [A_n,A,R_d,ContA] = Collatz_C4wTemp(V_max,Q,P_atm,c_i,T_l)


% PART 1: CALCULATE LIGHT-LIMITED (RUBP-LIMITED) CO2 ASSIMILATION RATE
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Quantum efficiency [mol CO2/mol photon] - (2)
alpha_q = 0.05;        

% Light-limited (RuBP-limited) CO2 assimilation rate [micromoles CO2/m^2 LA/s]
A_j = alpha_q*4.6*Q;


% PART 2: CALCULATE RUBISCO-LIMITED CO2 ASSIMILATION RATE
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Temperature dependence parameter - (2)
Q_10 = 2;

% Calculate V_max temperature correction factor
[f_T,f_TH,f_TR] = f_T_corr(T_l);

% Corrected max Rubisco carboxylation rate [micromoles CO2 /m^2 LA/s]
V_max = V_max*(Q_10^((T_l + 273.15 - 298.15)/10))/(f_T*f_TH);

% Rubisco limited CO2 assimilation rate [micromoles CO2/m^2 LA/s]
A_c = V_max;


% PART 3: CALCULATE PRODUCT-LIMITED CO2 ASSIMILATION RATE
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Initial slope of the C4 CO2 response curve - (2)
K_p = V_max*20000;

% Calculate temperature correction for K_p
K_p = K_p*Q_10^((T_l + 273.15 - 298.15)/10); 

% Product-limited assimilation rate [micromoles CO2/m^2 LA/s]
A_p = K_p*c_i/P_atm;                    


% PART 4: CALCULATE NET ASSIMILATION RATE USING COLIMITATION
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Empirical curvature parameters - (2)
Theta_cj = 0.8; 
Theta_ip = 0.95;

% Calculate co-limited CO2 assimilation rate [micromoles CO2/m^2 LA/s] - (2)
A_i = min(roots([Theta_cj,-(A_c + A_j),A_c*A_j]));
A = min(roots([Theta_ip,-(A_i + A_p),A_i*A_p]));     

% Dark respiration rate corrected for temperature [micromoles CO2/m^2 LA/s]
R_d = 0.015*V_max*(Q_10^((T_l + 273.15 - 298.15)/10))/f_TR;                         

% Net CO2 assimilation rate [micromoles CO2/m^2 LA/s]
A_n = A - R_d;
% Remove negative A_n values
A_n = max(A_n,0);
% Save which mechanism control assimilation
[~,ContA] = min([A_c,A_j,A_p]); 

end

% PART 6: TEMPERATURE CORRECTION FUNCTION
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

function [f_T,f_TH,f_TR] = f_T_corr(T_l)

% Temperature correction constants and functions based on Eq. 8.12-8.13 in
% Reference (2).
s1 = 0.3; s2 = 313.15; s3 = 0.2;  s4 = 288.15; s5 = 1.3; s6 = 328.15;
f_T = 1 + exp(s1*(T_l + 273.15 - s2));   
f_TH = 1 + exp(s3*(s4 - T_l - 273.15));  
f_TR = 1 + exp(s5*(T_l + 273.15 - s6)); 
end

