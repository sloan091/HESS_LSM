% =========================================================================
% Name   : calcBetaDyn.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This function creates the dynamic beta downregulation factor using a
% Weibull-type function as in Sloan et al. (2021). The dynamic beta couples
% the effects of soil water stress and atmospheric moisture demand on
% transpiration downregulation by adding an additional dependence on
% well-watered transpiration rate to the original static beta. See Sloan et
% al. (2021) for further discussion.
%
% INPUTS
%   psi_s          - Soil water potential [MPa]
%   T_ww           - Well-watered transpiration rate [mm/day]
%   psi_s_50_slope - Slope parameter for linear dependence of psi_s_50 on
%                    T_ww [MPa*day/mm].
%   psi_s_50_int   - Intercept parameter for linear dependence of psi_s_50 on
%                    T_ww [MPa].
%   b_s_slope      - Slope parameter for linear dependence of b_s on 
%                    T_ww [day/mm].
%   b_s_int        - Intercept parameter for linear dependence of b_s on 
%                    T_ww [day/mm].
%
% OUTPUTS
%   beta - Empirical beta downregulation factor
% =========================================================================

function beta = calcBetaDyn(psi_s,T_ww,psi_s_50_slope,psi_s_50_int,...
                            b_s_slope,b_s_int)

% Linear dependence of dynamic beta parameters on the well-watered
% transpiration rate. See Sect. S6.2 of Sloan et al. (2021) for details.

% Soil water potential at 50% loss of stomatal conductance [MPa]
psi_s_50 = psi_s_50_slope*T_ww + psi_s_50_int;

% Stomatal sensitivity parmaeter
b_s = b_s_slope*T_ww + b_s_int;

% Weibull formulation of empirical beta factor shown in Eq. 16 of Sloan et
% al. (2021) with the dependence on well-watered transpiration built into
% psi_s_50 and b_s.
beta = 2^(-(psi_s/psi_s_50)^b_s);

% Make sure beta is bounded between 0 and 1
beta = min(max(beta,0),1);

end