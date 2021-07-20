% =========================================================================
% Name   : calcBeta.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This function creates the emprical beta downregulation factor using a
% Weibull-type function as in Sloan et al. (2021). This code corresponds to
% the beta_s and beta_2L publications in said publication.
%
% INPUTS
%   psi_s - Soil water potential [MPa]
%   psi_s_50 - Soil water potential at 50% stomatal conductance loss[MPa]
%   b_s      - Stomatal sensitivity parameter
%
% OUTPUTS
%   beta - Empirical beta downregulation factor
% =========================================================================

function beta = calcBeta(psi_s,psi_s_50,b_s)

% Weibull formulation of empirical beta factor shown in Eq. 16 of Sloan et
% al. (2021) without the dependence on well-watered transpiration.  That
% dependence is reserved for the dynamic beta.
beta = 2^(-(psi_s/psi_s_50)^b_s);

% Make sure beta is bounded between 0 and 1
beta = min(max(beta,0),1);

end