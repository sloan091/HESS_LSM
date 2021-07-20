% =========================================================================
% Script Name   : CC.m
% Author        : Brandon Sloan
% Start Date    : Mar 18, 2019
% Last Updated  : Mar 18, 2019
%
% Description   : This function calculates the saturation vapor pressure
% and slope of the saturation vapor pressure-temperature curve using the
% Clausius-Clapeyron relationship.
%
%   INPUTS:
%   T - Temperature [degrees C]
%   
%
%   OUTPUTS:
%   e_sat - Saturated vapor pressure at T [Pa]
%   de_dT - Slope of Clausius Clapeyron [Pa/C]
% =========================================================================

function [e_sat,de_dT] = CC(T)
a = 611; b = 17.27; c = 273.3;
e_sat = a*exp(b*T./(T+c)); 
de_dT = ((a*b*c)/(T + c).^2).*exp(b*T./(T + c));
end
