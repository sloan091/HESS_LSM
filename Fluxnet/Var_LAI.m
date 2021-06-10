% =========================================================================
% Script Name   : Var_LAI.m
% Author        : Brandon Sloan
% Start Date    : Aug 12, 2019
% Last Updated  : Aug 12, 2019
%
% Description   : This function creates a time series of leaf area index
% interpolated from sparse measurements at the flux tower sites.
%
%   INPUTS:
%   FluxData - A MATLAB table of the selected flux tower data.
%
%   OUTPUTS:
%   Plant      - Updated plant structure with variable LAI.
% =========================================================================
function Plant = Var_LAI(LAI_file,FluxData,Plant)

% Load and interpolate
load(LAI_file)
LAI_ts = interp1(LAI(:,1),LAI(:,2),FluxData.TIMESTAMP_START,'linear','extrap');

% Replace any negative value with the previous value
if sum(LAI_ts<=0) >0
    ir = find(LAI_ts<=0);
    for nnn = 1:length(ir)
        LAI_ts(ir(nnn)) = LAI_ts(ir(nnn)-1);
    end
else
end
% Update the LAI value to a vector
Plant.LAI = LAI_ts;
end