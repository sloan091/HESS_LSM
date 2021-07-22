% =========================================================================
% Name   : LSM_Serial.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This script runs the LSM from Sloan et al. (2021) with each time step
% serially.
% =========================================================================

clc
clear


% Load observed flux data to force the LSM.
ObsFile = 'FLX_US_Me2_LSM_Forcing.mat';           
FluxData = unpacker(ObsFile);

% Load LSM parameter set.
load Parameters_HESS_PHM_and_beta_s.mat

% Run serial LSM
[LSM_Results_WW_check,LSM_Results_Beta,LSM_Results_PHM] = ...
    runLSMSerial(Const,Flag,Plant,Soil,FluxData);

% Save LSM results
save('LSM_HESS_Resuls_PHM_and_beta_s','LSM_Results_WW','LSM_Results_Beta',...
    'LSM_Results_PHM');

