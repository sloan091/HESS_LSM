% =========================================================================
% Name   : LSM_Parallel.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This script runs the LSM from Sloan et al. (2021) with each time step
% in parallel mode.  This is possible because we assume steady state at
% each time step and each can be calculated independently.
% =========================================================================

clc
clear

% Add the necessary paths and name run
addpath(genpath('/panfs/roc/groups/10/feng/sloan091/HESS_Codes'))

% Generate the parallel pool. Sometimes the parallel pool will not start,
% so this attempts to run it again if it fails.
workers = 24; parpool(workers)
if  isempty(gcp('nocreate'))
     parpool(workers)
     if isempty(gcp('nocreate'))
         error('Could not get all workers')
     else
     end
else
end


% Load observed flux data to force the LSM.
ObsFile = 'FLX_US_Me2_LSM_Forcing.mat';          
FluxData = unpacker(ObsFile);

% Load LSM parameter set.
load Parameters_HESS_PHM_and_beta_s.mat

% Run serial LSM
[LSM_Results_WW,LSM_Results_Beta,LSM_Results_PHM] = ...
    runLSMParallel(Const,Flag,Plant,Soil,FluxData);

% Save LSM results
save('LSM_HESS_Resuls_PHM_and_beta_s','LSM_Results_WW','LSM_Results_Beta',...
    'LSM_Results_PHM');

% Delete parallel pool and release workers
delete(gcp);
