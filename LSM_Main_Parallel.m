% Run bad LSM simulation to debug PHM
clc
clear

% Add the necessary paths and name run
addpath(genpath('/panfs/roc/groups/10/feng/sloan091/HESS_Codes'))

% Generate the parallel pool
workers = 24; parpool(workers)
if  isempty(gcp('nocreate'))
     parpool(workers)
     if isempty(gcp('nocreate'))
         error('Could not get all workers')
     else
     end
else
end


% Observation Files
ObsFile = 'FLX_US_Me2_LSM_Forcing.mat';          
FluxData = unpacker(ObsFile);

% Date constraints
year  = 2013:2014;
month = 5:8;
days  = 1:31;
hours = 8:20;

dates    = datevec(FluxData.TIMESTAMP_START);
DateLog  = ismember(dates(:,1),year) & ismember(dates(:,2),month)...
    & ismember(dates(:,3),days) & ismember(dates(:,4),hours);
FluxData = FluxData(DateLog,:);

% Parameter sets
load Parameters_HESS_Final.mat

% Loading parameters updates
Flag.Betaflag = 0;

[LSM_Results_WW,LSM_Results_Beta,LSM_Results_PHM] = ...
    runLSMParallel(Const,Flag,Plant,Soil,FluxData);

% Save outputs
save('LSM_HESS_Final','LSM_Results_WW','LSM_Results_Beta','LSM_Results_PHM');

delete(gcp);
