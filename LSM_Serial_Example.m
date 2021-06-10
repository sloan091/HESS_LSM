% Run bad LSM simulation to debug PHM
clc
clear

% Add the necessary paths and name run
addpath(genpath('/panfs/roc/groups/10/feng/sloan091/My_LSM'))
run_name = 'Me2_LM_wF';


% Observation Files
ObsFile = 'FLX_US_Me2_LSM_Forcing.mat';       
%ObsFile = 'FLX_US_Me2.mat';    
FluxData = unpacker(ObsFile);

% Date constraints
year  = 2014;
month = 8;
days  = 10:12;
hours = 8:20;

dates    = datevec(FluxData.TIMESTAMP_START);
DateLog  = ismember(dates(:,1),year) & ismember(dates(:,2),month)...
    & ismember(dates(:,3),days) & ismember(dates(:,4),hours);
FluxData = FluxData(DateLog,:);

% Parameter sets
load Parameters_HESS_Final.mat

% % New HESS parameters
% Soil.b      = 3.88;
% Plant.psi_l_50 = -2.5;
% Plant.bl = 4.36;
% Plant.psi_x_50 = -2;
% Plant.a =  1.85;
% Plant.k_sap = 0.6*10^(-7)*Plant.h_v*Const.rho_w;
% Soil.psi_sat = -10^(-1.9);


% Loading parameters updates
% Flag.PHMflag = 1;
% Flag.DownRegMethod = 1;
% Flag.DynBetaFlag = 0;
% Flag.SMflag = 1;
Flag.Betaflag = 0;

[LSM_Results_WW_check,LSM_Results_Beta,LSM_Results_PHM] = ...
    runLSMSerial(Const,Flag,Plant,Soil,FluxData);

% % Save outputs
% save('LSM_HESS_check2','LSM_Results_WW','LSM_Results_Beta','LSM_Results_PHM');
% 
% delete(gcp);
