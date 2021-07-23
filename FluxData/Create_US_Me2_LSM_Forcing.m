% =========================================================================
% Name   : Create_US_Me2_LSM_Forcing.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This script selects a subset of the combined FLUXNET2015 (1) and 
% (2) Ameriflux datasets for US-Me2 used as the LSM forcing in Sloan 
% et al. (2021). The data source references for US-Me2 are listed below.  
% Sloan et al. (2021) focused on May-August 2013-2014 as this time period 
% had all the requisite measurements to force the LSM.  Please see Sect. 
% S5 of Sloan et al. (2021) for full details.
%
% REFERENCES
%   (1) Pastorello, G. et al. (2020). The FLUXNET2015 dataset and the
%   ONEFlux processing pipeline for eddy covariance data. 
%   Scientific Data, 7(1), 225. https://doi.org/10.1038/s41597-020-0534-3
%   (2) Law, B. E. (2021). AmeriFlux US-Me2 Metolius mature ponderosa pine.
%   Ver. 16-5. (Dataset). Ameriflux AMP.
%   https://doi.org/doi.org/10.17190/AMF/1246076
% =========================================================================

clc
clear

% Unpack the combined FLUXNET2015/Ameriflux US-Me2 dataset
Data = unpacker('FLX_US_Me2.mat');

% Specify dates of interest
years  = 2013:2014;
months = 5:8;
days   = 1:31;
hours  = 8:20;

% Remove flux measurements 12 hours after a rainfall event
% Find rainfall events
Plog = Data.P_F > 0;
Pind = find(Data.P_F>0);

% Set how many timesteps to cancel after a rain event 
nn = 24;
rainday = ones(nn,1);
for i = 1:length(Pind)
    Plog(Pind(i):Pind(i)+nn-1,1) = rainday;
end
% Remove these post rain timesteps
Data = Data(~Plog,:);    
    
% Select the months, hours and years I want to predict
dates = datevec(Data.TIMESTAMP_START);
logdates = ismember(dates(:,1),years) & ismember(dates(:,2),months)...
    & ismember(dates(:,3),days) & ismember(dates(:,4),hours);

% Removes any timesteps that have missing or undesired measurements.
lognan = isnan(Data.G_F_MDS) | isnan(Data.H_F_MDS) | isnan(Data.LE_F_MDS)...
    | isnan(Data.TA_F) | isnan(Data.RH) | isnan(Data.PA_F) |...
    isnan(Data.USTAR) | isnan(Data.TS_F_MDS_1) | Data.SWC_1_5_1==-9999 |...
    Data.NIGHT == 1 | isnan(Data.PPFD_IN) | isnan(Data.PPFD_DIF) |...
    isnan(Data.LW_OUT) | isnan(Data.SW_OUT) | isnan(Data.WS_F) |...
    Data.GPP_DT_VUT_REF < 0;

% Creates the finalized data product and save
FluxData = Data(logdates & ~lognan,:);
datelog = logdates & ~lognan;

save('FLX_US_Me2_LSM_Forcing','FluxData')
