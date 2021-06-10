% =========================================================================
% Script Name   : Flux_Select.m
% Author        : Brandon Sloan
% Start Date    : Aug 12, 2019
% Last Updated  : Aug 12, 2019
%
% Description   : This function selects the appropriate FLUXNET or
% Ameriflux data for use as forcings in the simulations.
%
%   INPUTS:
%   years      - A vector of the years to select
%   months     - A vector of the months to select from the dataset
%   days       - A vector of the days to select from the dataset
%   hours      - A vector of the hours to select from the dataset
%   SiteFile   - Name of the .mat file containing the flux tower data
%
%   OUTPUTS:
%   FluxData   - A MATLAB table of the selected flux tower data.
% =========================================================================

function [FluxData,datelog] = Flux_Select(years,months,days,hours,SiteFile)

% Load FLUXNET dataset
Data = load(SiteFile);
nam = fieldnames(Data);
Data = Data.(nam{1,1});

% Find and remove all datapoints that are 24 hours after a rainfall event
Plog = Data.P_F>0;
Pind = find(Data.P_F>0);
nn = 24;
twod = ones(nn,1);
for i = 1:length(Pind)
    Plog(Pind(i):Pind(i)+nn-1,1) = twod;
end
Data = Data(~Plog,:);    
    
% Select the months, hours and years I want to predict
dates = datevec(Data.TIMESTAMP_START);
logseason = ismember(dates(:,1),years) & ismember(dates(:,2),months) & ismember(dates(:,3),days) & ismember(dates(:,4),hours);
% Removes any timesteps that have insufficient measurements
lognan = isnan(Data.G_F_MDS) | isnan(Data.H_F_MDS) | isnan(Data.LE_F_MDS) | isnan(Data.TA_F)...
    | isnan(Data.RH) | isnan(Data.PA_F) | isnan(Data.USTAR) | isnan(Data.TS_F_MDS_1) | Data.SWC_1_5_1==-9999 | Data.NIGHT == 1 ...
    | isnan(Data.PPFD_IN) | isnan(Data.PPFD_DIF) | isnan(Data.LW_OUT) | isnan(Data.SW_OUT) | isnan(Data.WS_F)| Data.GPP_DT_VUT_REF < 0;
% Creates the finalized data product
FluxData = Data(logseason & ~lognan,:);
datelog = logseason & ~lognan;
end