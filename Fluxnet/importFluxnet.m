function [dt,pname] = importFluxnet(foldername)
%---------- Function Name:importFluxnet.m--------------------------------%
% By: Brandon Sloan
% Created: 6/4/2018
% Description: This file reads in the half hourly fluxnet data and keeps
% only the variables I am interested in for the Max Entropy analysis.  For
% ease of use, the Matlab table formate was utilized so headers could be
% kept.
%------------Input/output--------------------------------------------------%
% foldername = folder where fluxnet .csv files are located
% dt = Fluxnet data table containing all relevant variables stored in
% MATLAB table format
% pname = this is the 10 digit unique site identifier name used by Fluxnet
%-------------------------------------------------------------------------%

% Begin function
fname = dir(['.\',foldername,'\*_HH_*']);
path = ['.\',foldername,'\',fname.name];
pname = fname.name(1:10);

% These variables are not relevant to our current study and are removed
vars2del = {''};
% vars2del = {'SW_IN_POT','SW_IN_F','SW_IN_F_QC','LW_IN_F','LW_IN_F_QC','WD'...
%     'PPFD_IN','CO2_F_MDS','CO2_F_MDS_QC','LE_CORR_25','LE_CORR_75','H_CORR_25','H_CORR_75'...
%     'NEE_VUT_REF','NEE_VUT_REF_QC','NEE_VUT_REF_RANDUNC','NEE_VUT_25','NEE_VUT_50','NEE_VUT_75'...
%     'NEE_VUT_25_QC','NEE_VUT_50_QC','NEE_VUT_75_QC','RECO_NT_VUT_REF','RECO_NT_VUT_25','RECO_NT_VUT_50'...
%     'RECO_NT_VUT_75','GPP_NT_VUT_REF','GPP_NT_VUT_25','GPP_NT_VUT_50','GPP_NT_VUT_75','RECO_DT_VUT_REF'...
%     'RECO_DT_VUT_25','RECO_DT_VUT_50','RECO_DT_VUT_75','GPP_DT_VUT_REF','GPP_DT_VUT_25','GPP_DT_VUT_50','GPP_DT_VUT_75','RECO_SR','RECO_SR_N'};
opts = detectImportOptions(path);
opts = setvaropts(opts,'TreatAsMissing','-9999');
tempv = opts.VariableNames;
Lia = ismember(tempv,vars2del);
opts.SelectedVariableNames = tempv(~Lia);
dt = readtable(path,opts);

% This portion uploads the PATCH QC flags discussed here: 
% http://fluxnet.fluxdata.org/data/fluxnet2015-dataset/known-issues/
patchfile = dir(['C:\Users\sloan091\Documents\PhD\Current Research Projects\Fluxnet_Analysis\Data\Raw\FluxnetPatch\',pname,'*']);
path2 = ['C:\Users\sloan091\Documents\PhD\Current Research Projects\Fluxnet_Analysis\Data\Raw\FluxnetPatch\',patchfile.name];
patchtable = readtable(path2);
tempd = patchtable.TIMESTAMP_START;
Lia = ismember(tempd,dt.TIMESTAMP_START);
dt.H_F_MDS_QC = patchtable.H_F_MDS_QC(Lia);
dt.LE_F_MDS_QC = patchtable.LE_F_MDS_QC(Lia);

% This portion converts everything to MATLAB datenumbers
dates = num2str(dt.TIMESTAMP_START);
dt.TIMESTAMP_START = datenum(dates,'yyyymmddHHMM');
dates = num2str(dt.TIMESTAMP_END);
dt.TIMESTAMP_END = datenum(dates,'yyyymmddHHMM');

% Need to replace hyphen in case I want to store as a structure with this
% name
pname(7)='_';
end

