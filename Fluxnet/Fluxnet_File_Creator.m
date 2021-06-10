%--------------------------------------------------------------------------
% Script Name: Fluxnet_File_Creator.m
% By:  Brandon Sloan
% Date Created: 6/4/2018
% Description:  This code imports the variables of interest from the 166
% tier 1 fluxnet data sites and stores the data in tables.
%--------------------------------------------------------------------------
clc
clear
cd('C:\Users\sloan091\Documents\PhD\Current Research Projects\Fluxnet_Analysis\Data\Raw\Fluxnet2015')
folders = dir('FLX*');

for i = 141%154:length(folders)
    [dt,pname] = importFluxnet(folders(i).name);
    eval([pname,'=dt;'])
    save(['C:\Users\sloan091\Documents\PhD\Current Research Projects\Fluxnet_Analysis\Data\MATLAB\',pname],pname)
    clearvars -except folders i
end

%% This loop can be used if I ever want to store the data in structure form
% for i = 1:10
%     [dt,pname] = importFluxnet(folders(i).name);
%     Fluxnet2015.(pname) = dt;
% end
