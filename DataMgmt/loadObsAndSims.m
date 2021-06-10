function [FluxData,Sims] = loadObsAndSims(ObsFile,SimFiles,M,D,H,Y)

% Load Ameriflux data
FluxData = unpacker(ObsFile);
year     = Y;
month    = M;
days     = D;
hours    = H;
dates    = datevec(FluxData.TIMESTAMP_START);
DateLog  = ismember(dates(:,1),year) & ismember(dates(:,2),month)...
    & ismember(dates(:,3),days) & ismember(dates(:,4),hours);

% Load simulations
for i = 1:length(SimFiles)
    Temp = unpacker(SimFiles{i});
    % Create logic vector to remove erroneous timesteps
    if isequal(i,1)
        ErrorLog = abs(Temp.EB_g) > 1 | abs(Temp.EB_sl) > 1 ...
            | abs(Temp.EB_sh) > 1;
    else
        ErrorLog = ErrorLog | abs(Temp.EB_g) > 1 | abs(Temp.EB_sl) > 1 ...
            | abs(Temp.EB_sh) > 1;
    end
    
    TempSim{i,1} = Temp;
end

if isequal(size(FluxData,1),size(TempSim{1},1))
    
    RemoveLog = DateLog & ~ErrorLog;
    % Only keep wanted dates and remove erroneous timesteps
    FluxData = FluxData(RemoveLog,:);
    % For any number of simulation files
    for i = 1:length(SimFiles)
        Sims{i} = TempSim{i,1}(RemoveLog,:);
    end
    
else
    % Only keep wanted dates and remove erroneous timesteps
    FluxData = FluxData(DateLog,:);
    % RemoveLog = DateLog & ~ErrorLog;
    FluxData = FluxData(~ErrorLog,:);
    % For any number of simulation files
    for i = 1:length(SimFiles)
        Sims{i} = TempSim{i,1}(~ErrorLog,:);
    end
end

end