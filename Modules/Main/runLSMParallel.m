% =========================================================================
% Name   : runLSMParallel.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This is the main function for running the land surface model from Sloan
% et. al (2021). This script takes in soil and plant parameters as well as
% subsurface and atmospheric boundary conditions and calculates the surface
% energy budget components of a dual source two big-leaf steady state land
% surface model similar to CLM v5. Section S2 of Sloan et al. (2021)
% contains all the model formulation details as well as parameter sets
% used to model the Ameriflux US-Me2 Ponderosa pine site. This version of
% the LSM runs in parallel using a parfor loop.
%
% INPUTS
%   Const    - Physical constants
%   Flag     - Simulation settings
%   Plant    - Plant-specific parameters
%   Soil     - Soil-specific parameters
%   FluxData - Environmental forcing data
%
% OUTPUTS
%   LSM_Results_  - Resulting LSM fluxes and states using the following
%                   transpiration downregulation schemes:
%
%       WW   - No transpiration downregulation (i.e., well-watered)
%       Beta - Empirical beta function
%       PHM  - Plant hydraulics model
% =========================================================================

function [LSM_Results_WW,LSM_Results_Beta,LSM_Results_PHM] = ...
    runLSMParallel(Const,Flag,Plant,Soil,FluxData)


% PART 1: INITIALIZE MODEL PARAMETERS, BOUNDARIES, AND SETTINGS
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Pre-allocate memory for output results and create error vector
LSM_Results_WW   = createOutputs(size(FluxData,1));
LSM_Results_PHM  = createOutputs(size(FluxData,1));
LSM_Results_Beta = createOutputs(size(FluxData,1));
Error            = repmat(-9999,1,size(LSM_Results_WW,2));

% Set optimization criteria used in model
opts.lsqnl         = optimset(@lsqnonlin);
opts.lsqnl.Display = 'none';
opts.lsqnl.TolX    = 1e-6;
opts.lsqnl.TolFun  = 1e-6;
opts.fz            = optimset(@fzero);
opts.fz.Display    = 'none';
opts.fz.TolX       = 1e-6;


% PART 2: RUN LSM FOR EACH TIMESTEP OF FORCING DATA
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

parfor step = 1:size(FluxData,1)
    
    % PART 3: RUN THE LSM UNDER WELL-WATERED CONDITIONS
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    % Create boundary conditions file for the ii
    BC = createBCs(FluxData(step,:),Soil,Const);
    
    % Initial conditions of the outer solver decision variables [degrees C]
    %              T_l_sl   |  T_l_sh  |  T_g
    OS_DVs_i = [1.1*BC.T_a,1.1*BC.T_a,1.1*BC.T_a];
    
    % Set boundaries for outer solver [degrees C]
    lb = [0,0,0];
    ub = [80,80,80];
    
    % Clears the OS_Solver persistant variables for each iteration
    clearSolver();
    
    % Define anonymous function to be used MATLAB's solver
    OS_Resid = @(OS_DVs)Outer_Solver_WW(OS_DVs,Const,BC,...
        Plant,Soil,opts,Flag);
    
    % Start the try-catch to ensure sim runs even if some steps fail
    try
        
        % Run solver and find optimum SEB temperatures
        [OS_DVs_opt,~,~,exitflag,~] = lsqnonlin(OS_Resid,OS_DVs_i,...
            lb,ub,opts.lsqnl);
        
        % Calculate results for optimum solution
        [~,Rad,StepResults] = Outer_Solver_WW(OS_DVs_opt,Const,BC,...
            Plant,Soil,opts,Flag);
        
        % Store well-watered LSM results
        LSM_step_WW = storeOutputs(Rad,StepResults,exitflag);
        LSM_Results_WW(step,:) = LSM_step_WW;
        WWerror = 0;
        
    catch ME
        
        % Displays error and sets time step outputs to -9999
        disp(ME.message);
        LSM_Results_WW{step,:} = Error;
        WWerror = 1;
        
    end
    
    % PART 4A: USE A PHM DOWNREGULATION SCHEME
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    if isequal(Flag.PHMflag,1)
        
        if isequal(WWerror,0)
            % If well-watered solution worked, proceed
            % Set boundaries for solver [degrees C]
            lb = [0,0,0,0,0];
            ub = [80,80,80,LSM_step_WW.g_s_sl*2,LSM_step_WW.g_s_sh*2];
            
            % Initial conditions of decision variables [degrees C | mol/m^2/s H2O]
            OS_DVs_i = [LSM_step_WW.T_l_sl,LSM_step_WW.T_l_sh,LSM_step_WW.T_g,...
                LSM_step_WW.g_s_sl,LSM_step_WW.g_s_sh];
            
            % Define function to be used by solver
            OS_Resid = @(OS_DVs)Outer_Solver_PHM(OS_DVs,LSM_step_WW,Const,BC,...
                Plant,Soil,opts,Flag);
            
            % Try-catch for LSM with PHM downregulation
            try
                
                % Run solver and find optimum SEB temperatures and sunlit
                % and shaded leaf stomatal conductances
                [OS_DVs_opt,~,~,exitflag,~] = lsqnonlin(OS_Resid,OS_DVs_i,...
                    lb,ub,opts.lsqnl);
                
                % Calculate results for optimum solution
                [~,Rad,StepResults] = Outer_Solver_PHM(OS_DVs_opt,...
                    LSM_step_WW,Const,BC,Plant,Soil,opts,Flag);
                
                
                % Store PHM downregulated LSM results
                LSM_step_PHM = storeOutputs(Rad,StepResults,exitflag);
                LSM_Results_PHM(step,:) = LSM_step_PHM;
                
            catch ME
                
                % Displays error and sets time step outputs to -9999
                disp(ME.message);
                LSM_Results_PHM{step,:} = Error;
                
            end
        else
            % Sets time step outputs to -9999 if error occurs.
            LSM_Results_PHM{step,:} = Error;
            
        end
        
    else
    end
    
    % PART 4B: USE A BETA DOWNREGULATION SCHEME
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    if isequal(Flag.Betaflag,1) || isequal(Flag.Betaflag,2)
        
        if isequal(WWerror,0)
            % If well-watered solution worked, proceed
            % Set boundaries for solver [degrees C]
            lb = [0,0,0,0,0];
            ub = [80,80,80,LSM_step_WW.g_s_sl*2,LSM_step_WW.g_s_sh*2];
            
            % Initial conditions of decision variables, which are the
            % well-watered states [degrees C | mol/m^2/s H2O]
            OS_DVs_i = [LSM_step_WW.T_l_sl,LSM_step_WW.T_l_sh,...
                LSM_step_WW.T_g,LSM_step_WW.g_s_sl,LSM_step_WW.g_s_sh];
            
            % Calculate beta factor using either static or dynamic beta
            if isequal(Flag.Betaflag,1)
                
                % Calculate sunlit and shaded static beta factors
                BC.beta_sl = calcBeta(BC.psi_s,Plant.psi_s_50_sl,...
                    Plant.b_s_sl);
                BC.beta_sh = calcBeta(BC.psi_s,Plant.psi_s_50_sh,...
                    Plant.b_s_sh);
                
            elseif isequal(Flag.Betaflag,2)
                
                % Convert sunlit and shaded well-watered transpiration rate
                % from W/m^2 to mm/day
                T_ww_sl = LSM_step_WW.E_sl*86400/Const.L_vap;
                T_ww_sh = LSM_step_WW.E_sh*86400/Const.L_vap;
                
                % Calculate sunlit and shaded dynamic beta factors
                BC.beta_sl = calcBetaDyn(BC.psi_s,T_ww_sl,...
                    Plant.psi_s_50_slope_sl,Plant.psi_s_50_int_sl,...
                    Plant.b_s_slope_sl,Plant.b_s_int_sl);
                BC.beta_sh = calcBetaDyn(BC.psi_s,T_ww_sh,...
                Plant.psi_s_50_slope_sh,Plant.psi_s_50_int_sh,...
                    Plant.b_s_slope_sh,Plant.b_s_int_sh);
                
            else
            end
            
            % Define function to be used by solver
            OS_Resid = @(OS_DVs)Outer_Solver_Beta(OS_DVs,LSM_step_WW,...
                Const,BC,Plant,Soil,Flag);
            
            try
                
                % Run solver and find optimum SEB temperatures
                [OS_DVs_opt,~,~,exitflag,~] = lsqnonlin(OS_Resid,OS_DVs_i,...
                    lb,ub,opts.lsqnl);
                
                % Calculate results for optimum solution
                [~,Rad,StepResults] = Outer_Solver_Beta(OS_DVs_opt,...
                    LSM_step_WW,Const,BC,Plant,Soil,Flag);
                
                % Store beta downregulated LSM results
                LSM_step_Beta = storeOutputs(Rad,StepResults,exitflag);
                LSM_Results_Beta(step,:) = LSM_step_Beta;
                
            catch ME
                % Displays error and sets time step outputs to -9999
                disp(ME.message);
                LSM_Results_Beta{step,:} = Error;
                
            end
            
        else
            % Sets time step outputs to -9999 if error occurs
            LSM_Results_Beta{step,:} = Error;
            
        end
        
    else
    end
    
end


% PART 5: POST-PROCESS RESULTS
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% For some reason MATLAB may randomly add an imaginary 0 components to
% results. This has no impact on the solution, but needs to be removed.
% The next three lines do this.
names = LSM_Results_WW.Properties.VariableNames;
LSM_Results_WW = varfun(@real,LSM_Results_WW);
LSM_Results_WW.Properties.VariableNames = names;

names = LSM_Results_PHM.Properties.VariableNames;
LSM_Results_PHM = varfun(@real,LSM_Results_PHM);
LSM_Results_PHM.Properties.VariableNames = names;

names = LSM_Results_Beta.Properties.VariableNames;
LSM_Results_Beta = varfun(@real,LSM_Results_Beta);
LSM_Results_Beta.Properties.VariableNames = names;


end



