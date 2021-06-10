% =========================================================================
% Name   : Scale_V_J.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This function scales the maximum carboxylation capacity
% (V_max) and electron tranport rate (J_max) based on leaf nitrogen
% content and the exponential light extinction using the methods laid
% out in De Pury and Farquhar (1997) and Dai et al. (2004). See Sect. S2.4.3 of
% Sloan et al. (2021) for full details. Note, I opted not to use Dai et al.
% (2004) for V_max because there appears to be a mistake in their Eq. 37a
% and 37b, where only k_b should be multiplied by LAI as in Eq. 22 of De
% Pury and Farquhar (1997). However, De Pury and Farquhar (1997) do not
% provide scaling for J_max, therefore, the formulations in Dai et al.
% (2004) were used for J_max and appear correct.
%
% INPUTS
%   V_max - Maximum Rubisco carboxylation rate assumed to be at top of 
%           canopy for scaling purposes [micromoles CO2 /m^2 leaf area/s] 
%   J_max - Maximum electron transport rate assumed to be at top of 
%           canopy for scaling purposes [micromoles electrons /m^2 leaf area/s] 
%   k_n    - Exinction coefficient of leaf nitrogen content through the 
%            through the canopy [-]
%   k_b    - Extinction rate of direct beam radiation through the 
%            through the canopy [-]
%   k_d    - Extinction rate of scattered diffuse radiation through the 
%            through the canopy [-]
%   LAI    - Leaf area index [m^2 leaf area/m^2 ground area]
%
% OUTPUTS
%   V_max_sl - Scaled maximum Rubisco carboxylation rate for sunlit big
%              leaf [micromoles CO2/m^2 ground area/s]
%   V_max_sh - Scaled maximum Rubisco carboxylation rate for shaded big
%              leaf [micromoles CO2/m^2 ground area/s]
%   J_max_sl - Scaled maximum electron transport rate for sunlit big
%              leaf [micromoles electrons/m^2 ground area/s]
%   J_max_sh - Scaled maximum electron transport rate for shaded big
%              leaf [micromoles electrons/m^2 ground area/s]
%
% REFERENCES
%   (1) De Pury, D. G. G., & Farquhar, G. D. (1997). Simple scaling of 
%   photosynthesis from leaves to canopies without the errors of big-leaf 
%   models. Plant, Cell and Environment, 20, 537–557. 
%   https://doi.org/10.1111/j.1365-3040.1997.00094.x
%
%   (2) Dai, Y. et al. (2004). A Two-Big-Leaf Model for Canopy Temperature,
%   Photosynthesis, and Stomatal Conductance. Journal of Climate, 17(12), 
%   2281–2299. https://doi.org/10.1175/1520-0442(2004)017<2281:ATMFCT>2.0.CO;2
%
% =========================================================================

function [V_max_sl,V_max_sh,J_max_sl,J_max_sh] = Scale_V_J(V_max,J_max,k_n,k_b,k_d,LAI)

% Calculate sunlit and shaded V_max scaled for exponential leaf nitrogen
% and light profiles [micromoles CO2/m^2 ground area/s] - (1)
V_max_sl = LAI*V_max*(1 - exp(-k_n - k_b*LAI))/(k_n + k_b*LAI);
V_max_sh = LAI*V_max*(1 - exp(-k_n))/k_n - V_max_sl; 

% Calculate sunlit and shaded J_max scaled for exponential light profiles
% [micromoles electrons/m^2 ground area/s] - (2)
J_max_sl = J_max*(1 - exp(-(k_d + k_b)*LAI))*1/(k_d + k_b);
J_max_sh = J_max*(1 - exp(-k_d*LAI))/k_d - J_max_sl; 

end


