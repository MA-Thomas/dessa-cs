function [CDF_F_list,distances,numVariances,var_list] = preCompute_propensityCDF_curves(...
    contact_sphere_radiusSqd,T_HARDLIMIT_DIFFUSE,maxRxDist_T_HARD,Da,Db,...
    earliestPropensityTime)

% % COMPUTE INPUTS
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% % Grid of points.

% % FIXED PARAMETERS
% n_sigma = 5; % used to calc max diffusion distance corresponding to max diffusion time
% Dfast = 1; 
% contact_sphere_radiusSqd = .8;
% T_HARDLIMIT_DIFFUSE = 5; % units [s]
% Da = 1; Db = 1; % units [um^2 / s]
% maxRxDist_T_HARD = n_sigma*sqrt(6*Dfast*T_HARDLIMIT_DIFFUSE); % units [um]
% % FIXED PARAMETERS


% variance_min = n_sigma*sqrt(6*Da*1e-10) + n_sigma*sqrt(6*Db*1e-10);
% variance_max = n_sigma*sqrt(6*Da*T_HARD_Diffuse_Limit) + n_sigma*sqrt(6*Db*T_HARD_Diffuse_Limit);

numDistances = 700;%100;
numVariances = 1000;%600;%200;

d_min = contact_sphere_radiusSqd;
d_max = 2*maxRxDist_T_HARD; % each particle can diffuse a max distance <maxRxDist_T_HARD>. Thus max separation dist is twice that.

d_list = logspace( log10(d_min),  log10(d_max), numDistances )';
distances = d_list;
d_list = repmat(d_list,1,numVariances);
d = d_list'; 
d_list = d(:);

times_list = logspace( log10(earliestPropensityTime),  log10(T_HARDLIMIT_DIFFUSE), numVariances )';
var_list = (6*Da*times_list) + (6*Db*times_list);

var_list = repmat(var_list,1,numDistances);
var_list = var_list(:);
times_list = repmat(times_list,1,numDistances);
times_list = times_list(:);
gridPoints = [d_list,var_list];

% % COMPUTE OBJECTIVE VALUES OF INPUTS

% Given inputs, compute objective function values
% (obj value = Integrated Propensity = rate_const*CDF_F(t) = 
%                              rate_const* F(var) integrated from 
%                              var 0 to var(t)
%                              where F is Pr( Q(X) < R_contact^2 )

% rate_constant will be multiplied during the simulation (because the value
% of the rate constant depends on the particular assemblies considered).
[CDF_F_list] = computeObjFunction(...
    gridPoints,numVariances,contact_sphere_radiusSqd,times_list);
CDF_F_list(CDF_F_list < 1e-12) = 0;

end

