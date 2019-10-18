function [wait_time,integPropValue,Var,Pk] = computeWaitTime_preComp_CDF(...
    IntegratedPropensity_curve, variance_list, Da, Db, t_offset_b)%,...
    %maxAllowed_t_wait)

% This file assumes there is no nonzero t_offset_a for particle A.

% Note:  IntegratedPropensity_curve = rate_constant * CDF_F;

% CDF_list size is [2,numVariances]
% The first row is CDF(var | d_lower) for the d_lower just below d in the
% CDF data set.
% The second row is CDF(var | d_upper) for the d_upper just above d in the
% CDF data set.

% times_list size is [1,numVariances]

% TODO: 
%{
0. Sample Pk
1. determine variance at which IntegatedPropensity(var_at_reaction|d) == Pk + IntegatedPropensity(var_at_t_offset|d)
2. determine the wait_time (depends on whether the t_offset is non zero)

%}

if t_offset_b > 0
    var_at_t_offset = 6*Db*t_offset_b;
    [~,ind_var_at_t_offset] = min(abs(variance_list - var_at_t_offset));
    IntegratedPropensity_at_var_at_t_offset = IntegratedPropensity_curve(ind_var_at_t_offset);
else
    IntegratedPropensity_at_var_at_t_offset = 0;
end

Pk = log( 1/rand );
% Pk = 0.3e-3;%0.1
inds = find(IntegratedPropensity_curve >= Pk + IntegratedPropensity_at_var_at_t_offset);



if ~isempty(inds)
    
    Var = variance_list(inds(1));
    integPropValue = IntegratedPropensity_curve( inds(1) );
    
    % Var == 6*Da*t + 6*Db*(t + t_offset) = t*(6Da + 6Db) + t_offset*(6Db)
    wait_time = (Var - t_offset_b*6*Db) / (6*Da + 6*Db);
    assert(wait_time>0)
else
    % No reaction in allowed max diffusion time
    wait_time = -1;
    integPropValue = -1;
    Var = -1;

end

% % Do not consider diffusion times greater than <maxAllowed_t_wait>
% % which is set by the distance to the nearest relevant boundary.
% if wait_time > maxAllowed_t_wait; wait_time = -1; end









end

