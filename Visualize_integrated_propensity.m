% Visualize integrated propensity

% Setting this too low can lead to non-convergence in
% Pr_dist2Gaussians_powerSeries_alt.m
earliestPropensityTime = 1e-4;

n_sigma = 5; % used to calc max diffusion distance corresponding to max diffusion time
Dfast = 1; 
contact_sphere_radiusSqd = .2;%3; % units [um]

% SETS THE UPPER LIMIT ON HOW LONG AN ASSEMBLY CAN DIFFUE BEFORE ITS NEXT
% EVENT. (Other factors like proximity to a boundary, or only considering "nearby"
% partners can lower this upper limit).
T_HARDLIMIT_DIFFUSE = 16; % units [s]
Da = 1; Db = 1; % units [um^2 / s]
maxRxDist_T_HARD = n_sigma*sqrt(6*Dfast*T_HARDLIMIT_DIFFUSE); % units [um]
% FIXED PARAMETERS


[CDF_F_list,distances,numVariances,var_list] = preCompute_propensityCDF_curves(...
    contact_sphere_radiusSqd,T_HARDLIMIT_DIFFUSE,maxRxDist_T_HARD,Da,Db,...
    earliestPropensityTime);

%%
d = 8;

[~,d_index] = min( abs(d - distances) );
begin_ind = (d_index-1)*numVariances+1;
end_ind = (d_index-1)*numVariances+numVariances;
IntegratedPropensity_curve = rate_constant * CDF_F_list(begin_ind:end_ind);
variances = var_list(begin_ind:end_ind);

t_offset_b = 0;
[wait_time,integPropValue,var_at_WaitTime,Pk] = computeWaitTime_preComp_CDF(IntegratedPropensity_curve, variances,...
    Da, Db, t_offset_b);

figure
plot(variances,IntegratedPropensity_curve)

hold on
if var_at_WaitTime>0; scatter(var_at_WaitTime, integPropValue);end
Pk_vector = Pk*ones(1,length(variances));
plot(variances,Pk_vector,'--')
text(10,Pk_vector(1)+0.005, ['\downarrow exponentially distributed R.N. (=',num2str(Pk),')'])
ylabel('Integrated Reaction Propensity','fontsize',16);
xlabel('Variance [um^2]','fontsize',16);
title('Reactant Pair Separation: d=8um','fontsize',16);
ylim([0 5e-4])
ax = gca;
ax.FontSize = 16; 
hold off

%% Plot Integrated Propensity at multiple distance values, d
figure
indices = [240:10:400, 408:12:500];
for d_index = indices%1:10:length(distances)

    begin_ind = (d_index-1)*numVariances+1;
    end_ind = (d_index-1)*numVariances + numVariances;
    IntegratedPropensity_curve = rate_constant * CDF_F_list(begin_ind:end_ind);
    variances = var_list(begin_ind:end_ind);
   
    plot(variances,IntegratedPropensity_curve)
    hold on
    
end

toString_distances = distances(indices);
legend( string( round(toString_distances,2) ), 'fontsize',16 )
ax = gca;
ax.FontSize = 16; 
ylabel('Integrated Propensity','fontsize',16)
xlabel('Variance (um^2)','fontsize',16)
title('Integrated Propensities at Increasing Reactant Pair Separation',...
    'fontsize',18)
hold off