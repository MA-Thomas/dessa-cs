function [CDF_F] = computeObjFunction(...
    gridPoints,numVariances,contact_sphere_radiusSqd,t)

% This function is based on computeWaitTime______v2_alt.m

numInputs = size(gridPoints,1);
F = zeros(numInputs,1);
CDF_F = zeros(numInputs,1);
dist_counter = 1;


Dist = gridPoints(1,1);
[x,y,z] = sph2cart(2*pi*rand,pi*rand,Dist);
mu_colVec = [x,y,z]';
Var = gridPoints(1,2);
[F(1),convergence] = Pr_dist2Gaussians_powerSeries_alt(mu_colVec,Var,...
contact_sphere_radiusSqd);
CDF_F(1) = 0;

% Sometimes the molecules are just too far apart (variance is too low given
% the available propensity integration max duration). This can lead to
% convergence problems. When these occur, just set F=0;
if convergence~=1
    F(1) = 0;
end
% assert(convergence==1)


for i = 2:numInputs
    
    % For the current distance, loop over all propensity integration time
    % upper limits.
    % I.e. loop over all variances at the current distance.
    
    Var = gridPoints(i,2);    
    Dist = gridPoints(i,1);   
    
    [x,y,z] = sph2cart(2*pi*rand,pi*rand,Dist);
    mu_colVec = [x,y,z]';

    [F(i),convergence] = Pr_dist2Gaussians_powerSeries_alt(mu_colVec,Var,...
    contact_sphere_radiusSqd);
    if convergence~=1
        F(1) = 0;
    end
    % assert(convergence==1)   
    
    dist_counter = dist_counter + 1;       

    if dist_counter == numVariances +1
       % We've finished integrating the probability(Var) at the
       % current pairwise reactant distance.  
       % Start integrating from Var=0 at next pairwise reactant distance.
       CDF_F(i) = 0;
       dist_counter = 1;
    else        
        CDF_F(i) = CDF_F(i-1) + trapz( [t(i-1),t(i)],[F(i-1),F(i)] );
    end   


end

end