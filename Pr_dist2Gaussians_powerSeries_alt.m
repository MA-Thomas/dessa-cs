function [F,convergence,mP] = Pr_dist2Gaussians_powerSeries_alt(mu_colVec,Var,...
    contact_sphere_radiusSqd)

% Output: F - Probability sample from N(mu_colVec,E) has magnitude less than contact_sphere_radius
% F is Pr(Q(X)<R) and thus 0 <= F <= 1.

% Inputs (3 dimensional space):
% vector of variances along each dimension: E = [sigma^2_NEW(1),sigma^2_NEW(2),sigma^2_NEW(3)]
% (Cov matrix: E = sigma^2_NEW * IdentityMatrix)

% contact_sphere_radiusSqd: the distance below which two assemblies are
% considered to be "in contact". (distance squared, that is)

%{
% THIS MAY NOT BE THE BEST WAY TO THINK ABOUT THE PROBABILITIES.
% EVEN THOUGH PHYSICALLY THERE MAY BE A MINIMAL ALLOWED DISTANCE BETWEEN
% DISTINCT REACTANTS, THIS DOES NOT IMPLY I NEED TO ACCOUNT FOR THIS BY
% SUBTRACTING OFF THE PROBABILITY MASS ASSOCIATED WITH dist<r_excluded.
% THAT WOULD NECESSITATE A RENORMALIZATION OF THE PROBABILITIES
% dist>=r_excluded. INSTEAD, I CAN JUST ASSUME THE PROBABILITY dist<r_contact
% IS CONCENTRATED IN A BAND BETWEEN r_excluded and r_contact.

% rSqd_excluded: the smallest allowed distance between assemblies such that
% they are not overlapping. (distance squared, that is)

% contact_sphere_radiusSqd > rSqd_excluded
%}
%{

Gaussian molecule A: a ~ N(mu_A, sigma^2_A)
Gaussian molecule B: b ~ N(mu_B, sigma^2_B)

Gaussian for difference vector: 
z = a-b ~ N( [mu_A - mu_B], [sigma^2_A + sigma^2_B] )
  = N(mu_NEW, sigma^2_NEW)

We're interested in the distribution of squared distances, i.e. the 
quadratic form:
 Q(z) = z(1)^2 + z(2)^2 + z(3)^2
 
More specifically, we're interested in the CDF(y) of Q(z), Pr{ Q(z) < y }

Also, the variances are time dependent, AND, there may be time
offsets.
e.g.
sigma^2_NEW(t) = sigma^2_A(t) + sigma^2_B(t+offset)
OR
sigma^2_NEW(t) = sigma^2_A(t+offset) + sigma^2_B(t)

%}


% F( lambda;b;y) Theorem 4.2b.1 Mathai-Provost p.95

% To investigate the frequencies of each case.
% Ensure that we first run the following before running the simulation:
% chunk1=0;chunk2=0;chunk3=0;chunk4=0;chunk5to20=0;chunk21plus=0;
% case1 = 0; case2 = 0; case3 = 0; case4 = 0; case5 = 0; save('pr_dist2Gaussians_powerSeries__CASES_STORAGE.mat','case1','case2','case3','case4','case5','chunk1','chunk2','chunk3','chunk4','chunk5to20','chunk21plus') 
% % load('pr_dist2Gaussians_powerSeries__CASES_STORAGE.mat','case1','case2','case3','case4','case5',...
% %     'chunk1','chunk2','chunk3','chunk4','chunk5to20','chunk21plus')    

% Propensity integration cannot start from exactly t==0 and t_offset==0.
assert(Var ~= 0); % This would lead to inf/NaN results.

E = Var*[1 1 1];
p = 3; % 3d space
% contact_sphere_radiusSqd = 4; % the distance squared below which 2 molecules react
y = contact_sphere_radiusSqd; % F is CDF(y)
% % % rE = rSqd_excluded; % F is CDF(rE)

% INTERESTING NOTE: as y increases, so must the number of terms included in
% the power series to maintain an accurate approximation
% See variable: <maxPowerSeries>.

%{
% Steps:
0. compute b_vector
1. compute d(k) for k = 1:5
2. compute c(k-1)
 
Follow steps on p.96

A = Identity, eye(3)
[P,D] = eig(CovMatrix)

Since the covariance matrix is diagonal, 
P <- 3x3 Identity
D <- Diagonal matrix of eigenVals. All eigenVals equal to sigma^2_NEW
i.e. THE EIGENVALUES ARE JUST THE VARIANCES ALONG EACH DIMENSION AND IN OUR
ISOMETRIC CASE, ALL VARIANCES ARE EQUAL

%}
% Write all the formulas using E (variances).
% In the final formula, replace E with time dependent variance.
% B = chol( E.*eye(3) ); % Matlab v2018
% B = chol( E(1).*eye(3) ); % Matlab v2016 temp
B = ( 1/sqrt(E(1)) ) * eye(3); % Since E = Var*[1 1 1], we're taking the chol
%                        decomp of a multiple of the identity matrix.

% b_colVec = P' * inv(B) * mu   
% But P == Identity
% b_colVec = B \ mu_colVec; %inv(B) * mu_colVec;
b_colVec = B * mu_colVec; % Since B is a scalar multiple of identity, and inv(I)==I

c0 = exp(-.5*(b_colVec'*b_colVec)) * (2*E(1))^(-1/2) * ...
    (2*E(2))^(-1/2) * (2*E(3))^(-1/2);

numTermsInApprox = 150; % Include up to this many terms in approximation. 
d = zeros(1,numTermsInApprox);
c = zeros(1,numTermsInApprox); % Does not include c0. Starts with c1.
% F_list = zeros(1,numTermsInApprox);

F = -1;
convergence = 0; % default: no convergence


sign_terms = ones(1,numTermsInApprox);
sign_terms(1:2:end) = -1*sign_terms(1:2:end);

y_exponents = (p/2)+1:(p/2)+numTermsInApprox;
y_list = y.^y_exponents;
%{
% Keep evaluating F with more terms in its approximation until convergence.
% Congergence: no change in approx of F (to 5 places after decimal) for
%              10 successive approximations.
% If no convergence after considering sums with as many as <numTermsInApprox>
% terms, throw error OR return <convergence==0>.


% for mP = 1:numTermsInApprox
%     
%     maxPowerSeries = mP;   
% 
%     for k = 1:maxPowerSeries
%         % Compute d vector
%         d(k) = (1/2) * ( ( 1-k*(b_colVec(1))^2 )/( (2*E(1))^k ) + ...
%             ( 1-k*(b_colVec(2))^2 )/( (2*E(2))^k ) + ...
%             ( 1-k*(b_colVec(3))^2 )/( (2*E(3))^k )  );
% 
%         % Compute c vector (Eq. 4.2b.5)
%         if k == 1; c(k) = d(k)*c0;else
%             c_inds = 1:k-1; 
%             d_inds = k-1:-1:1;
%             c(k) = (1/k)* (c0*d(k) + sum( d(d_inds).*c(c_inds) ) );     
%         end
%         
% % %         sign_terms = ones(1,maxPowerSeries);
% % %         sign_terms(1:2:end) = -1*sign_terms(1:2:end);
% % % 
% % %         y_exponents = (p/2)+1:(p/2)+maxPowerSeries;
% % %         y_list = y.^y_exponents;
%     end
% 
%         
%     F_list(mP) = (c0*y^(3/2 + 0)/gamma(3/2 + 1)) + ...
%         sum( sign_terms(1:maxPowerSeries).*c(1:maxPowerSeries).*...
%         y_list(1:maxPowerSeries)./gamma(y_exponents(1:maxPowerSeries)+1) );
% 
%     if mP >= 10
%         % Check last 10 terms for congergence
%         F_last10 = F_list(maxPowerSeries-9:maxPowerSeries);
%         if all( abs( diff(F_last10) )  ./ F_last10(1:9)  < 1e-5 ) || all( F_last10 < 1e-10 )
%             convergence = 1;
%             F = F_list(mP);
%             return
%         end
%         
%     end
%     
% end
%}

%chunks = 1:5:numTermsInApprox+1; % i.e. [1    31    61    91   121   151]
step = 3;
chunks = [1, 11, 22:step:numTermsInApprox+1]; % i.e. [1 21 22 23 ... 151]
w = ceil(10/step); w_counter = 0; ind_begin = 1;
Fs = zeros(1,numTermsInApprox);
k=1;
d(k) = (1/2) * ( ( 1-k*(b_colVec(1))^2 )/( (2*E(1))^k ) + ...
    ( 1-k*(b_colVec(2))^2 )/( (2*E(2))^k ) + ...
    ( 1-k*(b_colVec(3))^2 )/( (2*E(3))^k )  );
c(k) = d(k)*c0;

for mP = 1:length(chunks)-1 % i.e. [1 2 3 4 5]
%     chunks
%     mP
    
    % % Evaluate c and d for the current chunk.
    for k = chunks(mP) : chunks(mP+1)-1 % i.e. 1:20,21:21,22:22,23:23
        % Compute d vector
        d(k) = (1/2) * ( ( 1-k*(b_colVec(1))^2 )/( (2*E(1))^k ) + ...
            ( 1-k*(b_colVec(2))^2 )/( (2*E(2))^k ) + ...
            ( 1-k*(b_colVec(3))^2 )/( (2*E(3))^k )  );

        c_inds = 1:k-1; 
        d_inds = k-1:-1:1;
        c(k) = (1/k)* (c0*d(k) + sum( d(d_inds).*c(c_inds) ) );

        % Fill in the next chunk of F evaluations     
        Fs(k) = (c0*y^(3/2 + 0)/gamma(3/2 + 1)) + ...
            sum( sign_terms(1:k).*c(1:k).*...
            y_list(1:k)./gamma(y_exponents(1:k)+1) );
    end
    
    w_counter = w_counter + 1;
    if w_counter == w && mP > 2
       ind_begin = ind_begin + 10; % when ind_end has moved at least 10, update ind_begin forward by 10
       w_counter = 0;
    end
    ind_end = chunks(mP+1)-1; % ind_end changes with each loop iteration
    F_list = Fs( ind_begin:ind_end );

    numConsecutive = 9;
    % % Check for congergence (<numConsecutive> consecutive terms) in 3 ways:
    % % 1.) examining whether successive elements are basically 0
    % % 2.) comparing relative (ratio) changes in successive elements.
    % % 3.) examining whether rescaled successive elements are basically 0 ( this is useful when there are large oscillations)
    % %     In this case, the terms satisfying convergence will not be
    % %     from the 'large' oscillation region of F, but may still in fact
    % %     have oscillations larger than we'd like. In this case, take the
    % %     geometric mean of the abs of the terms as a measure of
    % %     central tendency.
    
    diff_list = abs( diff(F_list)  ./ F_list(1:end-1) ); % compare the relative changes in successive elements.

    logicals_list1 = F_list < 1e-10; % successive elements are basically 0
    logicals_list2 = diff_list < 1e-6; % relative changes in successive elements are small

    converg_start_index_1 = strfind(logicals_list1,ones(1,numConsecutive));
    converg_start_index_2 = strfind(logicals_list2,ones(1,numConsecutive));

    % This may not be successful if we haven't really seen the true 
    % max(F_list) yet. When this test is applied again below, after all 
    % <numTermsInApprox> terms in the approximation of the inf sum have
    % been computed, it should def work.
    F_list_rescaled = F_list / max(abs(F_list));
    logicals_list3 = F_list_rescaled < 1e-6;
    converg_start_index_3 = strfind(logicals_list3,ones(1,numConsecutive));
    

    if ~isempty(converg_start_index_1)
        ind = converg_start_index_1(1);
        F = abs( mean(F_list(ind:ind+(numConsecutive-1) )) ); %F = abs( F_list(ind) );
        convergence = 1;
%         if F <= 1 % as a proper probability should be
%             convergence = 1;
%         else
%             convergence = 0;
%         end
        
%         case1 = case1 + 1;
%         if mP==1;chunk1 = chunk1 + 1;end
%         if mP==2;chunk2 = chunk2 + 1;end
%         if mP==3;chunk3 = chunk3 + 1;end
%         if mP==4;chunk4 = chunk4 + 1;end
%         if mP>=5 && mP<=20;chunk5to20 = chunk5to20 + 1;end
%         if mP>20; chunk21plus = chunk21plus + 1;end
%         save('pr_dist2Gaussians_powerSeries__CASES_STORAGE.mat','case1','case2','case3','case4','case5',...
%             'chunk1','chunk2','chunk3','chunk4','chunk5to20','chunk21plus')  
        return
    elseif ~isempty(converg_start_index_2)
        ind = converg_start_index_2(1);
        F = abs( mean(F_list(ind:ind+(numConsecutive-1) )) ); %F = abs( F_list(ind) );
        convergence = 1;
%         if F <= 1 % as a proper probability should be
%             convergence = 1;
%         else
%             convergence = 0;
%         end
        
%         case2 = case2 + 1;
%         if mP==1;chunk1 = chunk1 + 1;end
%         if mP==2;chunk2 = chunk2 + 1;end
%         if mP==3;chunk3 = chunk3 + 1;end
%         if mP==4;chunk4 = chunk4 + 1;end
%         if mP>=5 && mP<=20;chunk5to20 = chunk5to20 + 1;end
%         if mP>20; chunk21plus = chunk21plus + 1;end
%         save('pr_dist2Gaussians_powerSeries__CASES_STORAGE.mat','case1','case2','case3','case4','case5',...
%             'chunk1','chunk2','chunk3','chunk4','chunk5to20','chunk21plus') 
        return
    elseif ~isempty(converg_start_index_3)
        ind = converg_start_index_3(1);
        F = geomean( abs(F_list(ind:ind+(numConsecutive-1) )) ); %F = abs( F_list_rescaled(ind) );
        if F <= 1 % as a proper probability should be
            convergence = 1;
        else
            convergence = 0;
        end
        %plot(Fs)
        
%         case3 = case3 + 1;
%         if mP==1;chunk1 = chunk1 + 1;end
%         if mP==2;chunk2 = chunk2 + 1;end
%         if mP==3;chunk3 = chunk3 + 1;end
%         if mP==4;chunk4 = chunk4 + 1;end
%         if mP>=5 && mP<=20;chunk5to20 = chunk5to20 + 1;end
%         if mP>20; chunk21plus = chunk21plus + 1;end
%         save('pr_dist2Gaussians_powerSeries__CASES_STORAGE.mat','case1','case2','case3','case4','case5',...
%             'chunk1','chunk2','chunk3','chunk4','chunk5to20','chunk21plus') 
        return
    end

end


    
%plot(Fs)


% If there is still no convergence, it could be the case that
% for medium values of <maxPowerSeries> there are extreme oscillations. But
% for small and large values of <maxPowerSeries> there is *rough*
% convergence - i.e. nothing like the 5decimal agreement of successive
% approximations to F we usually require, but close agreement compared to 
% the oscillations.
% FOR EXAMPLE:
% % mu_colVec = [0.431655350886963;-0.086185260500000;-0.051773981079100];
% % E = [0.023242383257147   0.023242383257147   0.023242383257147]; 
% This behavior holds for any mu_colVec with the same norm:
% i.e. [x,y,z] = sph2cart(2*pi*rand,pi*rand,0.443209641359006); mu_colVec = [x,y,z]'
% See investigating_Pr_dist2Gaussians_powerSeries.m
% A more systematic way to address these issues can be found in 
% Gideon&Gurland(1976)
F_list = Fs;
F_list_rescaled = F_list / max(abs(F_list));
% for f = 1:length(F_list_rescaled)-9
% 
%         % Check last 10 terms for congergence
%         F_last10 = F_list_rescaled(f:f+9);
%         if all( abs( diff(F_last10) )  ./ F_last10(1:9)  < 1e-5 ) || all( F_last10 < 1e-10 )
%             convergence = 1;
%             F = F_list_rescaled(f+9);
%             return
%         end    
%     
% end
% Check for congergence (10 consecutive terms)
diff_list = abs( diff(F_list_rescaled)  ./ F_list_rescaled(1:end-1) );

logicals_list1 = F_list_rescaled < 1e-11;
logicals_list2 = diff_list < 1e-6; % elements

converg_start_index_1 = strfind(logicals_list1,ones(1,numConsecutive));
converg_start_index_2 = strfind(logicals_list2,ones(1,numConsecutive));
if ~isempty(converg_start_index_1)
    ind = converg_start_index_1(1);
    F = geomean( abs(F_list(ind:ind+(numConsecutive-1) )) ); %F = F_list_rescaled(ind);
    convergence = 1;
    
%     case4 = case4 + 1;
%     save('pr_dist2Gaussians_powerSeries__CASES_STORAGE.mat','case1','case2','case3','case4','case5',...
%         'chunk1','chunk2','chunk3','chunk4','chunk5to20','chunk21plus')    

elseif ~isempty(converg_start_index_2)
    ind = converg_start_index_2(1);
    F = geomean( abs(F_list(ind:ind+(numConsecutive-1) )) ); %F = F_list_rescaled(ind);
    convergence = 1;
    
%     case5 = case5 + 1;
%     save('pr_dist2Gaussians_powerSeries__CASES_STORAGE.mat','case1','case2','case3','case4','case5',...
%         'chunk1','chunk2','chunk3','chunk4','chunk5to20','chunk21plus')    
else
end



% assert(0==1)

end






