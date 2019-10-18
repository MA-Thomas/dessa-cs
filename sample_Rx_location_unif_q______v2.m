function [ rA, thetaA,counter_rA,counter_thetaA ] = ...
    sample_Rx_location_unif_q______v2( t_elapsed_a,t_elapsed_b,...
    Da_actual,Db_actual,d,n_sigma )
% REJECTION SAMPLING USED TO OBTAIN rA_sample, thetaA_sample
% The proposal distribution q is UNIFORM DISTRIBUTION 

% Inputs: the A particle is the reactant with the SMALLER diffusion sphere
%         at reaction time. The B particle is that with the larger diff
%         sphere at reaction time.
%         <d> is the distance between the last known positions of A and B.
%         <t_elapsed_a> is the time elapsed between A's last state update and
%         the current reaction time. Analogous for B.

% Why use t_elapsed_a and t_elapsed_b instead of the waiting time?
% There are 2 situations to consider.
% Situation 1: A and B's states were last updated at the same time. They
%              both began diffusing at the same time. 
% Situation 2: Reaction A+B was evaluated when one of the particles had
%              already been diffusing (one particle has t_offset > 0, the
%              other has t_offset == 0)
% At most one of the particles can have a positive offset since we only
% consider reactions (to add to the PQ) involving a particle that was just 
% updated.
% In situation 2, the lesser of t_elapsed_a and t_elapsed_b is the 
% wait_time.


% 1. Compute diffusion sphere radii Ra,Rb at moment of bimolecular reaction
%    i.e. the "observed" diffusion sphere radii
Ra = n_sigma.*(6.*Da_actual.*t_elapsed_a)^(1./2);  
Rb = n_sigma.*(6.*Db_actual.*t_elapsed_b)^(1./2);


% "A" labeled the particle with the smaller diffusion sphere just before
% the wait_time was sampled. "B" labeled the particle with the larger
% diffusion sphere.
% HOWEVER, now the time each particle has been diffusing, t_elapsed_x
% includes that wait_time and the labels may no longer be appropriate.
% Switch them if necessary.
if Ra > Rb
   Da_temp = Da_actual;
   Da_actual = Db_actual;
   Db_actual = Da_temp;
   
   t_elapsed_temp = t_elapsed_a;
   t_elapsed_a = t_elapsed_b;
   t_elapsed_b = t_elapsed_temp;
   
%     Da_actual
%     Db_actual
%     t_elapsed_a
%     t_elapsed_b
    Ra = n_sigma*(6*Da_actual*t_elapsed_a)^(1/2);  
    Rb = n_sigma*(6*Db_actual*t_elapsed_b)^(1/2);
end

assert(Ra<=Rb);
    
    
% 2. Using d, Ra, Rb, infer Case (i.e. 1, 2, 3/4, or 5).
if d > Ra && d > Rb
    caseNum = 1; % Overlap Volume contains neither initial position of A nor B.
    
elseif d > Ra && d <= Rb && 2*Ra > Rb
    caseNum = 2; % Overlap Volume contains initial position of A, but not of B.
                 % Also, diameter of A is greater than radius of B.
elseif d > Ra && d <= Rb && 2*Ra <= Rb
    caseNum = 3; % Overlap Volume contains initial position of A, but not of B.
                 % Diameter of A is now less than radius of B (B grows fast).    
elseif d <= Ra && d <= Rb && Ra + d > Rb
    caseNum = 4; % Overlap Volume contains initial positions of A and of B.
                 % Part of the A sphere is still outside the B sphere.
else
    caseNum = 5;
end
    
% 3. If the diffusion spheres of molecules A and B started growing at the
%    same absolute time (no time offset, both states were last updated at 
%    the same time), use their actual Da and Db values.
%    However, if there was some offset, the actual Da,Db values would not 
%    lead to the observed Ra,Rb values. Therefore, infer temporary Da,Db 
%    values from the observed Ra,Rb at reaction time.

%    ********
%    THE INFERRED Da OR Db ENFORES THE ASSUMPTION THAT BOTH PARTICLES HAVE 
%    DIFFUSED FOR THE SAME DURATION. THAT ASSUMPTION IS NECESSARY FOR 
%    SAMPLING A REACTION LOCATION BY MEANS OF CASE, I.E. OF SPHERE-SPHERE 
%    OVERLAP STAGE.
%    ********


% Situation 1. No inference of Da or Db needed.
if t_elapsed_a == t_elapsed_b
    Da = Da_actual;
    Db = Db_actual;
    
    t = t_elapsed_a; % no need to infer diffusion time. It's just the  
                     % elapsed time since A (and B) was last updated.
% Situation 2.    
else
    % Infer Da (Db) if particle A (B) was diffusing longer.
    [t_inferred,Ind_a_or_b] = min([t_elapsed_a, t_elapsed_b]);
    if Ind_a_or_b == 1 % Particle B had an offset. I.e. it was diffusing for longer
        Da = Da_actual;
        Db = Rb^2 / (6 * t_inferred * n_sigma^2); % Inferred Db

    elseif Ind_a_or_b == 2 % Particle A had an offset. I.e. it was diffusing for longer
        Db = Db_actual;
        Da = Ra^2 / (6 * t_inferred * n_sigma^2); % Inferred Da
    end
    
    t = t_inferred; % Inferred diffusion time assuming both A and B began
                    % diffusing at the same instant.  
                    % This variable is used below.
end


% -------------------------------------------------------------------------
% Code below is from the original version (sample_Rx_location_unif_q.m)
% that assumed both A and B have diffused for the same duration
% -------------------------------------------------------------------------

Xint_B = (d^2 - Ra^2 + Rb^2)./(2.*d); % This is the distance from B to the intersection point
h_Bcap = Rb - Xint_B;
Xint_A = d-Xint_B; % This is the distance from A to the intersection point.
% t = wait_time;
counter_rA = 0;
counter_thetaA = 0;

if caseNum == 1    
    % % SAMPLE rA ----
    
    reject = 1;
    accepted = 0;
    while reject == 1        
        counter_rA = counter_rA+1;
        
        r_lb = d-Rb;%Xint_A - h_Bcap;
        r_ub = Ra;
        
        if r_lb < 0 % OLD: Likely because Db/Da > 4 and loaded propensity CDFs are for Db_Da_Ratios under 4.
            d
            Rb
            assert(r_lb>0)
        end
        
        % Use simple uniform proposal distribution q. It should be easy to
        % determine a constant c which guarantees that for any rA,
        % c*q(rA) >= p_rewieghted(rA) 
        
        % SAMPLE rA ~ q
        rA = (r_ub - r_lb)*rand(1,1) + r_lb;% q = uniform[r_lb,r_ub]
        
        % Compute true probability of rA.
        % Start with un-weighted probability of seeing rA 
        P_rA_given_wait_time = P_r__C1(Da, t, rA, Ra, Xint_A,h_Bcap);        

        denom1 = (1/4*d)*(-Rb^2 + 6*Da*t + d^2)*erf(r_ub/(12*Da*t)^(1/2)) - ...
            (1/4*d)*(-Rb^2 + 6*Da*t + d^2)*erf(r_lb/(12*Da*t)^(1/2));
        
        denom2 = (((3*Da*t/pi)^(1/2)) / (2*d)) * exp(-r_ub^2/(12*Da*t)) * (r_ub-2*d) - ...
            (((3*Da*t/pi)^(1/2)) / (2*d)) * exp(-r_lb^2/(12*Da*t)) * (r_lb-2*d);
        w_rA = rA*( ((rA^2 + d^2 - Rb^2)/(2*rA*d)) - 1 )/ (denom2-denom1);
                   
        P_rA_re_weighted = P_rA_given_wait_time * w_rA;
        
        % Sample Y|rA from uniform( 0,c*q(rA) ) 
        p_MAX_in_OV = P_r__C1(Da, t, r_lb, Ra, Xint_A,h_Bcap);
        c = w_rA * p_MAX_in_OV;
        Y = (c*w_rA)*rand(1,1);
        

        % ACCEPT or REJECT
        % If Y is less than the actual probability of seeing rA, ACCEPT
        if Y <= P_rA_re_weighted; accepted = 1;end

        % TEMP
        if counter_rA > 100; reject = 0; end
        
        % If accepted, exit while loop.
        if accepted == 1
           reject = 0; 
        end
       
    end
    
    % % SAMPLE thetaA ---- rA has been set.
    
    reject = 1;
    accepted = 0;
    while reject == 1        
        counter_thetaA = counter_thetaA+1; 
        
        theta_lb = 0;
        theta_ub = acos((rA^2 + d^2 - Rb^2) / (2*rA*d));
        
        % Use simple uniform proposal distribution q. It should be easy to
        % determine a constant c which guarantees that for any rA,
        % c*q(rA) >= P_thetaA(rA) 
        
        % SAMPLE rA ~ q
        thetaA = (theta_ub - theta_lb)*rand(1,1) + theta_lb;
        
        % Compute true probability of thetaA (i.e. pr(rB)*ring_circumf ).
        % Start with un-weighted probability of seeing rB
        rB = (rA + d^2 - 2*rA*d*cos(thetaA))^(1/2);
        P_rB_given_wait_time = P_r__C1(Db, t, rB, Ra, Xint_A,h_Bcap);        

        ring_Circumf = 2*pi*rA*sin(thetaA);
                   
        P_thetaA = P_rB_given_wait_time * ring_Circumf;
        
 
        K = 12*Db*t;
        Z = 2*rB*d;
        W = 1 / sqrt(12*pi*Db*t);        
        % P(theta|r,t) = constant*exp(-(W^2 - Zcos(theta))/K)*sin(theta).
        % We want the theta value that maximizes this probability in the
        % range [0,pi/2] since we know the OV in Case 1 will never include
        % theta beyond pi/2.
        thetaMAX = acos( (sqrt(K^2 + 4*Z^2) - K)/(2*Z) ); 
        
        % The argmax of P(theta|r,t) within [0,pi/2] may be a theta outside
        % the OV. If so, replace it with the nearest theta in the OV.
        if thetaMAX > theta_ub; thetaMAX=theta_ub; end
        
        rMAX = (rA + d^2 - 2*rA*d*cos(thetaMAX))^(1/2);
        MAXp = P_r__C1(Da, t, rMAX, Ra, Xint_A,h_Bcap)*2*pi*rA*sin(thetaMAX);
        
        % Sample Y|rA from uniform( 0,cq(rA)=MAX )
        Y = (MAXp)*rand(1,1);
        

        % ACCEPT or REJECT
        % If Y is less than the actual probability of seeing rA, ACCEPT
        if Y <= P_thetaA; accepted = 1;end
 
        % TEMP
        if counter_thetaA > 100; reject = 0; end
        
        % If accepted, exit while loop.
        if accepted == 1
           reject = 0; 
        end
       
    end
    
elseif caseNum == 2
    % % SAMPLE rA ----    
    
    reject = 1;
    accepted = 0;
    while reject == 1
        counter_rA = counter_rA+1;
        r_lb = 0;
        r_ub = Ra;     
        r_mb = (Rb-d);
        
        % Use simple uniform proposal distribution q. It should be easy to
        % determine a constant c which guarantees that for any rA,
        % c*q(rA) >= p_rewieghted(rA) 
        
        % SAMPLE rA ~ q
        rA = (r_ub - r_lb)*rand(1,1) + r_lb;% q = uniform[r_lb,r_ub]
        
        % Compute true probability of rA.
        % Start with un-weighted probability of seeing rA 
        P_rA_given_wait_time = P_r__C2_5(Da, t, rA);        

        r_max = max(rA,(Rb-d));
        
        denom1 = (12*Da*t/pi)^(1/2)*(exp(-r_mb^2 / (12*Da*t)) - exp(-r_lb^2 / (12*Da*t)));
        
        denom2 = (1/4*d)*(-Rb^2 + 6*Da*t + d^2)*erf(r_ub/(12*Da*t)^(1/2)) - ...
            (1/4*d)*(-Rb^2 + 6*Da*t + d^2)*erf(r_lb/(12*Da*t)^(1/2));
        
        denom3 = (((3*Da*t/pi)^(1/2)) / (2*d)) * exp(-r_ub^2/(12*Da*t)) * (r_ub-2*d) - ...
            (((3*Da*t/pi)^(1/2)) / (2*d)) * exp(-r_mb^2/(12*Da*t)) * (r_mb-2*d);
        
        w_rA = -rA*( ((r_max^2 + d^2 - Rb^2)/(2*r_max*d)) - 1 )/...
            (-denom1-denom2+denom3);
                   
        w_rA = abs(w_rA); %MT (1/15/2019): sometimes w_rA is neg and this is bad. not sure why.
        
        P_rA_re_weighted = P_rA_given_wait_time * w_rA;

        % Sample Y|rA from uniform( 0,c*q(rA) ) 
        p_MAX_in_OV = P_r__C2_5(Da, t, r_lb);
        c = w_rA * p_MAX_in_OV;
        Y = (c)*rand(1,1);        
        assert(c>=P_rA_re_weighted) %ensure that envolope distrib always greater than or equal to actual distrib.
        
        % ACCEPT or REJECT
        % If Y is less than the actual probability of seeing rA, ACCEPT
        if Y <= P_rA_re_weighted; accepted = 1;end
        
        % If accepted, exit while loop.
        if accepted == 1
           reject = 0; 
        end
       
        % TEMP
        if counter_rA > 100; reject = 0; end
        
    end

    % % SAMPLE thetaA ---- rA has been set.    
    reject = 1;
    accepted = 0;
    while reject == 1        
        counter_thetaA = counter_thetaA+1;
        
        theta_lb = 0;
        if rA <= (Rb-d)
            theta_ub = pi;
        else
            theta_ub = acos((rA^2 + d^2 - Rb^2) / (2*rA*d));

            assert(Rb>d) % this should be true in Case2
            %assert(rA<Rb-d) % If you want the code to stop here, uncomment.
        end
        
        % Use simple uniform proposal distribution q. It should be easy to
        % determine a constant c which guarantees that for any rA,
        % c*q(rA) >= P_thetaA(rA) 
        
        % SAMPLE rA ~ q (uniformly)
        thetaA = (theta_ub - theta_lb)*rand(1,1) + theta_lb;
        
        % Compute true probability of thetaA (i.e. pr(rB) ).
        % Start with un-weighted probability of seeing rB
        rB = (rA + d^2 - 2*rA*d*cos(thetaA))^(1/2);
        P_rB_given_wait_time = P_r__C2_5(Db, t, rB);
        

        ring_Circumf = 2*pi*rA*sin(thetaA);                
        P_thetaA = P_rB_given_wait_time * ring_Circumf;
        
 
        K = 12*Db*t;
        Z = 2*rB*d;
        W = 1 / sqrt(12*pi*Db*t);        
        % P(theta|r,t) = constant*exp(-(W^2 - Zcos(theta))/K)*sin(theta).
        % We want the theta value that maximizes this probability in the
        % range [0,pi/2] since we know the OV in Case 1 will never include
        % theta beyond pi/2.
        thetaMAX = acos( (sqrt(K^2 + 4*Z^2) - K)/(2*Z) ); 
        
        % The argmax of P(theta|r,t) within [0,pi/2] may be a theta outside
        % the OV. If so, replace it with the nearest theta in the OV.
        if thetaMAX > theta_ub; thetaMAX=theta_ub; end
        
        rMAX = (rA + d^2 - 2*rA*d*cos(thetaMAX))^(1/2);
        MAXp = P_r__C2_5(Da, t, rMAX)*2*pi*rA*sin(thetaMAX);
        
        % Sample Y|rA from uniform( 0,cq(rA)=MAX )
        Y = (MAXp)*rand(1,1);
        

        % ACCEPT or REJECT
        % If Y is less than the actual probability of seeing rA, ACCEPT
        if Y <= P_thetaA; accepted = 1;end
        
        % If accepted, exit while loop.
        if accepted == 1
           reject = 0; 
        end
       
        % TEMP
        if counter_thetaA > 100; reject = 0; end
    end

    
    
elseif caseNum == 3
    % % SAMPLE rA ----
    
    % NO REJECTION SAMPLING NEEDED.
    % SAMPLE rA ~ [original probabiilty density]
    rA = abs( normrnd(0, (6*Da*t)^(1/2)) );  
    
    % Ensure sample is within diffusion intersection volume
    if rA > min(Ra,Rb)
       rA = 0.97 * min(Ra,Rb);
       
       % For testing purposes, throw error. After testing comment this out
       %assert(0==1)
    end

    % % SAMPLE thetaA ---- rA has been set.    
    reject = 1;
    accepted = 0;
    while reject == 1        
        counter_thetaA = counter_thetaA+1;
        
        theta_lb = 0;
        theta_ub = pi;
        
        % Use simple uniform proposal distribution q. It should be easy to
        % determine a constant c which guarantees that for any rA,
        % c*q(rA) >= P_thetaA(rA) 
        
        % SAMPLE rA ~ q (uniformly)
        thetaA = (theta_ub - theta_lb)*rand(1,1) + theta_lb;
        
        % Compute true probability of thetaA (i.e. pr(rB) ).
        % Start with un-weighted probability of seeing rB
        rB = (rA + d^2 - 2*rA*d*cos(thetaA))^(1/2);
        P_rB_given_wait_time = P_r__C2_5(Db, t, rB);

        ring_Circumf = 2*pi*rA*sin(thetaA);                   
        P_thetaA = P_rB_given_wait_time * ring_Circumf;
        
 
        K = 12*Db*t;
        Z = 2*rB*d;
        W = 1 / sqrt(12*pi*Db*t);        
        % P(theta|r,t) = constant*exp(-(W^2 - Zcos(theta))/K)*sin(theta).
        % We want the theta value that maximizes this probability in the
        % range [0,pi/2] since we know the OV in Case 1 will never include
        % theta beyond pi/2.
        thetaMAX = acos( (sqrt(K^2 + 4*Z^2) - K)/(2*Z) ); 
        
        % The argmax of P(theta|r,t) within [0,pi/2] may be a theta outside
        % the OV. If so, replace it with the nearest theta in the OV.
        if thetaMAX > theta_ub; thetaMAX=theta_ub; end
        
        rMAX = (rA + d^2 - 2*rA*d*cos(thetaMAX))^(1/2);
        MAXp = P_r__C2_5(Da, t, rMAX)*2*pi*rA*sin(thetaMAX);
        
        % Sample Y|rA from uniform( 0,cq(rA)=MAX )
        Y = (MAXp)*rand(1,1);
        

        % ACCEPT or REJECT
        % If Y is less than the actual probability of seeing rA, ACCEPT
        if Y <= P_thetaA; accepted = 1;end
        
        % If accepted, exit while loop.
        if accepted == 1
           reject = 0; 
        end
        
        % TEMP
        if counter_thetaA > 100; reject = 0; end       
    end

    
    
elseif caseNum == 4    
    % % SAMPLE rA ---- EXACTLY AS IN CASE 2
        
    reject = 1;
    accepted = 0;
    while reject == 1
        counter_rA = counter_rA+1;
        
        r_lb = 0;
        r_ub = Ra;       
        r_mb = (Rb-d);
        
        % Use simple uniform proposal distribution q. It should be easy to
        % determine a constant c which guarantees that for any rA,
        % c*q(rA) >= p_rewieghted(rA) 
        
        % SAMPLE rA ~ q
        rA = (r_ub - r_lb)*rand(1,1) + r_lb;% q = uniform[r_lb,r_ub]
        
        % Compute true probability of rA.
        % Start with un-weighted probability of seeing rA 
        P_rA_given_wait_time = P_r__C2_5(Da, t, rA);        

        r_max = max(rA,(Rb-d));
        
        denom1 = (12*Da*t/pi)^(1/2)*(exp(-r_mb^2 / (12*Da*t)) - exp(-r_lb^2 / (12*Da*t)));
        
        denom2 = (1/4*d)*(-Rb^2 + 6*Da*t + d^2)*erf(r_ub/(12*Da*t)^(1/2)) - ...
            (1/4*d)*(-Rb^2 + 6*Da*t + d^2)*erf(r_lb/(12*Da*t)^(1/2));
        
        denom3 = (((3*Da*t/pi)^(1/2)) / (2*d)) * exp(-r_ub^2/(12*Da*t)) * (r_ub-2*d) - ...
            (((3*Da*t/pi)^(1/2)) / (2*d)) * exp(-r_mb^2/(12*Da*t)) * (r_mb-2*d);
        
        w_rA = abs(-rA*( ((r_max^2 + d^2 - Rb^2)/(2*r_max*d)) - 1 )/...
            (-denom1-denom2+denom3) );
                   
        P_rA_re_weighted = P_rA_given_wait_time * w_rA;
        
        % Sample Y|rA from uniform( 0,c*q(rA) ) 
        p_MAX_in_OV = P_r__C2_5(Da, t, r_lb);
        c = w_rA * p_MAX_in_OV;
        Y = (c)*rand(1,1);        
        assert(c>=P_rA_re_weighted) %ensure that envolope distrib always greater than or equal to actual distrib.
        
        % ACCEPT or REJECT
        % If Y is less than the actual probability of seeing rA, ACCEPT
        if Y <= P_rA_re_weighted; accepted = 1;end
        
        % If accepted, exit while loop.
        if accepted == 1
           reject = 0; 
        end
        
        % TEMP
        if counter_rA > 100; reject = 0; end       
    end    

    % % SAMPLE thetaA ---- rA has been set.    
    reject = 1;
    accepted = 0;
    while reject == 1        
        counter_thetaA = counter_thetaA+1;
        
        theta_lb = 0;
        if rA <= (Rb-d)
            theta_ub = pi;
        else
            theta_ub = acos((rA^2 + d^2 - Rb^2) / (2*rA*d));
        end
        
        % Use simple uniform proposal distribution q. It should be easy to
        % determine a constant c which guarantees that for any rA,
        % c*q(rA) >= P_thetaA(rA) 
        
        % SAMPLE rA ~ q (uniformly)
        thetaA = (theta_ub - theta_lb)*rand(1,1) + theta_lb;
        
        % Compute true probability of thetaA (i.e. pr(rB) ).
        % Start with un-weighted probability of seeing rB
        rB = (rA + d^2 - 2*rA*d*cos(thetaA))^(1/2);
        P_rB_given_wait_time = P_r__C2_5(Db, t, rB);

        ring_Circumf = 2*pi*rA*sin(thetaA);                   
        P_thetaA = P_rB_given_wait_time * ring_Circumf;
        
 
        K = 12*Db*t;
        Z = 2*rB*d;
        W = 1 / sqrt(12*pi*Db*t);        
        % P(theta|r,t) = constant*exp(-(W^2 - Zcos(theta))/K)*sin(theta).
        % We want the theta value that maximizes this probability in the
        % range [0,pi/2] since we know the OV in Case 1 will never include
        % theta beyond pi/2.
        thetaMAX = acos( (sqrt(K^2 + 4*Z^2) - K)/(2*Z) ); 
        
        % The argmax of P(theta|r,t) within [0,pi/2] may be a theta outside
        % the OV. If so, replace it with the nearest theta in the OV.
        if thetaMAX > theta_ub; thetaMAX=theta_ub; end
        
        rMAX = (rA + d^2 - 2*rA*d*cos(thetaMAX))^(1/2);
        MAXp = P_r__C2_5(Da, t, rMAX)*2*pi*rA*sin(thetaMAX);
        
        % Sample Y|rA from uniform( 0,cq(rA)=MAX )
        Y = (MAXp)*rand(1,1);
        

        % ACCEPT or REJECT
        % If Y is less than the actual probability of seeing rA, ACCEPT
        if Y <= P_thetaA; accepted = 1;end
        
        % If accepted, exit while loop.
        if accepted == 1
           reject = 0; 
        end
        
        % TEMP
        if counter_thetaA > 100; reject = 0; end       
    end

   
    
elseif caseNum == 5
    % % SAMPLE rA ----
    
    % NO REJECTION SAMPLING NEEDED.
    % SAMPLE rA ~ [original probabiilty density]
    rA = abs( normrnd(0, (6*Da*t)^(1/2)) ); 
    
    % Ensure sample is within diffusion intersection volume
    if rA > min(Ra,Rb)
       rA = 0.97 * min(Ra,Rb);
       
       % For testing purposes, throw error. After testing comment this out
       %assert(0==1)
    end
    

    % % SAMPLE thetaA ---- rA has been set.    
    reject = 1;
    accepted = 0;
    while reject == 1        
        counter_thetaA = counter_thetaA+1;
        
        theta_lb = 0;
        theta_ub = pi;
        
        % Use simple uniform proposal distribution q. It should be easy to
        % determine a constant c which guarantees that for any rA,
        % c*q(rA) >= P_thetaA(rA) 
        
        % SAMPLE rA ~ q (uniformly)
        thetaA = (theta_ub - theta_lb)*rand(1,1) + theta_lb;
        
        % Compute true probability of thetaA (i.e. pr(rB) ).
        % Start with un-weighted probability of seeing rB
        rB = (rA + d^2 - 2*rA*d*cos(thetaA))^(1/2);
        P_rB_given_wait_time = P_r__C2_5(Db, t, rB);

        ring_Circumf = 2*pi*rA*sin(thetaA);                   
        P_thetaA = P_rB_given_wait_time * ring_Circumf;
        
 
        K = 12*Db*t;
        Z = 2*rB*d;
        W = 1 / sqrt(12*pi*Db*t);        
        % P(theta|r,t) = constant*exp(-(W^2 - Zcos(theta))/K)*sin(theta).
        % We want the theta value that maximizes this probability in the
        % range [0,pi/2] since we know the OV in Case 1 will never include
        % theta beyond pi/2.
        thetaMAX = acos( (sqrt(K^2 + 4*Z^2) - K)/(2*Z) ); 
        
        % The argmax of P(theta|r,t) within [0,pi/2] may be a theta outside
        % the OV. If so, replace it with the nearest theta in the OV.
        if thetaMAX > theta_ub; thetaMAX=theta_ub; end
        
        rMAX = (rA + d^2 - 2*rA*d*cos(thetaMAX))^(1/2);
        MAXp = P_r__C2_5(Da, t, rMAX)*2*pi*rA*sin(thetaMAX);
        
        % Sample Y|rA from uniform( 0,cq(rA)=MAX )
        Y = (MAXp)*rand(1,1);
        

        % ACCEPT or REJECT
        % If Y is less than the actual probability of seeing rA, ACCEPT
        if Y <= P_thetaA; accepted = 1;end
        
        % If accepted, exit while loop.
        if accepted == 1
           reject = 0; 
        end
        
        % TEMP
        if counter_thetaA > 100; reject = 0; end       
    end

end
end

function P = P_r__C1(Da, wait_time, r, Ra, Xint_A, h_Bcap)
% Un-weighted original isotripic probability density, Case 1.
% Probability density is 0 until diffusion spheres overlap.

    if r >= Xint_A - h_Bcap && r <= Ra        
        
        P = (1 / (12*Da*wait_time*pi)^(1/2) )*...
            exp(-r^2 / (12*Da*wait_time));
    else
        P = 0;
    end
end

function P = P_r__C2_5(Da, wait_time, r)
% Un-weighted original isotripic probability density, Cases 2 to 5 

        P = (1 / (12*Da*wait_time*pi)^(1/2) )*...
            exp(-r^2 / (12*Da*wait_time));
end

