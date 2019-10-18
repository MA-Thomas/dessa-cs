% IF REACTION IS BIMOLECULAR    
% This is v2. PriorityQ stores absolute times of next reactions rather than
% wait time. 

bimolecular_count = bimolecular_count + 1;

% Set locations of reactants (to be used if reflective BC are needed)
xyz_Global_A = allAssembliesXYZ(:,row_col_logicals_index_a1);
xyz_Global_B = allAssembliesXYZ(:,row_col_logicals_index_a2);

% % Sample reaction location given waiting time.
d = pdist([xyz_Global_A';xyz_Global_B'],'euclidean');
    
% Da = (kB_T / (6*pi*eta*currentStokesRs(row_col_logicals_index_a1)) );   
% Db = (kB_T / (6*pi*eta*currentStokesRs(row_col_logicals_index_a2)) );
Da = Dmonomer;
Db = Dmonomer; % as in Chew et al.

[ rA, thetaA, ~,~ ] = sample_Rx_location_unif_q______v2( t_elapsed_a,t_elapsed_b,...
    Da,Db,d,n_sigma ); % caseNum calc within the function

% % Convert to XYZ coords.
% x-axis is defined by the line connecting particles A and B
x = real( rA * cos(thetaA) );

% Uniformly sample location on the equiprobability ring 
% (the ring lies on the y-z plane in the reference frame centered on A)
thetaRing = (2*pi)*rand(1,1); %unif sampling on [0,2pi]
ringradius = rA * sin(thetaA);
y = real( ringradius*cos(thetaRing) );
z = real( ringradius*sin(thetaRing) );

% % Update Assembly IDs to show A and B are now part of same assembly.
% Create a new assembly ID and assign it to both particles.
currentMaxID = currentMaxID + 1;
% 1) for each subunit involved as part of a reactant, assign it the new
% assemblyID
subunitIDsAssemblyIDs(subunitIDsAssemblyIDs==a1) = currentMaxID;
subunitIDsAssemblyIDs(subunitIDsAssemblyIDs==a2) = currentMaxID;
% 2) update currentAssemblyIDs (assemblies corresponding to cols of
%    assemblyDistMatrix).  assemblyDistMatrix is loosing 2
%    columns(/rows) worth of entries -> 1 column for each assembly that
%    reacted. It is also gaining 1 column's worth due to the reaction
%    product. It's arbitrary which column we choose for the reaction 
%    product and which we disregard until needed later. 
%    Choose to disregard the faster diffusing assembly.
currentAssemblyIDs(row_col_logicals_index_a1) = currentMaxID;
currentAssemblyIDs(row_col_logicals_index_a2) = inf;
% % WE CAN NO LONGER USE THE VARs "a1"/"a2" in the current loop.

% Update list of Stokes-Einstein Radii.
% % % currentStokesRs(row_col_logicals_index_a1) = ...
% % %     currentStokesRs(row_col_logicals_index_a1) +...
% % %     currentStokesRs(row_col_logicals_index_a2); % New Radius is sum of those of reactant assemblies
% -------------------------------------------------------------------------
% For Michaelis-Menten model, all assemblies have D ~ 1, therefore, all
% assemblies have roughly the same StokesEinsteinRadius.
% -------------------------------------------------------------------------
currentStokesRs(row_col_logicals_index_a1) = StokesEinsteinRadius + abs( normrnd(0,0.01*StokesEinsteinRadius) );
currentStokesRs(row_col_logicals_index_a2) = inf;

% Update absolute time, the time the most recent reaction occured.
t_absolute = t_reaction;   


% % % Update the position of the two reactants.
% % % One reactant's location becomes the product location, the other is ignored.

% See: https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d        
% 1) Perform a rotation of [x y z] so that it is describing the
% same point, but with respect to the global reference frame's axes.
% Thus, I need a matrix, RotMat, that rotates the vector AB into x-hat.
ABvec = allAssembliesXYZ(:,row_col_logicals_index_a2)' - allAssembliesXYZ(:,row_col_logicals_index_a1)';        
ABvec = ABvec / norm(ABvec);
xAxis = [1 0 0]';
v = cross(ABvec,xAxis);
v_skew_symm_matrix = [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];
s = norm(v);
c = dot(ABvec,xAxis);
RotMat = eye(3) + v_skew_symm_matrix + ...
    (v_skew_symm_matrix*v_skew_symm_matrix)*( (1-c)/(s^2) );

% 2) Rotate the reaction location vector 
xyz_aligned_with_Global_axes = RotMat*[x y z]';

% 3) Perform translation so that reaction point (x,y,z) shifts to
%    being described w.r.t. the true global reference frame's origin.
xyz_Global = allAssembliesXYZ(:,row_col_logicals_index_a1) + xyz_aligned_with_Global_axes;  

% Apply Reflective Boundary Conditions to ensure xyz_Global is in
% simulation box
if any(abs(xyz_Global) > containerLength/2)
    
    xyz_Global = testing_reflective_boundary_conditions_Bimol(...
        xyz_Global_A,xyz_Global_B,xyz_Global,containerLength);
end

% 4) Store updated position of the product assembly in the a1 col.
%    The other reactant assembly's col (a2 col) is ignored. 
%    However, a later monomolecular breaking reaction may make use
%    of the col.
allAssembliesXYZ(:,row_col_logicals_index_a1) = xyz_Global;
allAssembliesXYZ(:,row_col_logicals_index_a2) = inf;%xyz_Global;


% --------------------------------------------------
% Begin Update Michaelis Menten benchmark parameters
% --------------------------------------------------
% All bimolecular reactions are E + S -> ES
currentAssemblyTYPES(row_col_logicals_index_a1) = 3;
currentAssemblyTYPES(row_col_logicals_index_a2) = -1;

% ------------------------------------------------
% End Update Michaelis Menten benchmark parameters
% ------------------------------------------------


% Keep track of the absolute time this product was born, i.e. the
% time its state was last updated.
assemblies_lastUpdate_absTime(row_col_logicals_index_a1) = t_absolute; % updated
assemblies_lastUpdate_absTime(row_col_logicals_index_a2) = inf; % ignored








% % ----------------------------------------------------
% % Update the Priority Queue with potential new events.
% % ----------------------------------------------------
% 1. Bimolecular event involving the new product and 1 other assembly. 
%    If none, temporarily store position only updates info for this channel.
% FOR MICHAELIS-MENTEN, BIMOL EVENT NOT POSSIBLE HERE AS DIMER IS MAX
% OLIGOMER SIZE.

% Let rows cover the single new product.
% Let cols cover all other assemblies.
row_single_product = find(currentAssemblyIDs==currentMaxID);


for row = row_single_product
    %%{
    % Compute list of distances to all other molelcules
    assemblyDistVector = real( (bsxfun(@plus,dot(allAssembliesXYZ,allAssembliesXYZ,1)',...
        dot(allAssembliesXYZ(:,row),allAssembliesXYZ(:,row),1))-...
        2*(allAssembliesXYZ'*allAssembliesXYZ(:,row)) ).^(1/2) )';
    
    % If this reactant is not allowed to react with any others, skip. 
    if ~any(assemblyDistVector < inf);continue;end 

    
    % Iterate over upperTriangular part of assemblyDistMatrix only.
    % And, only consider columns with:
    % 1. distances less than <maxRxDistance>
    % 2. distances greater than <maxRxDistance> only if partner has a
    %    nonzero t_offset and d < maxRxDistance+R(t_offset)
    %    i.e. partner has been diffusing long enough to reach distance
    %    of <maxRxDistance>.
    columns = 1:length(assemblyDistVector);
    col_indices1 = assemblyDistVector(columns) <= maxRxDistance;
     
    Rs = real( n_sigma*(6*(kB_T ./ (6*pi*eta*currentStokesRs(columns)) )*...
        t_absolute - assemblies_lastUpdate_absTime(columns)).^(1/2) );    
    
    col_indices2 = assemblyDistVector(columns) <= maxRxDistance+Rs;    
    col_indices = col_indices1 | col_indices2;  

    % 3. Michaelis-Menten allows 1 bimol Rx: E + S -> ES, i.e. 1 + 2 -> ...
    if currentAssemblyTYPES(row) == 1
        col_indices_ReactionType = currentAssemblyTYPES(columns) == 2;
        col_indices = col_indices & col_indices_ReactionType;
    elseif currentAssemblyTYPES(row) == 2
        col_indices_ReactionType = currentAssemblyTYPES(columns) == 1;
        col_indices = col_indices & col_indices_ReactionType;
    else
        col_indices = col_indices & 0; % currentAssemblyTYPES(row) is not E or S
    end
    
   
    % ********** BEGIN NEW CODE ************
    cL = (containerLength/2);
    
    Cols = columns(col_indices); % i.e. [ 3 4 11 13 ]
    distList__currRow = assemblyDistVector(Cols); % i.e. [ 1.2 3.01 5.7 1.3 ]
    
    % Calculate distance to the nearest boundary.
    curr_XYZ_colVec = allAssembliesXYZ(:,row);
    assert( sum( abs(curr_XYZ_colVec) > cL*ones(size(curr_XYZ_colVec)) ) ==0) 
    
    partners_XYZ_colVecs =  allAssembliesXYZ(:,Cols); 
    
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    % For the current product and its partners, we need to evaluate the
    % distances to their nearest relevant boundaries. The relevant boundaries might
    % be at the reflective boundary or imaginary boundary 1 (d_inner1 from
    % reflective boundary) or imaginary boundary 2 (d_inner2 from
    % reflective boundary).
    three_smallest_distances_to_reflective_boundaries = sort(cL - abs(curr_XYZ_colVec));
    d_reflectiveBoundary = three_smallest_distances_to_reflective_boundaries(1); % Distance to current molecule's nearest reflective boundary
    
    
    % Particles can always exceed boundary.
    d_partners_nearestBoundary = min( cL - abs(partners_XYZ_colVecs) ) + ...
        exceedBoundaryDist;
    d_nearestBoundary =  d_reflectiveBoundary + exceedBoundaryDist;
    
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------  
    D = Dmonomer; % Michaelis Menten
    
    % Evaluate the max waiting time for current product to reach its
    % nearest relevant boundary.
    maxAllowed_t_wait_curr = (d_nearestBoundary)^2 / (6*D*n_sigma^2);

    
    % Evaluate the max waiting time for partners to reach their
    % nearest relevant boundary.
    T_LAST_partners = assemblies_lastUpdate_absTime(Cols); % time of last update
    t_offset_partners = t_absolute - T_LAST_partners; % past diffusing time (until now)
    D_partners = kB_T ./ (6*pi*eta*currentStokesRs(Cols));
    maxAllowed_t_wait_partners = (d_partners_nearestBoundary - ....
        6.*D_partners.*n_sigma^2.*t_offset_partners).^2 ./ (6*D*n_sigma^2);

    assert( all(maxAllowed_t_wait_partners > 0) );
    % For each pairing (i.e. product & partner), list the max waiting
    % time.
    maxAllowed_t_wait_list = min(maxAllowed_t_wait_curr,maxAllowed_t_wait_partners);
    inds_partner_sets_maxWait = maxAllowed_t_wait_curr > maxAllowed_t_wait_partners;
    
    % Which pairs of diff spheres intersect by the max allowed waiting time?
    % i.e., R_curr(max_t) + R_partner(max_t) >= d ??
    % First, remember that partners may have been diffusing already, so we
    % want to look at diffusion sphere intersections after *their* ELAPSED
    % diffusion times.
    t_offset_list = zeros(1, length(maxAllowed_t_wait_list));
    t_offset_list(inds_partner_sets_maxWait) = t_offset_partners(inds_partner_sets_maxWait);
    
    % When the partner sets the max waiting time, we need to consider its
    % total elapsed diffusion time in determining if there is diffusion 
    % sphere overlap by the max waiting time.
    maxAllowed_t_elapsed_list = maxAllowed_t_wait_list + t_offset_list; 
    
    % Get all pairs with overlapping diffusion spheres by their max waiting time.
    % Also, to be valid, the wait time for a given pair must be long enough
    % bc we integrate the propensity for a duration at most this wait time.
    % We don't want to integrate the propensity from 0 to 1e-9 for example
    % as that's just a waste of resources.
    valid_pair_indices = n_sigma*sqrt(6.*D.*maxAllowed_t_elapsed_list) + ...
        n_sigma*sqrt(6.*D_partners.*maxAllowed_t_elapsed_list) ...
        >= distList__currRow & maxAllowed_t_wait_list > minPropensityDuration;
    
    % SET HARD GLOBAL UPPER LIMIT ON TIME CURRENT ASSEMBLY CAN DIFFUSE.
    % THIS HARD LIMIT IS INDEPENDENT OF NEAREST BOUNDARIES AND INDEPENDENT
    % OF THE DEFINITION OF NEARBY.
    % THIS HARD GLOBAL UPPER LIMIT IS ALSO APPLIED TO POS-ONLY-UPDATES.
    % THIS HARD GLOBAL UPPER LIMIT ALSO APPLIES TO PROPENSITY INTEGRATION
    % IN computeWaitTime_____v2_alt.m
    valid_pair_indices_HARDLIMIT = maxAllowed_t_wait_list <= T_HARDLIMIT_DIFFUSE;  
    valid_pair_indices = valid_pair_indices & valid_pair_indices_HARDLIMIT;
    
    % Only keep valid entries. 
    maxAllowed_t_wait_list = maxAllowed_t_wait_list(valid_pair_indices);
    
    % Update <Cols> and <distList__currRow> for valid pairs
    Cols = Cols(valid_pair_indices); % i.e. [ 3 4 ]
    distList__currRow = assemblyDistVector(Cols); % i.e. [ 1.2 3.01 ]
    
    % For each inner loop, there are at most 2 event vectors that can be 
    % added to ( <lists1> and <lists2> ). The contents of both will be 
    % added to the PQ after the inner loop.
    numCols = length(Cols);
    lists1 = zeros(numCols, 7);    
    lists2 = zeros(numCols, 7); 
    

    % ********** END NEW CODE ************
    
    % INNER LOOP.
    if numCols > 0
    parfor c_ind = 1:numCols 

    d = distList__currRow(c_ind);    
    col = Cols(c_ind);
    if d < inf % This entry is a valid reaction pair.
        Pk_minus_Tk = log( 1/rand );
        
        T_LASTs = [assemblies_lastUpdate_absTime(row), assemblies_lastUpdate_absTime(col) ];        
        %Ds = [ (kB_T / (6*pi*eta*currentStokesRs(row)) ), (kB_T / (6*pi*eta*currentStokesRs(col)) ) ];
        Ds = [Dmonomer Dmonomer]; % to match Chew et. al
        diffSphereRadii = [ n_sigma*(6*Ds(1)* (t_absolute - T_LASTs(1)) )^(1/2),...
            n_sigma*(6*Ds(2)* (t_absolute - T_LASTs(2)) )^(1/2) ];      
        
% % % %         % Assign particle labels. A -> smaller diffusion sphere, B ->
% % % %         % larger diffusion sphere (assuming sphere's are diff sizes, 
% % % %         % otherwise A is smaller diffusion rate)        
% % % %         if diffSphereRadii(1) == diffSphereRadii(2)
% % % %            [~,DaIndex] = min(Ds);
% % % %            [~,DbIndex] = max(Ds); 
% % % %         else
% % % %             [~,DaIndex] = min( diffSphereRadii );         
% % % %             [~,DbIndex] = max( diffSphereRadii );
% % % %         end
% % % % 
% % % %         if DaIndex == DbIndex; DaIndex=1; DbIndex=2;end
        
        % Assign particle labels. A -> row (just updated particle).
        % B -> partner
        DaIndex=1; DbIndex=2;
        
        Da = Ds(DaIndex);
        Db = Ds(DbIndex);
        T_LASTa = T_LASTs(DaIndex); % time reactant A's state last updated.
        T_LASTb = T_LASTs(DbIndex); % time reactant B's state last updated.
        t_offset_a = t_absolute - T_LASTa; % Elapsed time since A update.
        t_offset_b = t_absolute - T_LASTb; % Elapsed time since B update.

        assert(t_offset_a>=0)
        assert(t_offset_b>=0)
        
        % ------------------------
        % BEGIN Compute wait_time
        % ------------------------
        % Find the appropriate CDF curve in the pre-computed data (the CDF
        % curve at the distance closest to our currently considered <d>.
        [~,d_index] = min( abs(d - distances) );
        begin_ind = (d_index-1)*numVariances+1;
        end_ind = (d_index-1)*numVariances+numVariances;
        IntegratedPropensity_curve = rate_constant * CDF_F_list(begin_ind:end_ind);
        variances = var_list(begin_ind:end_ind);
       
        assert(t_offset_a==0)
        wait_time = computeWaitTime_preComp_CDF(IntegratedPropensity_curve, variances,...
            Da, Db, t_offset_b);
        % ------------------------
        % END Compute wait_time
        % ------------------------
     

        % ADD [WAITTIME+t_absolute, slower diffusing ASSEMBLY, faster diffusing ASSEMBLY, updateType, t_elapsed_a,t_elapsed_b] 
        % TO PRIORITY QUEUE. updateType==1 means reaction&position(s). updateType==2 means                        % 
        % position(s) only. 
        % t_elapsed_a, t_elapsed_b are the durations each particle had been
        % diffusing.
        
        % -----------------------------------------------------------------
        if wait_time > 0 % Add Bimolecular Event to PQ
        % -----------------------------------------------------------------    
            % t_elapsed_x is the total duration the diffusion sphere has been
            % expanding, starting from when x's state was last updated, ending
            % when x undergoes its next reaction. 
            % These variables are needed for correctly sampling the reaction
            % LOCATION. Also needed for position only updates.
            t_elapsed_a = wait_time + t_offset_a; 
            t_elapsed_b = wait_time + t_offset_b;           

            % Add the bimolecular reaction to the PQ.

            % The smaller diffusion (A) sphere will be used as the origin
            % when sampling a reaction location.
            reactant_list = [currentAssemblyIDs(row),currentAssemblyIDs(col)];
            list = [wait_time + t_absolute,reactant_list(DaIndex),reactant_list(DbIndex),1,t_elapsed_a,t_elapsed_b,1];
            lists1(c_ind,:) = list;
                
        end
    end
    % -----------------------------------------------------------------    
    end     % ENDS PARFOR
    % -----------------------------------------------------------------
    % Add Bimolecular Events to PQ.
    [logical_inds, ~]=ismember(lists1,zeros(1,7),'rows');
    if sum(logical_inds)>0; lists1(logical_inds,:) = []; end
    PQueue = [PQueue;lists1];    
    
    end    
  
    %}
    
    % -----------------------------------------------------------------        
    % 2. Add POSITION ONLY UPDATE
    % -----------------------------------------------------------------

    % ---------------------------------------------------------------------    
    t_wait = T_HARDLIMIT_DIFFUSE;
    list = [t_wait + t_absolute,currentAssemblyIDs(row),-1,2,t_wait,-1,-1];
    PQueue(end+1,:) = list;


    % ---------------------------------------------------------------------
    % 3. Unimolecular events involving the new product.
    % ---------------------------------------------------------------------
    
    %    If  greater than threshold duration, temporarily store
    %    position only update for it.
    %     Add any possible UNIMOLECULAR reactions to the PriorityQ.
    %     Simple model of unimolecular reactions: No dependence on
    %     assembly size (number of bonds), so long as dimer or
    %     larger.
    %     When executing a break reaction (while loop), the molecule to break
    %     off is chosen randomly. Obviously this unimolecular model is not
    %     great but I'm going for simple right now. 


    currentAssemID = currentAssemblyIDs( row );
    currentAssemSize = sum(subunitIDsAssemblyIDs == currentAssemID);
    if currentAssemSize >= 2
      
    % Iterate over possible unimolecular reaction types (Michaelis-Menten).
    for rType = 2:3
        % ADD [REACTION TIME, ASSEMBLY to Break, -1, updateType, t_elapsed_a, -1, reactionTpe] TO PRIORITY QUEUE.
        % The slower diffusing assembly will be used as the origin
        % when sampling a reaction location.
        %wait_time = -k_UNI(rType-1)*log(rand(1,1));
        wait_time = (1 / k_UNI(rType-1))*log(1/rand(1,1));
        T_LAST = assemblies_lastUpdate_absTime(row);
        t_offset = t_absolute - T_LAST;
        assert(t_offset==0)

        % The amount of time A has diffused before it undergoes unimolecular reaction
        t_elapsed_a = wait_time;% + t_offset; 

        % UPDATE
        % Unimolecular event does not depend on diff sphere encountering
        % boundary        
        list_UNI = [wait_time + t_absolute, currentAssemID, -1, 1, t_elapsed_a -1, rType];
        PQueue(end+1,:) = list_UNI;   
       
    end    
    end    
    
[logical_inds_PQueue, ~]=ismember(PQueue,zeros(1,7),'rows');
PQueue(logical_inds_PQueue,:) = [];
% assert( sum(logical_inds_PQueue) == 0 )

end





