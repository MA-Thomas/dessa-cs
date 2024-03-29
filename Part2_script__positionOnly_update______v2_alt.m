% POSITION ONLY UPDATE - Parallel Version

%{
When a Pos-Only-Update (POU) appears in the PQ, keep removing next events
as long as they share the same reaction time. Execute these reactions in
parallel. 

Data representing significant parallelization overhead:
allAssembliesXYZ
assemblyDistMatrix
assemblies_lastUpdate_absTime

Next - after parallel event execution - select new events (normal
parallelized loop with outerloop over just-updated assemblies.

%}






% % ----------------------------------------------------
% % Update the Priority Queue with potential new events.
% % ----------------------------------------------------
% 1. Bimolecular event involving the assembly and 1 other assembly. 
%    If none, temporarily store position only updates info for this channel.

% Let rows cover the single new product ('new' means 'updated state').
% Let cols cover all other assemblies.
row_products = appropriateCols; % Computed in Part1_script_positionOnly_update______v2_alt.m

numProducts = length(row_products);
lists_POU = zeros(numProducts,7);
lP = 1; % for indexing into list_POU

for row_number = 1:numProducts
    row = row_products(row_number);
    
    
    % Compute list of distances to all other molelcules
    assemblyDistVector = real( (bsxfun(@plus,dot(allAssembliesXYZ,allAssembliesXYZ,1)',...
        dot(allAssembliesXYZ(:,row),allAssembliesXYZ(:,row),1))-...
        2*(allAssembliesXYZ'*allAssembliesXYZ(:,row)) ).^(1/2) )';
    
    % If this reactant is not allowed to react with any others, skip. 
    if ~any(assemblyDistVector < inf);continue;end 

    
    columns = 1:length(assemblyDistVector);
    
    % Only consider columns with:
    % 1. distances less than <maxRxDistance>
    % 2. distances greater than <maxRxDistance> only if partner has a
    %    nonzero t_offset and d < maxRxDistance+R(t_offset)
    %    i.e. partner has been diffusing long enough to reach distance
    %    of <maxRxDistance>.
        % ---------------------------------------------------------------------
        % ----------------------Uncommented below-------------------------------
    col_indices1 = assemblyDistVector(columns) <= maxRxDistance;
     
    Rs = real( n_sigma*(6*(kB_T ./ (6*pi*eta*currentStokesRs(columns)) )*...
        t_absolute - assemblies_lastUpdate_absTime(columns)).^(1/2) );    
    
    col_indices2 = assemblyDistVector(columns) <= maxRxDistance+Rs;
        % ---------------------------------------------------------------------
        % ----------------------Uncommented above-------------------------------
    
%     col_indices1 = logical(columns); % start with all columns
%     col_indices2 = logical(columns); % start with all columns
    col_indices = col_indices1 | col_indices2; 
    
    % Don't over count. I.e. don't consider the same pairs in successive
    % iterations of the for loop.
    col_indices( row_products(1:row_number) ) = false;

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
    
   

    cL = (containerLength/2);
    
    Cols = columns(col_indices); % i.e. [ 3 4 11 13 ]

    partners_XYZ_colVecs =  allAssembliesXYZ(:,Cols); 
    bad_indices = find(all(isnan(partners_XYZ_colVecs),1));
    partners_XYZ_colVecs( :,bad_indices) = [];
    Cols(bad_indices) = [];
    
    distList__currRow = assemblyDistVector(Cols);
    
    
    % Calculate distance to the nearest boundary.
    
    curr_XYZ_colVec = allAssembliesXYZ(:,row);
    assert( sum( abs(curr_XYZ_colVec) > cL*ones(size(curr_XYZ_colVec)) ) ==0)    

    three_smallest_distances_to_reflective_boundaries = sort(cL - abs(curr_XYZ_colVec));
    d_reflectiveBoundary = three_smallest_distances_to_reflective_boundaries(1); % Distance to current molecule's nearest reflective boundary
    
  
    % Particles can always exceed boundary. no need for d_inner1/d_inner2
    d_partners_nearestBoundary = min( cL - abs(partners_XYZ_colVecs) ) + ...
        exceedBoundaryDist;
    d_nearestBoundary =  d_reflectiveBoundary + exceedBoundaryDist;
    
    D = Dmonomer; % Michaelis-Menten
    
    % Evaluate the max waiting time for current product to reach its
    % nearest relevant boundary.
    maxAllowed_t_wait_curr = (d_nearestBoundary)^2 / (6*D*n_sigma^2);

    
    % Evaluate the max waiting time for partners to reach their
    % nearest relevant boundary.
    T_LAST_partners = assemblies_lastUpdate_absTime(Cols); % time of last update
    t_offset_partners = t_absolute - T_LAST_partners; % past diffusing time (until now)
    D_partners = ones(1,length(Cols));%kB_T ./ (6*pi*eta*currentStokesRs(Cols));
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
    
    % This will be an input to computeWaitTime function.
    % Only keep valid entries. 
    maxAllowed_t_wait_list = maxAllowed_t_wait_list(valid_pair_indices);    
    
    % Update <Cols> and <distList__currRow> for valid pairs
    Cols = Cols(valid_pair_indices); % i.e. [ 3 4 ]
    
    distList__currRow = assemblyDistVector(Cols);
    
    numCols = length(Cols);
    lists1 = zeros(numCols, 7); 
    
    
    % INNER LOOP.
    % Iterate over all valid pairs of reactant molecules.
    if numCols > 0
%     parfor c_ind = 1:numCols
    for c_ind = 1:numCols

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

        assert(t_offset_a==0)
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
     

        % -----------------------------------------------------------------
        if wait_time > 0 % Add Bimolecular Event to PQ
        % -----------------------------------------------------------------    
            % t_elapsed_x is the total duration the diffusion sphere has been
            % expanding, starting from when x's state was last updated, ending
            % when x undergoes its next reaction. 
            t_elapsed_a = wait_time + t_offset_a; 
            t_elapsed_b = wait_time + t_offset_b;

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
    
    
    % -----------------------------------------------------------------        
    % 2. Add POSITION ONLY UPDATE
    % -----------------------------------------------------------------
    t_wait = min([maxAllowed_t_wait_curr,T_HARDLIMIT_DIFFUSE]);
    list = [t_wait + t_absolute,currentAssemblyIDs(row),-1,2,t_wait,-1,-1];
    lists_POU(lP,:) = list;
    lP = lP + 1;        

end

PQueue = [PQueue; lists_POU];

[logical_inds_PQueue, ~]=ismember(PQueue,zeros(1,7),'rows');
PQueue(logical_inds_PQueue,:) = [];
% assert( sum(logical_inds_PQueue) == 0 )

% UPDATE: DON'T ADD NEW UNIMOLECULAR EVENTS AFTER EXECUTING
%         POSITION-ONLY-UPDATES.
