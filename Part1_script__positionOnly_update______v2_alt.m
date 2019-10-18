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
positionOnlyUpdate_count = positionOnlyUpdate_count + numParallel;
assert( size(event_Matrix,1) == numParallel )

new_locations = zeros(3,numParallel);
appropriateCols = zeros(1,numParallel);

% parfor pf = 1:numParallel
for pf = 1:numParallel
    
event_vector = event_Matrix(pf,:);
t_reaction = event_vector(1);
a1 = event_vector(2); % This is the assemblyID a1 (smaller diffusion sphere at t_reaction if bimolecular)
a2 = event_vector(3); % This is the assemblyID a2 (larger diffusion sphere at t_reaction if bimolecular)
updateType = event_vector(4); % 1 means reaction&postion, 2 means only position
t_elapsed_a = event_vector(5); % time elapsed since A started diffusing
t_elapsed_b = event_vector(6); % time elapsed since B started diffusing
reactionRule = event_vector(7); % which (of possibly many) reaction rules to apply within bi/uni update script

row_col_logicals_index_a1 = currentAssemblyIDs==a1;
appropriateCols(pf) = find(row_col_logicals_index_a1);

Da = (kB_T / (6*pi*eta*currentStokesRs( find(row_col_logicals_index_a1) )) );    
Da = Dmonomer; % to match Chew et al.

% [ rA, thetaA, ~,~ ] = sample_Rx_location_unif_q( wait_time,Da,Db,d,n_sigma ); % caseNum calc within the function
rA = sample_POS_ONLY_update( t_elapsed_a,Da,n_sigma ); 

% % -----------------------------------------------------------------------
% % PARTICLE A ------------------------------------------------------------
% % Convert, using random angular variables, to cartesian displacements
% we want sph2cart to produce both positive and negative z values (at
% random)
if rand(1,1)>0.5; rA = -rA;end 

[xA,yA,zA] = sph2cart(2*pi*rand,pi*rand,rA);
XYZ_storage_POU1(xyz_counter,:) = [xA,yA,zA];
xyz_counter = xyz_counter+1;

% For use in applying reflective BC if necessary.
xyz_Global_A = allAssembliesXYZ(:,row_col_logicals_index_a1);

% Add displacement vector to original position vector.
xyz_Global = allAssembliesXYZ(:,row_col_logicals_index_a1) + [xA,yA,zA]';  

% Apply Reflective Boundary Conditions to ensure xyz_Global is in
% simulation box
if any(abs(xyz_Global) > containerLength/2)
  
    %xyz_Global = testing_reflective_boundary_conditions(xyz_Global,containerLength);
    xyz_Global = testing_reflective_boundary_conditions(xyz_Global,xyz_Global_A,containerLength);
end

% Store new positions so we can update <allAssembliesXYZ> outside the 
% parfor loop because MATLAB won't allow it within parfor. 
new_locations(:,pf) = xyz_Global; 

% % END PARTICLE A --------------------------------------------------------
% % -----------------------------------------------------------------------


end
    
% % % t_reaction = event_Matrix(1,1);
% % % a1_list = event_Matrix(:,2); % This is the assemblyID a1 (smaller diffusion sphere at t_reaction if bimolecular)
% % % t_elapsed_a_list = event_Matrix(:,5); % time elapsed since A started diffusing
% % % 
% % % % Locb is the size of <currentAssemblyIDs>. It's entries are indices to
% % % % a1_list.
% % % [Lia,Locb] = ismember(currentAssemblyIDs,a1_list);
% % % Da_list = kB_T / (6*pi*eta*currentStokesRs( Lia) );  
% % % 
% % % 
% % % % We want event_Matrix rows to be in the same order as the assemblies
% % % % appear in currentAssemblyIDs. This allows us to use <Lia> to update
% % % % the correct entries in variables below.
% % % eventMatrix_ordering = Locb(Locb > 0);
% % % assert(length(eventMatrix_ordering) == size(event_Matrix,1));
% % % event_Matrix = event_Matrix(eventMatrix_ordering,:);
% % % a1_list = event_Matrix(:,2); % This is the assemblyID a1 (smaller diffusion sphere at t_reaction if bimolecular)
% % % t_elapsed_a_list = event_Matrix(:,5); % time elapsed since A started diffusing
% % % 
% % % 
% % % rA_colVec = sample_POS_ONLY_update( t_elapsed_a_list,Da_list,n_sigma ); 
% % % 
% % % % % -----------------------------------------------------------------------
% % % % % BEGIN MORE AND REFLECT ------------------------------------------------
% % % % % Convert, using random angular variables, to cartesian displacements
% % % [xA_colVec,yA_colVec,zA_colVec] = sph2cart(2*pi*rand(size(rA_colVec)),...
% % %     pi*rand(size(rA_colVec)),rA_colVec);
% % % 
% % % % For use in applying reflective BC if necessary.
% % % xyz_Global_pre_POU = allAssembliesXYZ(:,Lia);
% % % 
% % % % Add displacement vector to original position vector.
% % % xyz_Global_post_POU = allAssembliesXYZ(:,Lia) + [xA_colVec,yA_colVec,zA_colVec]';  
% % % 
% % % % Determine which cols (which assemblies) of xyz_Global need reflection.
% % % if any(abs(xyz_Global_post_POU) > containerLength/2)
% % %     
% % %     colsToReflect_logical = sum (abs(xyz_Global_post_POU) > containerLength/2, 1);
% % %     colsToReflect = find(colsToReflect_logical);
% % %     % Apply Reflective Boundary Conditions to ensure xyz_Global is in
% % %     % simulation box
% % %     for i = 1:colsToReflect
% % %     xyz_Global_post_POU(:,i) = testing_reflective_boundary_conditions(...
% % %         xyz_Global_post_POU(:,i),...
% % %         xyz_Global_pre_POU(:,i),...
% % %         containerLength);
% % %     end
% % % end
% % % 
% % % new_locations = xyz_Global_post_POU; 
% % % appropriateCols = find(Lia);
% % % % % END MOVE AND REFLECT --------------------------------------------------------
% % % % % -----------------------------------------------------------------------




% Update absolute time, new value used to update Tk.
t_absolute = event_Matrix(1,1); 

% Store new positions in the appropriate columns of <allAssembliesXYZ>.
allAssembliesXYZ(:,appropriateCols) = new_locations;


% Keep track of the absolute time these states were last updated.
assemblies_lastUpdate_absTime(appropriateCols) = t_absolute; % updated
