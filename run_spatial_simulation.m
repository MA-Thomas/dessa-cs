clear all
numReactants = 1000;%[100,200,400,800,1600,3200,6400,12000,24000];

trials = length(numReactants);%4;
bimolecular_count_trials = zeros(1,trials);
unimolecular_count_trials = zeros(1,trials);
positionOnlyUpdate_count_trials = zeros(1,trials);
mainLoop_time_trials = zeros(1,trials);
maxRxDistance_trials = zeros(1,trials);
t_end_trials = zeros(1,trials);
initialization_time_trials = zeros(1,trials);


% Michaelis-Menten
numE_trials = cell(1,trials); % to hold # at the end of each simulation
numS_trials = cell(1,trials);
numES_trials = cell(1,trials);
numP_trials = cell(1,trials);
timePoints_trials = cell(1,trials);
stored_events_trials = cell(1,trials);


%%
for trial = 1:trials
%%    

% % ---------- HOUSEKEEPING -----------------------------------------------
numTotal = numReactants(trial);

containerLength = 90.9;%30; % units [um]. V=30um^3 like Chew et al. diffusion test. 90.9 for reaction test

lower = -45;%-13; %-35
upper = 45;%13;  %35 
allAssembliesXYZ = lower + (upper-lower).*rand(3,numTotal); % from box

unimolecular_count = 0; % count number of times a break reaction happens.
bimolecular_count = 0;
positionOnlyUpdate_count = 0;

StokesEinsteinRadius = 1e-8; % For monomer. Units: m (i.e. 1e-8m is 10nm) 
% ** The resulting inferred D will have units um^2*s^-1 **

% NO LONGER USED IN COMPUTING PROPENSITY ---------
eta = 1e-3; % viscosity of water at room temp is ~ 10e-3 kg *m^-1 *s^-1
kB_T = 4.16e-9; % Units: kg*um^2*s^-2 / Units J = kg*m^2*s^-2 =>  4.16e-21;
% NO LONGER USED IN COMPUTING PROPENSITY ---------


t_absolute = 0;

PQueue = [];

% Initialize stored events vector
stored_events = zeros(1,7);

currentMaxID = numTotal; 

% distanceMatrix = [dist(A1,A1) dist(A1,A2) ... dist(A1,BM)
%                   dist(A2,A1)
%                     ...                                   ]
% Upper triangular part of FULL DISTANCE MATRIX - for all assemblies.
% The size of this matrix never changes, though as the number of assemblies
% in the system decreases, more and more rows/columns will fill with 'inf'
% entries.
assemblyDistMatrix = real(triu( (bsxfun(@plus,dot(allAssembliesXYZ,allAssembliesXYZ,1)',dot(allAssembliesXYZ,allAssembliesXYZ,1))-2*(allAssembliesXYZ'*allAssembliesXYZ) ).^(1/2) ));

% Each assemblyID can appear more than once.
subunitIDsAssemblyIDs = 1:size(assemblyDistMatrix,2); %index=subunitID, element(index)=assemblyID

% Each assemblyID can appear at most once. Remaining entries may be 'inf'.
currentAssemblyIDs = 1:size(assemblyDistMatrix,2); %index=col in assemblyDistMatrix, element(index)=assemblyID or inf

% Keep track of when each assembly was last updated.
% This allows us to compute for how long an assembly has been diffusing.
assemblies_lastUpdate_absTime = zeros(1,length(currentAssemblyIDs));



% For a given reactant pair, at least this much time must be available for
% propensity function integration for the algorithm to call the function
% that samples a waiting time. 
minPropensityDuration = 1e-4; % units [s] 

% Setting this too low can lead to non-convergence in
% Pr_dist2Gaussians_powerSeries_alt.m
earliestPropensityTime = 1e-4;

% -------------------------------------------------------------------------
% Michaelis Menten benchmark parameters BEGIN
% -------------------------------------------------------------------------
% There are 4 species: E,S,ES,P, and Null. I.e. species 1,2,3,4, and -1.
% Null just means the assembly no longer exists in the simulation.
numS = floor(0.9 * numTotal);
numE = numTotal - numS;
currentAssemblyTYPES = [ones(1,numE),2*ones(1,numS)];


numE = zeros(1,1001); % to hold vector of # within single simulation
numS = zeros(1,1001);
numES = zeros(1,1001);
numP = zeros(1,1001);
timePoints = zeros(1,1001);
currIntervalDone = zeros(1,1001);
nCounter = 1;

% Unimolecular rate constants
k_ES_to_E_S = 0.1;
k_ES_to_E_P = 0.1;
k_UNI = [k_ES_to_E_S,k_ES_to_E_P];

% Bimolecular rate constant (intrinsic)
rate_constant = 540;

maxSimTime = 100; % Do not simulate beyond t_absolute = maxSimTime [seconds]


% Assemblies can be roughly this distance beyond boundary before
% reflective boundary conditions are enforced.
exceedBoundaryDist = 0.3 * containerLength;

% Reaction-Rules:
% 1 : E + S -> ES
% 2 : ES -> E + S
% 3 : ES -> E + P

% Dmonomer = kB_T / (6*pi*eta*min(currentStokesRs)); % fastest monomer in the simulation
Dmonomer = 1; % Diffusion Coefficient Michaelis-Menten, Chew et al benchmark
Dh = Dmonomer; % units um^2*s^-1 

% StokesEinsteinRadius = 1e-8; % For monomer. Units: m (i.e. 1e-8m is 10nm) 
StokesEinsteinRadius = kB_T / (6*pi*eta); % This guarantees D = kB_T / (6*pi*eta*StokesEinsteinRadius) = 1
% ** The resulting inferred D will have units um^2*s^-1 **

% The Stokes-Einstein radius of each assembly. 
% Add random noise to these so that none of the inferred diffusion consts
% are identical and all are slower than Dh=Dmonomer .
currentStokesRs = StokesEinsteinRadius*ones(1,length(currentAssemblyIDs)) + ...
    abs( normrnd(0,0.01*StokesEinsteinRadius,[1,length(currentAssemblyIDs)]) ); %index=col in assemblyDistMatrix, element(index)=stokesR or inf

% -------------------------------------------------------------------------
% Michaelis Menten benchmark parameters END
% -------------------------------------------------------------------------


%--------------------------------------------------------------------------
% FIXED PARAMETERS for pre-computing propensityCDFs
n_sigma = 5; % used to calc max diffusion distance corresponding to max diffusion time
Dfast = Dmonomer; 
contact_sphere_radiusSqd = .2;%3; % units [um]
% This has to be small (0.8 seems to work, 3 sometimes doesn't) for the
% function that approximates the propensity function to be well-behaived.
% (Pr_dist2Gaussians_powerEries_alt.m)


% SETS THE UPPER LIMIT ON HOW LONG AN ASSEMBLY CAN DIFFUE BEFORE ITS NEXT
% EVENT. (Other factors like proximity to a boundary, or only considering "nearby"
% partners can lower this upper limit).
T_HARDLIMIT_DIFFUSE = 16; % units [s]
Da = Dmonomer; Db = Dmonomer; % units [um^2 / s]
maxRxDist_T_HARD = n_sigma*sqrt(6*Dfast*T_HARDLIMIT_DIFFUSE); % units [um]
% FIXED PARAMETERS
tic;
if trial == 1
[CDF_F_list,distances,numVariances,var_list] = preCompute_propensityCDF_curves(...
    contact_sphere_radiusSqd,T_HARDLIMIT_DIFFUSE,maxRxDist_T_HARD,Da,Db,...
    earliestPropensityTime);
end
preCompuTime = toc

% only sample potential reactions from reactants separated by less than
% maxRxDistance.
maxRxDistance = maxRxDist_T_HARD; %(1/8)*containerLength;%(1/6)*containerLength; % This should be low (e.g. 1/7) for fast initialization


tic;
% % ---------- INITIALIZATION ---------------------------------------------
% Helper Function replaces (/contains) FOR loop.
% Based on initial positions of assemblies and other factors defined above,
% add all possible reactions to the queue (PQ) and return PQ.
% The while loop below is the main part of the algorithm. It will execute
% the events in the queue, and add new possible events to PQ.
[ PQueue,initial_events,lists_cell ] = initializeRxQueue______v2_alt( assemblyDistMatrix,...
    currentAssemblyIDs,...
    PQueue,t_absolute,assemblies_lastUpdate_absTime,allAssembliesXYZ,...
    maxRxDistance,n_sigma,rate_constant,Dh,...
    containerLength,currentAssemblyTYPES,...
    exceedBoundaryDist,minPropensityDuration,T_HARDLIMIT_DIFFUSE,...
    var_list,numVariances,distances,CDF_F_list);

%--------------------------------------------------------------------------
initialization_time = toc
display('Initialization Complete')
pause(2)

reaction_runtimes = 0;

XYZ_storage_POU1 = zeros(100000,3);
xyz_counter = 1;

mainLoop_startTime = tic;
%%
% ---------------------- Start algorithm loop -----------------------------

continue_algo = 1;
while_counter = 1;

while continue_algo
    
    % % GET NEXT EVENT FROM QUEUE
    [~, ind_t] = min(PQueue(:,1));
    event_vector = PQueue(ind_t,:);
    PQueue(ind_t,:) = [];
    
    t_reaction = event_vector(1);
    a1 = event_vector(2); % This is the assemblyID a1 (smaller diffusion sphere at t_reaction if bimolecular)
    a2 = event_vector(3); % This is the assemblyID a2 (larger diffusion sphere at t_reaction if bimolecular)
    updateType = event_vector(4); % 1 means reaction&postion, 2 means only position
    t_elapsed_a = event_vector(5); % time elapsed since A started diffusing
    t_elapsed_b = event_vector(6); % time elapsed since B started diffusing
    reactionRule = event_vector(7); % which (of possibly many) reaction rules to apply within bi/uni update script
    
    assert(a1 ~= a2);
    

    if updateType == 1
        
        row_col_logicals_index_a1 = currentAssemblyIDs==a1; % This is the row/col in the current distance matrix for assembly a1
        row_col_logicals_index_a2 = currentAssemblyIDs==a2; % This is the row/col in the current distance matrix for assembly a2

        % If reaction no longer possible, continue to next loop iteration.
        % E.g. if a previous reaction used up one of the assemblies of the
        % current reaction.
        if sum(row_col_logicals_index_a1) == 0

            % % If no events left, end.
            if isempty(PQueue)
                continue_algo = 0;
            end            
            continue
        end % assembly a1 no longer exists
        
        if a2>0 && sum(row_col_logicals_index_a2) == 0
            
            % % If no events left, end.
            if isempty(PQueue)
                continue_algo = 0;
            end
            continue
        end % bimolecular reaction no longer possible, a2 no longer exists
    
        %reaction_startTime = tic;
        % % ---------------------------- % %
        if a2<0 
            %disp('run uni script')
            script__unimolecular_update______v2_alt
        else
            %disp('run bi script')
            script__bimolecular_update______v2_alt
        end
        % % ---------------------------- % %
        numE(nCounter) = sum(currentAssemblyTYPES==1);
        numS(nCounter) = sum(currentAssemblyTYPES==2);
        numES(nCounter) = sum(currentAssemblyTYPES==3);
        numP(nCounter) = sum(currentAssemblyTYPES==4);
        timePoints(nCounter) = t_reaction;
        nCounter = nCounter+1;
        stored_events(end+1,:) = event_vector;
        %reaction_runtimes(end+1) = toc(reaction_startTime);
        
    elseif updateType == 2
        
        % <event_Matrix> is the set of all rows in PQueue with same
        % t_reaction
        event_Matrix = event_vector;
        [logical_contains,~] = ismember(PQueue(:,1),event_vector(1));   
        if sum(logical_contains) > 0
            event_Matrix = [event_vector; PQueue(logical_contains,:)]; 
        end
        numParallel = 1; % default
        
        if size(event_Matrix,1) > 0
            % Keep in <event_Matrix> only POU events for still existing
            % assemblies.
            [logical_existing_assemblies,~] = ismember(event_Matrix(:,2),...
                currentAssemblyIDs);            

            % It may be the case that a bimol event and a POU event occur
            % at the same t_reaction. Only keep the POUs.
            logical_existing_POUs = event_Matrix(:,4) == 2;
            
            logical_inds = logical_existing_assemblies & logical_existing_POUs;            
            event_Matrix = event_Matrix(logical_inds,:);            
            
            % All events are POUs (updateType==2)
            assert(all(event_Matrix(:,4) == 2))

            % Delete event_Matrix rows from PQueue.
            [logical_toDelete,~] = ismember(PQueue, event_Matrix,'rows');
            PQueue = PQueue(~logical_toDelete,:);
            numParallel = size(event_Matrix,1);
        end

        if numParallel > 0    
         
        % % ------------------- % %
        Part1_script__positionOnly_update______v2_alt % Parallelized executions
        Part2_script__positionOnly_update______v2_alt % Parallelized sampling new events
        % % ------------------- % %
        
        stored_events = [stored_events; event_Matrix];
        
        end
        
        
    end
    
    
    while_counter = while_counter + 1;


    % % If no events left, end.
    if isempty(PQueue)
        continue_algo = 0;
    end
    
    if t_absolute >= maxSimTime
        continue_algo = 0;
    end


end

bimolecular_count
unimolecular_count
positionOnlyUpdate_count
mainLoop_time = toc(mainLoop_startTime)
t_absolute

bimolecular_count_trials(trial) = bimolecular_count;
unimolecular_count_trials(trial) = unimolecular_count;
positionOnlyUpdate_count_trials(trial) = positionOnlyUpdate_count;
mainLoop_time_trials(trial) = mainLoop_time;
maxRxDistance_trials(trial) = maxRxDistance;
t_end_trials(trial) = t_absolute;
initialization_time_trials(trial) = initialization_time;

stored_events_trials{trial} = stored_events;

numE_trials{trial} = numE;
numS_trials{trial} = numS;
numES_trials{trial} = numES;
numP_trials{trial} = numP;
timePoints_trials{trial} = timePoints;
clearvars -except bimolecular_count_trials unimolecular_count_trials ...
    positionOnlyUpdate_count_trials mainLoop_time_trials maxRxDistance_trials trials trial t_end_trials ...
    initialization_time_trials numE_trials numS_trials numES_trials numP_trials ...
    stored_events_trials reaction_runtimes d_inner1 d_inner2 rate_constant ...
    timePoints_trials numReactants T_HARD_Diffuse_Limits XYZ_storage_POU1 ...
    CDF_F_list distances numVariances var_list


end



        