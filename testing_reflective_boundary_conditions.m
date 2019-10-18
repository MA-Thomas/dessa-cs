function xyz_Global = testing_reflective_boundary_conditions(xyz_Global,...
    xyz_Global_A,containerLength)
% Appropriate for cube simulation volume and Pos-Only or Unimol Update
% xyz_Global_A is the location before POU or unimolecular update.
% xyz_Global is the position after event execution (before reflective BC)

% xyz_Global is a column vector


dist_AtoProduct = norm(xyz_Global - xyz_Global_A);

xG = xyz_Global(1); yG = xyz_Global(2); zG = xyz_Global(3);

xG_A = xyz_Global_A(1); yG_A = xyz_Global_A(2); zG_A = xyz_Global_A(3);

xG_original = xyz_Global(1); yG_original = xyz_Global(2); zG_original = xyz_Global(3);

cL = containerLength/2;

% Which boundaries are candidates for the boundaries A and B first exited.
beyond_cL_logical = xyz_Global;
beyond_cL_logical(beyond_cL_logical>=-cL & ...
    beyond_cL_logical<=cL ) = 0;
beyond_cL_logical(beyond_cL_logical<-cL) = -1;
beyond_cL_logical(beyond_cL_logical>cL) = 1;

beyond_cL_indices = find(beyond_cL_logical);
beyond_cL = beyond_cL_logical(beyond_cL_indices).*beyond_cL_indices;
% e.g. beyond_cL = [1 -3] means the point is beyond cL in
% the x direction and beyond -cL in the z direction

% Determine the distance from A to each potential boundary, and also
% determine the point at which AP intersects each potential boundary.
d_Ia_list = zeros(1,length(beyond_cL));
Ia_list = zeros(3,length(beyond_cL));
for j = 1:length(beyond_cL)
    
    if beyond_cL(j) == 3 %zG > cL
        % Top
        del = cL - zG_A; % delZ is positive
        DEL = zG_original - zG_A; % DELZ is positive

    elseif beyond_cL(j) == -3 %-zG > cL
        % Bottom   
        del = -cL - zG_A; % delZ is negative
        DEL = zG_original - zG_A; % DELZ is negative

    elseif beyond_cL(j) == -2 %-yG > cL
        % Left
        del = -cL - yG_A; % delZ is negative
        DEL = yG_original - yG_A; % DELZ is negative

    elseif beyond_cL(j) == 2 %yG > cL
        % Right
        del = cL - yG_A; % delZ is positive
        DEL = yG_original - yG_A; % DELZ is positive

    elseif beyond_cL(j) == -1 %-xG > cL
        % Back
        del = -cL - xG_A; % delZ is negative
        DEL = xG_original - xG_A; % DELZ is negative

    elseif beyond_cL(j) == 1 %xG > cL
        % Front
        del = cL - xG_A; % delZ is positive
        DEL = xG_original - xG_A; % DELZ is positive  
    end    

    d_Ia_list(j) = dist_AtoProduct* (del/DEL); % distance from A to Ia 

    normalVec = ([xG_original;yG_original;zG_original] - [xG_A;yG_A;zG_A]) /... 
    norm([xG_original;yG_original;zG_original] - [xG_A;yG_A;zG_A]);

    Ia_list(:,j) = [xG_A;yG_A;zG_A] + d_Ia_list(j)*normalVec;

end

[distToFirstBoundaryCrossed_A,indexA] = min(d_Ia_list);
Ia = Ia_list(:,indexA);
firstBoundaryCrossed_A = beyond_cL(indexA);


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%            Preliminary Calculations Done, Start Reflecting
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


    
product_still_outside_boundary = 1;
% ------------ BEGIN While Loop ---------------------------------------
while product_still_outside_boundary

if firstBoundaryCrossed_A == 3 %zG > cL
    % Top
    %display('overshoots z')
    overshoot = zG - cL;
    zG = cL - overshoot;

elseif firstBoundaryCrossed_A == -3 %-zG > cL
    % Bottom
    %display('overshoots -z')
    overshoot = abs(zG) - cL;
    zG = -cL + overshoot;

elseif firstBoundaryCrossed_A == -2 %-yG > cL
    % Left
    %display('overshoots -y')
    overshoot = abs(yG) - cL;
    yG = -cL + overshoot;

elseif firstBoundaryCrossed_A == 2 %yG > cL
    % Right
    %display('overshoots y')
    overshoot = yG - cL;
    yG = cL - overshoot; 

elseif firstBoundaryCrossed_A == -1 %-xG > cL
    % Back
    %display('overshoots -x')
    overshoot = abs(xG) - cL;
    xG = -cL + overshoot;

elseif firstBoundaryCrossed_A == 1 %xG > cL
    % Front
    %display('overshoots x')
    overshoot = xG - cL;
    xG = cL - overshoot;    
end    
xyz_Global = [xG;yG;zG];

% % If the reflection of xyz_Global results in a new position outside
% % the simulation boundaries:
boundaryCrossed = cL < abs(xyz_Global);
if sum(boundaryCrossed) > 0    

    % % Set xyz_Global_A <-- Ia
    xyz_Global_A = Ia;

    xG = xyz_Global(1); yG = xyz_Global(2); zG = xyz_Global(3);
    xG_A = xyz_Global_A(1); yG_A = xyz_Global_A(2); zG_A = xyz_Global_A(3);

    % % Use xyz_Global_A and new product position, xyz_Global, to
    % % determine through which boundary the line
    % % [xyz_Global_A,xyz_Global] intersects first.
    % % Reflect xyz_Global about this boundary.        

    beyond_cL_logical = xyz_Global;
    beyond_cL_logical(beyond_cL_logical>=-cL & ...
        beyond_cL_logical<=cL ) = 0;
    beyond_cL_logical(beyond_cL_logical<-cL) = -1;
    beyond_cL_logical(beyond_cL_logical>cL) = 1;

    beyond_cL_indices = find(beyond_cL_logical);
    beyond_cL = beyond_cL_logical(beyond_cL_indices).*beyond_cL_indices;
    % e.g. beyond_cL = [1 -3] means the point is beyond cL in
    % the x direction and beyond -cL in the z direction

    % Determine the distance from A to each potential boundary, and also
    % determine the point at which AP intersects each potential boundary.
    d_Ia_list = zeros(1,length(beyond_cL));
    Ia_list = zeros(3,length(beyond_cL));
    for j = 1:length(beyond_cL)
        if beyond_cL(j) == 3 %zG > cL
            % Top
            del = cL - zG_A; % delZ is positive
            DEL = zG - zG_A; % DELZ is positive

        elseif beyond_cL(j) == -3 %-zG > cL
            % Bottom   
            del = -cL - zG_A; % delZ is negative
            DEL = zG - zG_A; % DELZ is negative

        elseif beyond_cL(j) == -2 %-yG > cL
            % Left
            del = -cL - yG_A; % delZ is negative
            DEL = yG - yG_A; % DELZ is negative

        elseif beyond_cL(j) == 2 %yG > cL
            % Right
            del = cL - yG_A; % delZ is positive
            DEL = yG - yG_A; % DELZ is positive

        elseif beyond_cL(j) == -1 %-xG > cL
            % Back
            del = -cL - xG_A; % delZ is negative
            DEL = xG - xG_A; % DELZ is negative

        elseif beyond_cL(j) == 1 %xG > cL
            % Front
            del = cL - xG_A; % delZ is positive
            DEL = xG - xG_A; % DELZ is positive  
        end
        d_Ia_list(j) = dist_AtoProduct* (del/DEL); % distance from A to Ia
        normalVec = ([xG;yG;zG] - [xG_A;yG_A;zG_A]) /... 
        norm([xG;yG;zG] - [xG_A;yG_A;zG_A]);
        Ia_list(:,j) = [xG_A;yG_A;zG_A] + d_Ia_list(j)*normalVec;
    end        
    [distToFirstBoundaryCrossed_A,indexA] = min(d_Ia_list);
    Ia = Ia_list(:,indexA);
    firstBoundaryCrossed_A = beyond_cL(indexA);



else
    product_still_outside_boundary = 0; % We are done reflecting.
end

% ------------ END While Loop -----------------------------------------
end
% ------------ END While Loop -----------------------------------------

% Assert that xyz_Global is now in the main simulation box. 
assert( sum( abs(xyz_Global) > cL*ones(size(xyz_Global)) ) ==0) 
    

















%{

Particles are allowed to interact with the nearest (a single) reflective
boundary. The simulation box has 6 reflective sides: 
Top,Bottom,Left,Right,Back,Front

%%}


if zG > cL
    % Top
    overshoot = zG - cL;
    zG = cL - overshoot;
 
elseif -zG > cL
    % Bottom   
    overshoot = abs(zG) - cL;
    zG = -cL + overshoot;
    
elseif -yG > cL
    % Left
    overshoot = abs(yG) - cL;
    yG = -cL + overshoot;
    
elseif yG > cL
    % Right
    overshoot = yG - cL;
    yG = cL - overshoot; 
    
elseif -xG > cL
    % Back
    overshoot = abs(xG) - cL;
    xG = -cL + overshoot;
    
elseif xG > cL
    % Front
    overshoot = xG - cL;
    xG = cL - overshoot;    
end

xyz_Global = [xG;yG;zG];

% Assert that xyz_Global is now in the main simulation box. 
assert( sum( abs(xyz_Global) > cL*ones(size(xyz_Global)) ) ==0) 


%}
end