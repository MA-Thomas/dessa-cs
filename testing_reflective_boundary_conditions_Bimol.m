function xyz_Global = testing_reflective_boundary_conditions_Bimol(...
    xyz_Global_A,xyz_Global_B,xyz_Global,containerLength)
% Appropriate for cube simulation volume and Bimolecular Update
% xyz_Global is a column vector for bimol reaction product
% xyz_Global_A is col vec for first reactant - pre reaction
% xyz_Global_B is col vec for second reactant - pre reaction


%{
% General Principle for performing reflection:

Each reactant has a straight line path to the product location which
crosses the reflective boundary once. Call these crossing points Ia and Ib.
The line IaIb is the axis about which we reflect the product back into the
simulation volume. 

Post reflection, all straight line paths must exist within the boundaries,
i.e., they cannot cross a boundary.

With a curved boundary, this reflection axis needs to be calculated.
With a linear boundary, this reflection axis might be just the boundary
itself or it might need to be calculated.


This file is for a linear boundary:
CASE 1: Lines [Ia,product] and [Ib,product] both cross the same face of the
cube. 
In this case, we simply need to reflect line [Ia,product] about the linear
boundary - i.e. the cube face. (Or the line [Ib,product])

If the product itself is reflected to a location outside the boundary
(possible at corner of the cube, for example), reflect it again.

CASE 2: Lines [Ia,product] and [Ib,product] cross different faces of the
cube.
In this case we need to calculate the reflection axis as the line IaIb.
Reflection ~ rotation about this axis by 180deg. Use the Rodrigues formula.


%}
dist_AtoProduct = norm(xyz_Global - xyz_Global_A);
dist_BtoProduct = norm(xyz_Global - xyz_Global_B);

xG = xyz_Global(1); yG = xyz_Global(2); zG = xyz_Global(3);

xG_A = xyz_Global_A(1); yG_A = xyz_Global_A(2); zG_A = xyz_Global_A(3);
xG_B = xyz_Global_B(1); yG_B = xyz_Global_B(2); zG_B = xyz_Global_B(3);
xG_original = xyz_Global(1); yG_original = xyz_Global(2); zG_original = xyz_Global(3);

cL = containerLength/2;

% % % Assert that at most one coordinate of xyz_Global is outside reflective
% % % boundary
% % assert( sum( abs(xyz_Global) > cL*ones(size(xyz_Global)) ) <= 1)


% % % % Along which dimension(s) was a boundary crossed
% % % boundaryCrossed_logical = abs(xyz_Global) > cL;

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
% -------------------------------------------------------------------------
% Determine the distance from B to each potential boundary, and also
% determine the point at which BP intersects each potential boundary.
d_Ib_list = zeros(1,length(beyond_cL));
Ib_list = zeros(3,length(beyond_cL));
for j = 1:length(beyond_cL)
    
    % Determine through which boundary the line AP intersects first.
    if beyond_cL(j) == 3 %zG > cL
        % Top
        del = cL - zG_B; % delZ is positive
        DEL = zG_original - zG_B; % DELZ is positive

    elseif beyond_cL(j) == -3 %-zG > cL
        % Bottom   
        del = -cL - zG_B; % delZ is negative
        DEL = zG_original - zG_B; % DELZ is negative

    elseif beyond_cL(j) == -2 %-yG > cL
        % Left
        del = -cL - yG_B; % delZ is negative
        DEL = yG_original - yG_B; % DELZ is negative

    elseif beyond_cL(j) == 2 %yG > cL
        % Right
        del = cL - yG_B; % delZ is positive
        DEL = yG_original - yG_B; % DELZ is positive

    elseif beyond_cL(j) == -1 %-xG > cL
        % Back
        del = -cL - xG_B; % delZ is negative
        DEL = xG_original - xG_B; % DELZ is negative

    elseif beyond_cL(j) == 1 %xG > cL
        % Front
        del = cL - xG_B; % delZ is positive
        DEL = xG_original - xG_B; % DELZ is positive  
    end    

    d_Ib_list(j) = dist_BtoProduct* (del/DEL); % distance from A to Ia 

    normalVec = ([xG_original;yG_original;zG_original] - [xG_B;yG_B;zG_B]) /... 
    norm([xG_original;yG_original;zG_original] - [xG_B;yG_B;zG_B]);

    Ib_list(:,j) = [xG_B;yG_B;zG_B] + d_Ib_list(j)*normalVec;

end

[distToFirstBoundaryCrossed_A,indexA] = min(d_Ia_list);
Ia = Ia_list(:,indexA);
firstBoundaryCrossed_A = beyond_cL(indexA);

[distToFirstBoundaryCrossed_B,indexB] = min(d_Ib_list);
Ib = Ib_list(:,indexB);
firstBoundaryCrossed_B = beyond_cL(indexB);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%            Preliminary Calculations Done, Start Reflecting
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


%%{

% Case 1: We can simply reflect straight line path of A 
% (don't worry about B)
if firstBoundaryCrossed_A == firstBoundaryCrossed_B
    
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
        %xG_B = xyz_Global_B(1); yG_B = xyz_Global_B(2); zG_B = xyz_Global_B(3); % not needed here
        
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
    
% Case 2    
else 
    
    product_still_outside_boundary = 1;
    while product_still_outside_boundary
    %display('else case')
    
    % Axis of reflection/rotation given by the 2 points Ia and Ib.
    % Use the Rodrigues formula to rotate vector [Ia,P] about the
    % reflection axis by 180deg.
    % (Equivalently, we could rotate vector [Ib,P] about the same 
    % reflection axis).
    
    % The Ridgrigues formula assume the vector to be rotated and 
    % unit vector k are rooted at the origin. We therefore need to first
    % translate our space such that Ia is the origin.
    % Then perform the rotation. Then translate back.
    Ia_translated = Ia - Ia;
    Ib_translated = Ib - Ia;
    
    % We want to reflect the portion of A's path to xyz_Global that is
    % beyond Ia. Equivalently, reflect portion of B's path beyond Ib.
    pre_refl_vector_translated = xyz_Global - Ia;     

    % unit vector defining axis of rotation
    k = (Ib_translated - Ia_translated) / norm(Ib_translated - Ia_translated); 
    angl_rad = pi;
    pre_refl_vector_translated = pre_refl_vector_translated*cos(angl_rad) + ...
        cross(k,pre_refl_vector_translated)*sin(angl_rad)...
        + k*dot(k,pre_refl_vector_translated)*(1-cos(angl_rad));
    
    % Translate back.
    xyz_Global = pre_refl_vector_translated + Ia;

    % % If the reflection of xyz_Global results in a new position outside
    % % the simulation boundaries:
    boundaryCrossed = cL < abs(xyz_Global);
    if sum(boundaryCrossed) > 0
        % % Set xyz_Global_A <-- Ia
        % % Set xyz_Global_B <-- Ib
        xyz_Global_A = Ia;
        xyz_Global_B = Ib;
        
        xG = xyz_Global(1); yG = xyz_Global(2); zG = xyz_Global(3);
        xG_A = xyz_Global_A(1); yG_A = xyz_Global_A(2); zG_A = xyz_Global_A(3);
        xG_B = xyz_Global_B(1); yG_B = xyz_Global_B(2); zG_B = xyz_Global_B(3);
        
        % % Use xyz_Global_A and xyz_Global_B to reflect new product position, 
        % % xyz_Global, about the new axis defined by the line
        % % IaIb. This means we need to evaluate the new Ia and Ib 
        % % boundary intersection points.
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
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
        % -------------------------------------------------------------------------
        % Determine the distance from B to each potential boundary, and also
        % determine the point at which BP intersects each potential boundary.
        d_Ib_list = zeros(1,length(beyond_cL));
        Ib_list = zeros(3,length(beyond_cL));
        for j = 1:length(beyond_cL)
            % Determine through which boundary the line AP intersects first.
            if beyond_cL(j) == 3 %zG > cL
                % Top
                del = cL - zG_B; % delZ is positive
                DEL = zG - zG_B; % DELZ is positive

            elseif beyond_cL(j) == -3 %-zG > cL
                % Bottom   
                del = -cL - zG_B; % delZ is negative
                DEL = zG - zG_B; % DELZ is negative

            elseif beyond_cL(j) == -2 %-yG > cL
                % Left
                del = -cL - yG_B; % delZ is negative
                DEL = yG - yG_B; % DELZ is negative

            elseif beyond_cL(j) == 2 %yG > cL
                % Right
                del = cL - yG_B; % delZ is positive
                DEL = yG - yG_B; % DELZ is positive

            elseif beyond_cL(j) == -1 %-xG > cL
                % Back
                del = -cL - xG_B; % delZ is negative
                DEL = xG - xG_B; % DELZ is negative

            elseif beyond_cL(j) == 1 %xG > cL
                % Front
                del = cL - xG_B; % delZ is positive
                DEL = xG - xG_B; % DELZ is positive  
            end
            d_Ib_list(j) = dist_BtoProduct* (del/DEL); % distance from A to Ia 
            normalVec = ([xG;yG;zG] - [xG_B;yG_B;zG_B]) /... 
            norm([xG;yG;zG] - [xG_B;yG_B;zG_B]);

            Ib_list(:,j) = [xG_B;yG_B;zG_B] + d_Ib_list(j)*normalVec;
        end
        [distToFirstBoundaryCrossed_A,indexA] = min(d_Ia_list);
        Ia = Ia_list(:,indexA);

        [distToFirstBoundaryCrossed_B,indexB] = min(d_Ib_list);
        Ib = Ib_list(:,indexB);


        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        
    else
        product_still_outside_boundary = 0; % We are done reflecting.        
        
    end
    % ------------ END While Loop -----------------------------------------
    end
    % ------------ END While Loop -----------------------------------------
    
    % Assert that xyz_Global is now in the main simulation box. 
    assert( sum( abs(xyz_Global) > cL*ones(size(xyz_Global)) ) ==0) 

end




%}
end
