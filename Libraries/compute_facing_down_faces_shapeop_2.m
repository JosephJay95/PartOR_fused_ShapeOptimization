function [vertex_scalar_field, overhang_FA, face_angles, facedown_faces,overhang_indices] = compute_facing_down_faces_shapeop_2(V, F)
    

    %% Notes
    % Overhang_FA - Angles of the facing down faces
    % face_angles - angles of all faces
    % vertex scalar field - Interpoalted overhang values usign Barycentric

    % Step 1: Compute face normals
    [N] = normals(V, F);

    % Step 2: Z-axis (print direction)
    build_d = [0, 0, 1];
    
    % Step 3: Compute angles between face normals and Z-axis
    face_angles = zeros(size(F, 1), 1); % Preallocate for angles
    for i = 1:size(F, 1)
        % Compute angle using dot product and cross product for robustness
        face_angles(i) = atan2d(norm(cross(N(i, :), build_d)), dot(N(i, :), build_d));
    end
    
    tol = min(V(:,3))+ 1;
    

    % Step 4: Identify "facing-down" faces.
    % Here we assume faces with an angle >= 120 degrees are "facedown".
    
    facedown_face_indices = find(face_angles >= 120);
    % Remove faces that are in the bounding box
    
    % Step 5: Identify "in-BBox" faces
    % in_bbox_id = exclude_faces_in_bbox(V, F, tol);
    % facedown_face_indices = setdiff(facedown_face_indices, in_bbox_id);
    

    facedown_faces = F(facedown_face_indices, :);
    
    % Step 6: Initialize vertex scalar field (one value per vertex)
    nV = size(V,1);
    vertex_scalar_field = zeros(nV, 1);
    % Preallocate overhang_FA array for facedown faces
    overhang_FA = zeros(length(facedown_face_indices), 1);

    % Step 7: Compute barycentric weights and assign scalar values
    for idx = 1:length(facedown_face_indices)
        face_idx = facedown_face_indices(idx);
        face = F(face_idx, :); % Current face vertices (indices)
        theta = face_angles(face_idx); % Angle for the current face
        
        overhang_FA(idx) = theta; 
        % Vertices of the face
        v1 = V(face(1), :);
        v2 = V(face(2), :);
        v3 = V(face(3), :);

        % Compute barycentric coordinates for each vertex
        % Use the centroid as a reference for uniform weights
        centroid = (v1 + v2 + v3) / 3;
        area = norm(cross(v2 - v1, v3 - v1)) / 2; % Face area

        % Assign weights proportionally
        w1 = norm(cross(v2 - centroid, v3 - centroid)) / (2 * area);
        w2 = norm(cross(v3 - centroid, v1 - centroid)) / (2 * area);
        w3 = norm(cross(v1 - centroid, v2 - centroid)) / (2 * area);

        % Accumulate scalar values
        vertex_scalar_field(face(1)) = vertex_scalar_field(face(1)) + theta * w1;
        vertex_scalar_field(face(2)) = vertex_scalar_field(face(2)) + theta * w2;
        vertex_scalar_field(face(3)) = vertex_scalar_field(face(3)) + theta * w3;
    end
     % Also, return the unique vertex indices that are part of facedown faces.
    overhang_indices = unique(facedown_faces(:));

end
