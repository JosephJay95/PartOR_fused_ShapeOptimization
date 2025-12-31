function [facedown_faces,facedown_indices] = compute_facing_down_and_overhang(V, F, tol)

    % Compute face normals
    %N = compute_F_normals(V, F); % Function to compute face normals, assume provided

    %GP toolbox function for per face normal vectors
    [ N ] = normals(V, F);

    % Z-axis (print direction)
    build_d = [0, 0, 1];
    
    % Compute angles between face normals and Z-axis
    face_angles = zeros(size(F, 1), 1); % Preallocate for angles
    for i = 1:size(F, 1)
        % Compute angle using dot product and cross product for robustness
        face_angles(i) = atan2d(norm(cross(N(i, :), build_d)), dot(N(i, :), build_d));
    end
    
    % Identify facing down faces
    % A face is "facing down" if the normal points opposite to the z-axis
    % i.e., the angle is greater than 90 degrees
    
    % Define overhang threshold (customizable, e.g., 45 degrees)
    overhang_threshold = 45;
    
    % Step 4: Identify "in-BBox" faces
    in_bbox_id= exclude_faces_in_bbox(V, F, tol);


    % Identify overhang faces
    % An overhang face has an angle larger than the threshold but less than 90
    facedown_indices = find(face_angles > overhang_threshold & face_angles >= 120);
    
    % Step 5: Remove "in-BBox" faces from facing-down and overhang lists
    facedown_indices = setdiff(facedown_indices, in_bbox_id);

    facedown_faces = F(facedown_indices,:);

    % no of facedownfaces
    num = numel(facedown_indices);
    
end



