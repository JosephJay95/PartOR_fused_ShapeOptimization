function support_proxy = compute_weighted_support_proxy5(V, F, angle_threshold)
 
 
    % Compute face area  
    v1 = V(F(:,1), :);
    v2 = V(F(:,2), :);
    v3 = V(F(:,3), :);

    % note this is 1/2 *a *b  (traingle area)
    face_areas = 0.5 * vecnorm(cross(v2 - v1, v3 - v1, 2), 2, 2);
    
    % so we can compute triangular prism volume = 1/2*a*b* H , H is prism height
    [face_normals] = normals(V, F);
    build_d = [0 0 1];
     

    % Ensure build_d is repeated across all faces
    build_d_repeated = repmat(build_d, size(face_normals, 1), 1);
    
    face_angles = atan2d( vecnorm(cross(face_normals, build_d_repeated, 2), 2, 2), ...
       dot(face_normals, build_d_repeated, 2));
    

    facedown_indices = find(face_angles > angle_threshold);
    
    %% Compute height of face centroids (Z-coordinate)
    face_centers = (v1 + v2 + v3) / 3;
    face_heights = face_centers(:, 3);  % Z-coordinates
    h = face_heights(facedown_indices);
    a = face_areas(facedown_indices);

    w1 = 1000 ; w2 = 50;

    support_proxy =  w1 * numel(facedown_indices) + sum(w2*a.*h)/numel(facedown_indices)  ;

end
 
