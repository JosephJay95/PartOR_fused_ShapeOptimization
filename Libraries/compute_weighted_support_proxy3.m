function support_proxy = compute_weighted_support_proxy3(V, F, angle_threshold)
 
    %% updated with the correct overhang face computation 
 

    % Compute face normals and areas
    v1 = V(F(:,1), :);
    v2 = V(F(:,2), :);
    v3 = V(F(:,3), :);
    face_areas = 0.5 * vecnorm(cross(v2 - v1, v3 - v1, 2), 2, 2);

    [face_normals] = normals(V, F);
    build_d = [0 0 1];
    h_ref = 0;

     % Ensure build_d is repeated across all faces
    build_d_repeated = repmat(build_d, size(face_normals, 1), 1);
    
    face_angles = atan2d( vecnorm(cross(face_normals, build_d_repeated, 2), 2, 2), ...
       dot(face_normals, build_d_repeated, 2));
 
    facedown_indices = find(face_angles > angle_threshold);
    

    % Support distance at each vertex in each face
    z1 = V(F(:,1),3);
    z2 = V(F(:,2),3);
    z3 = V(F(:,3),3);
    support_dists_sq = ((z1 - h_ref).^2 + (z2 - h_ref).^2 + (z3 - h_ref).^2) / 3;

    % Compute weighted proxy for overhanging faces
    % What is considered - the number of faces and the approx support
    % volume (V) ~ sum (area*(h)^2)
    w1 = 5; w2 = 0.01;
    support_proxy = w1* numel(facedown_indices) + round(w2* sum(face_areas(facedown_indices) .* support_dists_sq(facedown_indices)));
end
 
