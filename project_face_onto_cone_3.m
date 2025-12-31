function projected_face_vertices = project_face_onto_cone_3(X_f, angle_threshold, weights)
% PROJECT_FACE_ONTO_CONE_LEAST_SQUARES Projects a triangular face onto
% a constraint set using a least-squares rigid transformation.
%
% Input:
%   X_f: 3x3 matrix, where each row is a vertex [x, y, z] of the face.
%   angle_threshold: Maximum allowed angle (in degrees) between the face normal
%                    and the vertical direction [0 0 1].
%   weights: A vector of weights for least-squares fitting (default: [1 1 1]).
%
% Output:
%   projected_face_vertices: 3x3 matrix of new vertex positions after projection.

    if nargin < 3
        weights = ones(3,1); % Default weights if not provided
    end
    
    % Compute the centroid of the original face
    centroid_X = mean(X_f, 1);
    
    % Compute face normal
    edge1 = X_f(2, :) - X_f(1, :);
    edge2 = X_f(3, :) - X_f(1, :);
    normal_X = cross(edge1, edge2);
    normal_X = normal_X / norm(normal_X); % Normalize
    
    % Define the vertical direction
    vertical = [0, 0, 1];
    
    % Compute the current normal's angle with the vertical
    current_angle = acosd(dot(normal_X, vertical));
    
    % If within threshold, no projection is needed
    if current_angle <= angle_threshold
        projected_face_vertices = X_f;
        return;
    end
    
    %% The below process is if the face normal vector isn't within the threshold
    
    % Compute the target normal on the constraint cone
    target_cos_theta = cosd(angle_threshold);

    %now we have to project the face normal of the overhang face to the
    %constraint space (i.e., allowed overhang range)
    projected_normal = normal_X - ((dot(normal_X, vertical) - target_cos_theta) * vertical);
    projected_normal = projected_normal / norm(projected_normal);
    
    % Compute target vertex positions Y_f using the new normal
    % Find the plane passing through the centroid with this normal
    d = -dot(projected_normal, centroid_X);
    
    % Project each vertex onto this plane to get Y_f
    Y_f = zeros(3,3);
    for i = 1:3
        point = X_f(i, :);
        distance_to_plane = dot(projected_normal, point) + d;
        Y_f(i, :) = point - distance_to_plane * projected_normal;
    end
    
    % Compute the optimal rigid transformation using Horn’s method
    [R, t] = compute_optimal_rigid_transform(X_f, Y_f, weights);
    
    % Apply the rigid transformation
    projected_face_vertices = (R * X_f')' + t;
end

function [R, t] = compute_optimal_rigid_transform(X, Y, weights)
% COMPUTE_OPTIMAL_RIGID_TRANSFORM Finds the optimal rotation and translation 
% between two sets of 3D points using Horn's method.
%
% Input:
%   X: Original 3x3 vertex positions.
%   Y: Target 3x3 vertex positions.
%   weights: Weights for each point in least squares.
%
% Output:
%   R: 3x3 rotation matrix.
%   t: 1x3 translation vector.

    % Compute weighted centroids
    W = diag(weights);
    centroid_X = sum(W * X) / sum(weights);
    centroid_Y = sum(W * Y) / sum(weights);

    % Center the points around their centroids
    X_centered = X - centroid_X;
    Y_centered = Y - centroid_Y;

    % Compute weighted covariance matrix
    H = X_centered' * W * Y_centered;

    % Compute SVD of H
    [U, ~, V] = svd(H);

    % Compute rotation matrix
    R = V * U';
    
    % Ensure a right-handed coordinate system (det(R) should be +1)
    if det(R) < 0
        V(:,3) = -V(:,3);
        R = V * U';
    end

    % Compute translation
    t = centroid_Y - (R * centroid_X')';
end
