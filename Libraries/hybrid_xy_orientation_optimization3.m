function [opt_V] = hybrid_xy_orientation_optimization3(V, F, angle_threshold)
        % INPUTS:
    % V - Nx3 matrix of vertices
    % F - Mx3 matrix of face indices
    % N - number of directions to sample on the sphere
    % OUTPUT:
    % rotated_meshes - cell array of rotated vertex matrices
    
    N = 1000;
     
    best_val = Inf;

    % Canonical "up" direction
    z_axis = [0 0 1]';

    % Generate uniformly spaced directions using Fibonacci sampling
    dirs = fibonacci_sphere(N);  % Nx3 matrix of directions

    for i = 1:N
        dir = dirs(i, :)';  % Convert to column vector

        % Compute rotation axis and angle
        axis = cross(z_axis, dir);
        axis_norm = norm(axis);

        if axis_norm < 1e-8
            % Vectors are (almost) parallel
            if dot(z_axis, dir) > 0
                R = eye(3);  % No rotation
            else
                % 180° rotation around any axis perpendicular to z (e.g., x-axis)
                R = [-1  0  0;
                      0 -1  0;
                      0  0  1];
            end
        else
            axis = axis / axis_norm;
            angle = acos(dot(z_axis, dir));
            R = axis_angle_to_matrix(axis, angle);
        end

        % Rotate vertices
        V_rot = (R * V')';

        val = compute_weighted_support_proxy5(V_rot, F, angle_threshold);
        %val =compute_weighted_support_proxy6_volumeapproxONLY(V, F, angle_threshold);
        
        if val < best_val
            best_val = val;
            opt_V = V_rot;
        end
        
    end
end


function dirs = fibonacci_sphere(N)
    phi = pi * (3 - sqrt(5));  % Golden angle in radians
    dirs = zeros(N, 3);

    for i = 1:N
        y = 1 - (i - 0.5) * (2 / N);  % y ranges from 1 to -1
        r = sqrt(1 - y^2);            % Radius at height y
        theta = phi * i;

        x = cos(theta) * r;
        z = sin(theta) * r;

        dirs(i, :) = [x, y, z];  % Store unit direction
    end
end

function R = axis_angle_to_matrix(axis, angle)
    % Rodrigues' rotation formula
    K = [     0     -axis(3)  axis(2);
           axis(3)     0     -axis(1);
          -axis(2)  axis(1)     0    ];
    R = eye(3) + sin(angle)*K + (1 - cos(angle))*(K*K);
end
