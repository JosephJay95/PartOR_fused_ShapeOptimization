%% This code is to test the full pipeline  - part orientation and shape optimisation

% The maximum geometric deviation can be changed from shape_optimization3.m
% function (inside optimise_shape_and_orientation3.m)

clear; close all; clc;
addpath(genpath('Libraries'))

model = createpde(1);
importGeometry(model, 'wiresphere.stl'); 
b = generateMesh(model, 'Hmax', 2, 'GeometricOrder', 'linear');
V = b.Nodes';    % Vertex positions (Nx3)
T = b.Elements'; % Triangular elements (indices)


%Compute boundary faces and extract unique boundary vertices
F = boundary_faces(T);
boundary_p_indices = unique(F);
boundary_p = V(boundary_p_indices, :);
V = boundary_p;  
%[V,F] = readOBJ('new_venus.obj');
%[V,F] = readSTL('kitten_big.stl');

[saliency, ~, ~] = compute_saliency (V, F); 

max_iter = 3;  % Set number of iterations
angle_threshold = 120;  % Constraint: face normal must be within 120° of vertical
h_ref = min(V(:,3)); 

%% Base fixations - Kitten, Armadillo **********************************
z_min = 0; z_max = 2;
N_size = 0.1;


%z_min = -14.3; z_max = -13.5;

 selected_p = 4;

% z_min = min(V(:,3)); z_max = z_min+1; N_size = 0.1;
% 
% p_filtered = V(V(:,3) >= z_min & V(:,3) <= z_max, :);
% % Find intersection and then corresponding indices
% [commonPoints, ia, ib] = intersect(V, p_filtered, 'rows');
% [~, selected_p] = ismember(commonPoints, V, 'rows');

% p_filtered = V(V(:,2) >= z_min & V(:,2) <= z_max, :);
% % Find intersection and then corresponding indices
% [commonPoints, ia, ib] = intersect(V, p_filtered, 'rows');
% [~, selected_p] = ismember(commonPoints, V, 'rows');

% z_min = max(V(:,3))-2; z_max = max(V(:,3));
%  
% p_filtered = V(V(:,3) >= z_min & V(:,3) <= z_max, :);
% % Find intersection and then corresponding indices
% [commonPoints, ia, ib] = intersect(V, p_filtered, 'rows');
% [~, selected_p2] = ismember(commonPoints, V, 'rows');
% 
% selected_p = [selected_p1;selected_p2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[points] = compute_geodesic_neighbourhood(V, F, selected_p, N_size);
fixed_vertex_indices = cell2mat(points);

 
fixed_region = V(fixed_vertex_indices,:);
figure (10);
patch('Vertices', V, 'Faces', F, 'FaceColor', 'green', 'EdgeColor', 'none'); hold on;
scatter3 (fixed_region(:,1),fixed_region(:,2),fixed_region(:,3), 'r' );
axis equal; view(3); camlight; lighting gouraud;
title('Original mesh');

%% QUICK CHECK OF NO OF OVERHANG FACES
%[V,F] = readSTL('lambda_overhangfaces_weight50.stl');
[face_normals] = normals(V, F);
build_d = [0 0 1];
h_ref = 0;

 % Ensure build_d is repeated across all faces
 build_d_repeated = repmat(build_d, size(face_normals, 1), 1);
    
face_angles = atan2d( vecnorm(cross(face_normals, build_d_repeated, 2), 2, 2), ...
  dot(face_normals, build_d_repeated, 2));
    

facedown_indices = find(face_angles > angle_threshold);
numel(facedown_indices)    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART ORIENTATION PHASE CHECK

%just testing the hybrid orientation - multiscale+PS
% We have to set the overhang face percentage inside HYBRID_XY func
tic;
[opt_V] = hybrid_xy_orientation_optimization3(V, F, angle_threshold);
 
toc;
figure (3);
patch('Vertices', opt_V , 'Faces', F, 'FaceColor', 'green', 'EdgeColor', 'none');
axis equal; view(3); camlight; lighting gouraud;


%% new test (updated version)
tic;
current_sup = Inf;
[V_opt, best_support,saved_shapes] = optimise_shape_and_orientation3(V, F, fixed_vertex_indices, max_iter, angle_threshold, saliency,current_sup);
toc;

figure (2);
patch('Vertices', V_opt , 'Faces', F, 'FaceColor', 'green', 'EdgeColor', 'none');hold on;
fixed_region = V_opt(fixed_vertex_indices,:);
hold on;scatter3 (fixed_region(:,1),fixed_region(:,2),fixed_region(:,3), 'r' );

axis equal; view(3); camlight; lighting gouraud;
title('Optimized Mesh via Least-Squares Projection');
