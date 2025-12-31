
clear; close all; clc;
addpath(genpath('Libraries'))
 

model = createpde(1);
importGeometry(model, 'kitten_big.stl'); 
b = generateMesh(model, 'Hmax', 2, 'GeometricOrder', 'linear');
V = b.Nodes';    % Vertex positions (Nx3)
T = b.Elements'; % Triangular elements (indices)


%Compute boundary faces and extract unique boundary vertices
F = boundary_faces(T);
boundary_p_indices = unique(F);
boundary_p = V(boundary_p_indices, :);
V = boundary_p;  


V2=V;

%% Main Optimization Loop (Least-Squares Projection)
max_iter = 50;  % Set number of iterations
angle_threshold = 120;  % Constraint: face normal must be within 120° of vertical
h_ref = min(V(:,3)); 

% For simplicity, we update the mesh by averaging projections from all faces sharing a vertex.
% Create a cell array to store projections for each vertex.
nV = size(V,1);
vertexProjections = cell(nV,1);

% Variable V will be used in the optimisation loop, 
V1 = V; V_best = V;

[Mesh_saliency_multiscale, ~, ~] = compute_saliency (V2, F); 
N_size = 1; %neighbourhood size


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fixed region selection

% method 1 - Random point selection
%selected_p = 4;

% method 2 - height range point selection (REQUIRES THE GEODESIC FIELD LIBRARY
% FOR THIS)

z_min = 0; z_max = 1;
p_filtered = V(V(:,3) >= z_min & V(:,3) <= z_max, :);
% Find intersection and then corresponding indices
[commonPoints, ia, ib] = intersect(V, p_filtered, 'rows');
[~, selected_p] = ismember(commonPoints, V, 'rows');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[points] = compute_geodesic_neighbourhood(V, F, selected_p, N_size);
fixed_vertex_indices = cell2mat(points);
fixed_region = V(fixed_vertex_indices,:);
figure (10);
patch('Vertices', V2, 'Faces', F, 'FaceColor', 'green', 'EdgeColor', 'none'); hold on;
scatter3 (fixed_region(:,1),fixed_region(:,2),fixed_region(:,3), 'r' );
axis equal; view(3); camlight; lighting gouraud;
title('Original mesh');

[~, ~, ~, ~,overhangIndices] = compute_facing_down_faces_shapeop_2(V, F);

face_angles = zeros(size(F, 1), 1); % Preallocate for angles

% running with the new proxy
global_best_support = compute_weighted_support_proxy3(V, F, angle_threshold);


% support weight matrix
support_weights = zeros(nV, 1);
support_weights_wrong = zeros(nV, 1);

v1 = V(F(:,1), :);
v2 = V(F(:,2), :);
v3 = V(F(:,3), :);
face_areas = 0.5 * vecnorm(cross(v2 - v1, v3 - v1, 2), 2, 2);
face_centers = (v1 + v2 + v3)/3;        % face centroids
face_heights = face_centers(:,3);  
% Z height
face_weights = face_areas.*face_heights;

%face normals
[face_normals] = normals(V, F);

% 
V_prev = V ; E_prev = 0;
tolerance = 1e-5;
energy_tolerance = 1e-6;
 
max_allowed_deviation = 0.1;


tic;

%% OPTIMISATION LOOP
for iter = 1:max_iter
   
    % Clear previous projections for this iteration.
    for i = 1:nV
        vertexProjections{i} = [];
    end

 
    %[face_normals] = normals(V, F);
    build_d = [0 0 1];
    
    build_d_repeated = repmat(build_d, size(face_normals, 1), 1);
    
    face_angles = atan2d( vecnorm(cross(face_normals, build_d_repeated, 2), 2, 2), ...
       dot(face_normals, build_d_repeated, 2));
    
    facedown_indices = find(face_angles > angle_threshold);
    
    %% Support weight matrix computation
    for idx = 1:length(facedown_indices)
        f_idx = facedown_indices(idx);        % face index
        verts = F(f_idx, :);                  % 3 vertex indices
        a = face_areas(f_idx);
        h = face_heights(f_idx);              % centroid Z-height
        
        n = face_normals(f_idx, :);  % face normal
        nz = abs(n(3)); % dot(n, [0 0 1])

        % projected face area (equivalent to dot product between [001] and face normal) 
        a_proj = a * nz;

        support_weights(verts) = support_weights(verts) + (a_proj*h);  % accumulate

        %support_weights_wrong = support_weights_wrong(verts) + (a *h); 
    end
    
    support_weights = round(support_weights);
    
    % normalisation - remove if not needed

    facedown_faces = F(facedown_indices, :); % m x 3
    % Flatten to get all vertex indices involved in facedown faces
    facedown_vertex_indices = facedown_faces(:);   % (3m x 1) vector
    
    % Count how many times each vertex appears in facedown faces
    support_counts = accumarray(facedown_vertex_indices, 1, [nV, 1]);
    % Avoid divide by zero for vertices not involved in overhang
    support_counts(support_counts == 0) = 1;
    support_weights = support_weights ./ support_counts;

    % why this step? this is important to scale down the support weights to
    % avoid drastic changes/mesh collasping
    support_weights = support_weights / numel(facedown_indices);
    
    %OV_F = F (facedown_indices,:);
    %exceeding_angles = overhang_angles (facedown_indices,:);

    % Display progress and stop if no overhang exists
    if isempty(facedown_indices)
        fprintf('No overhanging faces remain at iteration %d.\n', iter);
        break;
    end
    
    %% PROJECTION STEP
    % For each overhanging face, compute the least-squares projection.
    for idx = 1:length(facedown_indices)
        fIdx = facedown_indices(idx);
        face_v = F(fIdx, :);  % indices of vertices for face fIdx
        current_face_vertices = V(face_v, :);
        
        % Compute the projection onto the constraint set (cone)
        projected_face_vertices = project_face_onto_cone_3(current_face_vertices, angle_threshold);
        
        % For each vertex in the face, store the projection.
        for j = 1:length(face_v)
            vertex_idx = face_v(j);
            vertexProjections{vertex_idx} = [vertexProjections{vertex_idx}; projected_face_vertices(j, :)];
        end
    end

    %% Now update each vertex as the average of its projected positions -This needs more clarification
    V_proj = zeros(size(V,1),3);  % initialize projected vertices
    %V_proj =V;

%     for i = 1:nV
%         if ~isempty(vertexProjections{i})
%             V_proj(i, :) = mean(vertexProjections{i}, 1);  
%         else
%             % If the vertex was not involved in any overhanging face, keep it unchanged
%             V_proj(i, :) = V(i, :);
%         end
%     end

    for i = 1:nV
        if ~isempty(vertexProjections{i})
            Vp_mean = mean(vertexProjections{i}, 1); % mean of projections for this vertex
            s = Mesh_saliency_multiscale(i);                         % saliency at this vertex (in [0,1])
    
            % Blend between projected target and original vertex
            V_proj(i, :) = (1 - s) * Vp_mean + s * V(i, :);
        else
            % If no projections (i.e., not incident to any overhanging face), keep it unchanged
            V_proj(i, :) = V(i, :);
        end
    end


    % Global update - currently we are not using regularisation weight
    % So, we simply update V = V_proj.
    V_temp = V_proj;
    
    % fixed weights
    %     beta = 0.4; 
    %     alpha = 0.3; 
    %     lambda_close = 0.1; % Closeness weight (NEW ADDITION)

    %adaptive weights
    %gamma_proj = 1;
    k_proj = 0.05;
    k_smooth = 0.000002;
    k_close = 0.0001;
    
    % Projection term: decreases over time
    w_proj   = 0.0005* exp(-k_proj * iter);           %alpha
    %w_proj = (1 - exp(-k_proj * iter)); 

    % Smoothness term: increases over time
    w_smooth =  (1 - exp(-k_smooth * iter));   %beta

    % Closeness term: increases over time
    w_close  =  (1 - exp(-k_close * iter));    %lambda

    w_sup  = exp(-0.1 * iter);  
    
    %S = spdiags(( Mesh_saliency_multiscale), 0, nV, nV);  % Diagonal saliency weights
    %L = cotmatrix_embedded(V, F);                      % Laplacian
    L = cotmatrix(V, F);

    %support_weights =   support_weights/ size(facedown_indices);
    Supp = spdiags(support_weights, 0, nV, nV);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

    
    % damping factor for models
    % 0.0001 - Kitten

    A = w_proj* speye(nV) + w_smooth * (L' * L)   ;
  
    %% updated with displacement smoothing (not surface smoothing)
    rhs_x = w_proj * V_proj(:,1)   + w_smooth * 0.001*(L' * L) * V2(:,1);
    rhs_y = w_proj * V_proj(:,2)   + w_smooth * 0.001*(L' * L) * V2(:,2);
    rhs_z = w_proj * V_proj(:,3)   + w_smooth * 0.001*(L' * L) * V2(:,3)- 0.0001*  support_weights;
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fixing vertices
    A(fixed_vertex_indices, :) = 0;
    A(sub2ind(size(A), fixed_vertex_indices, fixed_vertex_indices)) = 1;
    
    rhs_x(fixed_vertex_indices) = V(fixed_vertex_indices, 1);
    rhs_y(fixed_vertex_indices) = V(fixed_vertex_indices, 2);
    rhs_z(fixed_vertex_indices) = V(fixed_vertex_indices, 3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Inside your optimization loop, just before solving for V(:,1:3)
    % Enforce sparse LU factorization explicitly
    % Inside optimization loop, for each iteration
%     [L,U,P] = lu(A,'vector');
%     
%     V(:,1) = U \ (L \ rhs_x(P));
%     V(:,2) = U \ (L \ rhs_y(P));
%     V(:,3) = U \ (L \ rhs_z(P));


    % Solve for V directly (not displacements)
    V(:,1) = A \ rhs_x;
    V(:,2) = A \ rhs_y;
    V(:,3) = A \ rhs_z;
    
 


   % fprintf('Iteration %d: %d overhanging faces processed.\n', iter, length(facedown_indices));

   % --- Check matrix properties ---
    if issparse(A)
        fprintf('Matrix is sparse.\n');
    else
        fprintf('Matrix is dense.\n');
    end
    
    [nRows,nCols] = size(A);
    if nRows == nCols
        fprintf('Matrix is square.\n');
    else
        fprintf('Matrix is rectangular → MATLAB will use QR factorization.\n');
    end
    
    if issymmetric(A)
        fprintf('Matrix is symmetric.\n');
        try
            R = chol(A);  % Try Cholesky
            fprintf('Cholesky factorization succeeds → MATLAB would use Cholesky.\n');
        catch
            fprintf('Cholesky fails → MATLAB would fall back to LU factorization.\n');
        end
    else
        fprintf('Matrix is not symmetric → MATLAB will use LU (if square) or QR (if rectangular).\n');
    end
 
     % new proxy
    current_support = compute_weighted_support_proxy5(V, F, angle_threshold);
    fprintf('Iteration %d: %d support volume.\n', iter,  current_support);
  
    % === Check convergence ===
    delta = norm(V - V_prev, 'fro') / sqrt(numel(V));  % normalized change

    % === Compute energy (optional for energy-based stopping) ===
    E_proj = w_proj * sum(vecnorm(V - V_proj, 2, 2).^2);
    E_smooth = w_smooth * sum(vecnorm(L*(V - V2), 2, 2).^2);
    E_supp = w_sup * sum(support_weights .* V(:,3));  % only Z influences support

    E_curr = E_proj + E_smooth +E_supp ;
    energy_change = abs(E_curr - E_prev);

    fprintf('Iter %d | ΔV = %.6f | ΔE = %.6f \n', ...
        iter, delta, energy_change );
    
    % === Check deviation from initial geometry ===
    deviation = norm(V - V2, 'fro')^2 / nV;  % mean squared deviation

    fprintf('Iter %d | Deviation from V0 = %.6f\n', iter, deviation);

    if delta < tolerance || energy_change < energy_tolerance
        disp('Converged.');
        break;
    end

     if deviation > max_allowed_deviation
        warning('Mesh deviated too far from initial geometry. Terminating.');
        break;
    end

    E_prev = E_curr;
    V_prev = V;

    
    E_curr_saved (iter,:)=E_curr;

end

toc;


figure (2);
patch('Vertices', V  , 'Faces', F, 'FaceColor', 'green', 'EdgeColor', 'none');
axis equal; view(3); camlight; lighting gouraud;
title('Optimized Mesh via Least-Squares Projection');

% Final plotting of the optimized mesh.
figure (1);
patch('Vertices', V1, 'Faces', F, 'FaceColor', 'green', 'EdgeColor', 'none');
axis equal; view(3); camlight; lighting gouraud;
title('Original mesh');



