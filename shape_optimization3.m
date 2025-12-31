function [V ] = shape_optimization3(max_iter , V,  V_original, F, angle_threshold, saliency, fixed_vertex_indices )
    
% V_original - the best oriented vertex coordinates from the part orientation change
% V - variable in the shape optimisation
% fixed_vertex_indices - fixations ( no movement applied)

% Author - Joseph Jayakody

    %% settings
    V2 = V; V_prev = V ;

    E_prev = 0;
    %tolerance = 1e-5;
    energy_tolerance = 1e-6;     
    max_allowed_deviation = 0.0075;
    
    support_weights = zeros(size(V,1), 1);
    v1 = V(F(:,1), :);
    v2 = V(F(:,2), :);
    v3 = V(F(:,3), :);
    face_areas = 0.5 * vecnorm(cross(v2 - v1, v3 - v1, 2), 2, 2);
    face_centers = (v1 + v2 + v3)/3;        % face centroids
    face_heights = face_centers(:,3);  
     
    [face_normals] = normals(V, F);

    %% loop
    for iter = 1:max_iter   
        nV = size(V, 1); %changed
        vertexProjections = cell(nV, 1);
        
        
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
            nz = abs(n(3)); % dot(n, [0 0 1]
            % projected face area (equivalent to dot product between [001] and face normal) 
            a_proj = a * nz;

             support_weights(verts) = support_weights(verts) + (a_proj*h);  % accumulate
        
        end
        
        support_weights = round(support_weights);
    
        facedown_faces = F(facedown_indices, :); % m x 3
        % Flatten to get all vertex indices involved in facedown faces
        facedown_vertex_indices = facedown_faces(:);   % (3m x 1) vector
        
        % Number of appearances for each vertex in facedown faces
        support_counts = accumarray(facedown_vertex_indices, 1, [nV, 1]);
        % Avoid divide by zero for overhang-free vertices  
        support_counts(support_counts == 0) = 1;
        support_weights = support_weights ./ support_counts;
    
        % Scaling down the support weights to avoid drastic changes/mesh collasping
        support_weights = support_weights / numel(facedown_indices);
        
        if isempty(facedown_indices)
            fprintf('No overhanging faces remain at iteration %d.\n', iter);
            break;
        end
    
        %% ELEMENT PROJECTION STEP
        for idx = 1:length(facedown_indices)
            fIdx = facedown_indices(idx);
            face_v = F(fIdx, :);
            projected_face_vertices = project_face_onto_cone_3(V(face_v, :), angle_threshold);
            
            for j = 1:length(face_v)
                vertex_idx = face_v(j);
                vertexProjections{vertex_idx} = [vertexProjections{vertex_idx}; projected_face_vertices(j, :)];
            end
        end
    
        V_proj = zeros(size(V,1),3);
        
        %   NEW PROJECTION - saliency aware
    
        for i = 1:nV
            if ~isempty(vertexProjections{i})
                Vp_mean = mean(vertexProjections{i}, 1); % mean of projections for this vertex
                s = saliency(i);      % saliency at this vertex (in [0,1])
        
                % Blend between projected target and original vertex
                V_proj(i, :) = (1 - s) * Vp_mean + s * V(i, :);
            else
                % If no projections (i.e., not incident to any overhanging face), keep it unchanged
                V_proj(i, :) = V(i, :);
            end
        end
        
    
        %% Weight setting 
        %W_proj should decrease (lowering penalty for deviation from projection)
        % W_close, W_smooth should increase 
    
        k_proj = 0.05;
        k_smooth = 0.000002;
        %k_close = 0.0001;
        
       % Projection term: decreases over time
        w_proj   = 0.0005* exp(-k_proj * iter);           %alpha
        %w_proj = (1 - exp(-k_proj * iter)); 
    
        % Smoothness term: increases over time
        w_smooth =  (1 - exp(-k_smooth * iter));   %beta
    
        % Closeness term: increases over time
        %w_close  =  (1 - exp(-k_close * iter));    %lambda
        w_sup  = exp(-0.1 * iter); 
       
       %L = cotmatrix_embedded(V, F);    % Laplacian edge length-based
        L = cotmatrix(V, F);            % laplacian vertex position-based        
        
         
        A = w_proj* speye(nV) + w_smooth * (L' * L);
        %% updated with displacement smoothing (not surface smoothing)
        rhs_x = w_proj * V_proj(:,1)   + w_smooth * 0.001*(L' * L) * V2(:,1);
        rhs_y = w_proj * V_proj(:,2)   + w_smooth * 0.001*(L' * L) * V2(:,2);
        rhs_z = w_proj * V_proj(:,3)   + w_smooth * 0.001*(L' * L) * V2(:,3)- 0.0001*  support_weights;
        
        %% x(N) projection
   
%         A = 5*w_proj* speye(nV) + w_smooth * (L' * L)   ;
%       
%         %% updated with displacement smoothing (not surface smoothing)
%         rhs_x = 5*w_proj * V_proj(:,1)   + w_smooth * 0.001*(L' * L) * V2(:,1);
%         rhs_y = 5*w_proj * V_proj(:,2)   + w_smooth * 0.001*(L' * L) * V2(:,2);
%         rhs_z = 5*w_proj * V_proj(:,3)   + w_smooth * 0.001*(L' * L) * V2(:,3)- 0.0001*  support_weights;

    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fixing vertices
    
        % 1. Removes influence of energy terms for fixed vertices
        % setting the entire row to 0.
        A(fixed_vertex_indices, :) = 0;  
    
        % 2. Here diagonal elements of A (of fixed vertex rows) are set to 1
        A(sub2ind(size(A), fixed_vertex_indices, fixed_vertex_indices)) = 1;
        
        % The below lines force  that vertices are fixed to the original position
        rhs_x(fixed_vertex_indices) = V(fixed_vertex_indices, 1);
        rhs_y(fixed_vertex_indices) = V(fixed_vertex_indices, 2);
        rhs_z(fixed_vertex_indices) = V(fixed_vertex_indices, 3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
        % Solve for V directly (not displacements)
        V(:,1) = A \ rhs_x;
        V(:,2) = A \ rhs_y;
        V(:,3) = A \ rhs_z;    
        
        current_support = compute_weighted_support_proxy5(V, F, angle_threshold);
        fprintf('Iteration %d: %d support volume.\n', iter,  current_support);
      
        % === Check convergence ===
        
        delta = norm(V - V_prev, 'fro') / sqrt(numel(V));  % normalized change
    
        % === Compute energy (optional for energy-based stopping) ===
        E_proj = w_proj * sum(vecnorm(V - V_proj, 2, 2).^2);
        E_smooth = w_smooth * sum(vecnorm(L*(V - V2), 2, 2).^2);
        E_supp = w_sup * sum(support_weights .* V(:,3));  % only Z influences support
    
        E_curr = E_proj + E_smooth +  E_supp ;
        energy_change = abs(E_curr - E_prev);
    
        fprintf('Iter %d | ΔV = %.6f | ΔE = %.6f \n', ...
            iter, delta, energy_change );
        
        % === Check deviation from initial geometry ===

        if any(isnan(V), 'all') || any(isinf(V), 'all')
            error('Vertex coordinates became NaN or Inf during optimization.');
        end

        assert(size(V,1) == size(V_original,1), 'Vertex count mismatch between current and original mesh.');
        
        deviation = norm(V - V_original, 'fro')^2 / nV;  % mean squared deviation
    
        fprintf('Iter %d | Deviation from V0 = %.6f\n', iter, deviation);
    
        if energy_change < energy_tolerance
            disp('Converged.');
            break;
        end
    
        if deviation > max_allowed_deviation
            warning('Mesh deviated too far from initial geometry. Terminating.');       
            break;
        end
    
        E_prev = E_curr;
        V_prev = V;
    
    end
end
     
 