function [V_opt, best_support, saved_shapes] = optimise_shape_and_orientation3(V, F, fixed_v, max_iter, angle_threshold, saliency,current_sup)
% OPTIMIZE_SHAPE_AND_ORIENTATION
% Alternates between global orientation optimization and local shape optimization 

%% Variables to select - number of unit sphere samples, max deviation

% INPUTS:
%   V               - (Nx3) Initial vertex positions.
%   F               - (Mx3) Triangular face indices.
%  fixed_ v         - fixed vertex indices
%   max_iter        - Maximum number of alternating iterations.
%   angle_threshold - Overhang angle threshold (e.g., 120°).
 
% OUTPUTS:
%   V_opt - Optimized vertex positions.
%   best support - best support volume achieved.

     
 
    % no of part orientation optimising iterations
    temp_iter = 2;
    
    % initial support volume setting - Inf (for any model)
    %temp_best_vol = Inf;
    
    % current support volume is useful for certain models like Kitten
    temp_best_vol =current_sup;
    
    % Add this at the start of your function
    saved_shapes = {};  % Cell array to store updated V
   
    
    % Optimization Loop
    for iter = 1:max_iter


        fprintf('Initial Support Volume: %.3f\n', temp_best_vol);

        for j = 1: temp_iter
        %% --- STEP 1: GLOBAL ORIENTATION OPTIMIZATION ---
         
            % here the numer of samples of unit sphere sampling has to be
            % set (either 500 or 1000)
            [opt_V] = hybrid_xy_orientation_optimization3(V, F, angle_threshold);
             
            [rotate_support_V] = compute_weighted_support_proxy5(opt_V, F, angle_threshold);

            if rotate_support_V < temp_best_vol
                
                fprintf('Iteration %d: Orientation updated. New Support Volume = %.3f\n', iter, rotate_support_V);
                V = opt_V;     
                temp_best_vol = rotate_support_V;
                saved_shapes{end+1} = V;
                
            else
                fprintf('Iteration %d: Orientation is not updated\n', iter);
            end
        end
        
        V_original = V;

        %% --- STEP 2: LOCAL SHAPE OPTIMIZATION ---
        % Change the number of iterations based on the model 
        % Kitten ~ 20, Armadillo ~ 8-10
        shape_iter = 20;   
        [V_new ] = shape_optimization3(shape_iter, V, V_original, F, angle_threshold,  saliency, fixed_v);
        saved_shapes{end+1} = V_new;
        V = V_new;
    end

    V_opt = V; % Final optimized shape
    [best_support] = compute_weighted_support_proxy5(V_opt, F, angle_threshold);
end


