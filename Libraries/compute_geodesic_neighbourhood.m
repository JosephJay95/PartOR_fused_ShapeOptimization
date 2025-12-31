function [points] = compute_geodesic_neighbourhood(V, F, handle_indices, N_size)

    d_sigma = N_size;

    % we find the influence regions for all handles of the segment
    for i = 1:size(handle_indices, 1)
            
       % Compute geodesic distances from the current imp_point to all vertices
       [D, ~] = heat_geodesics_singlesource(V, F, handle_indices(i));

       % Compute verticeds within f-Geodesic radius , can extend to finding faces 
       points {i,:} = unique(find(D < d_sigma));
       %Inf_regions {i,:} = V(Inf_regions_indices,:);

    end

end