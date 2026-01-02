function result = applyFibrosisThreshold(noise_field, geometry, target_density, bz_cfg)
% Applies thresholding logic using the PRE-COMPUTED geometry.

    % 1. Determine Core Threshold (Bisection Method)
    % Only look at pixels inside the core
    core_values = noise_field(geometry.core_indices);

    % Bisection setup
    low = min(core_values);
    high = max(core_values);
    max_iters = 40;
    tol = 1e-3;

    for i = 1:max_iters
        thresh = (low + high) / 2;
        current_dens = sum(core_values >= thresh) / geometry.num_core_points;

        if abs(current_dens - target_density) < tol
            break;
        end

        if current_dens < target_density
            high = thresh; % Need lower threshold to get more density
        else
            low = thresh;
        end
    end
    core_threshold = thresh;
    actual_core_density = current_dens;

    % 2. Build Full Threshold Map
    % Initialize with a value > 1 (impossible to cross for normalized noise)
    threshold_map = ones(size(noise_field)) * 999;

    % Set Core
    threshold_map(geometry.core_mask) = core_threshold;

    % Set Border Zone Layers
    if geometry.total_layers > 0
        rate = bz_cfg.threshold_rate;
        for layer = 1:geometry.total_layers
            % Progression 0..1
            progression = (layer / geometry.total_layers) ^ rate;

            % Interpolate
            layer_thresh = core_threshold + (1.0 - core_threshold) * progression;

            % Vectorized assignment using logical indexing
            mask_layer = (geometry.bz_map == layer);
            threshold_map(mask_layer) = layer_thresh;
        end
    end

    % 3. Apply Threshold
    presence = (noise_field >= threshold_map);

    % 4. Reshape presence and threshold_map to dimensional matrices
    target_dims = size(geometry.core_mask);
    presence = reshape(presence, target_dims);
    threshold_map = reshape(threshold_map, target_dims);

    % 5. Package Result
    result.presence = presence;
    result.bz_layers = geometry.bz_map;
    result.fc_density = actual_core_density;
    result.threshold_map = threshold_map;
end
