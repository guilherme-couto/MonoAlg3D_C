function final_result = runCompositionLoop(mesh, geometry, noise_params, target_density, bz_cfg)
% Iteratively adds patterns until density is reached.
% Critical: It reuses 'geometry', avoiding recalculation of BZ layers.

    tolerance = 0.005;
    target = target_density;

    Nx = mesh.Nx;
    Ny = mesh.Ny;
    Nz = mesh.Nz;
    is3D = mesh.is3D;

    % Initialize
    current_density = 0;
    if is3D
        accumulated_presence = false(Ny, Nx, Nz);
    else
        accumulated_presence = false(Ny, Nx);
    end

    % Initial Threshold logic for composition
    req_density = target;

    iter = 0;
    max_tries = 10;
    tries = 0;

    while abs(current_density - target) > tolerance
        iter = iter + 1;

        % 1. Generate Noise (Fast)
        seed_iter = noise_params.seed + iter;
        noise_field = computeNoiseField(mesh, noise_params, seed_iter);

        % 2. Apply Threshold (Reuses geometry)
        if req_density >= 0.1 - tolerance
             eff_thresh_dens = 0.1;
        else
             eff_thresh_dens = req_density;
        end

        step_result = applyFibrosisThreshold(noise_field, geometry, eff_thresh_dens, bz_cfg);

        % 3. Accumulate
        temp_presence = accumulated_presence | step_result.presence;

        % 4. Check Density (Only inside Core)
        temp_density = sum(temp_presence(geometry.core_indices)) / geometry.num_core_points;

        tries = tries + 1;

        % Logic to accept or reject
        if temp_density < target + tolerance
            % Accept
            accumulated_presence = temp_presence;
            current_density = temp_density;
            tries = 0;

            % Update requirement
            req_density = abs(current_density - target);
        end

        % Reset if stuck
        if tries > max_tries
            if is3D
                accumulated_presence = false(Ny, Nx, Nz);
            else
                accumulated_presence = false(Ny, Nx);
            end
            current_density = 0;
            req_density = target;
            tries = 0;
        end

        if iter > 200
            warning('Composition loop limit reached.');
            break;
        end
    end

    final_result.presence = accumulated_presence;
    final_result.bz_layers = geometry.bz_map;
    final_result.fc_density = current_density;
end
