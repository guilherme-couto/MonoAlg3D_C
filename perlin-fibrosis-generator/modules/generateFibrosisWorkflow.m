function generateFibrosisWorkflow(domain_cfg, geom_cfg, bz_cfg, noise_params, run_cfg)
% Orchestrates the GNT (Geometry-Noise-Threshold) pipeline.

    %% 1. BUILD MESH
    % Assuming you have a buildMesh function that returns struct
    mesh = buildMesh(domain_cfg);

    %% 2. COMPUTE GEOMETRY (The "G" in GNT)
    % Runs ONCE. Heavily optimized with bwdist.
    geometry = computeFibrosisGeometry(mesh, geom_cfg, bz_cfg);

    if geometry.num_core_points == 0
        error('Geometry generation failed: Core is empty.');
    end

    %% 3. GENERATE PATTERN (The "N" and "T" in GNT)
    if strcmpi(run_cfg.mode, 'threshold')
        % Single Pass
        combined_noise = computeNoiseField(mesh, noise_params);

        final_result = applyFibrosisThreshold(combined_noise, geometry, run_cfg.target_density, bz_cfg);

    elseif strcmpi(run_cfg.mode, 'composition')
        % Iterative Loop (Uses the PRE-COMPUTED geometry)
        final_result = runCompositionLoop(mesh, geometry, noise_params, run_cfg.target_density, bz_cfg);

    else
        error('Unknown mode: %s', run_cfg.mode);
    end

    %% 4. SAVE RESULTS
    % No flips, no permutes. buildMesh geometry (Row 1 = Bottom).
    % The data is passed exactly as generated in geometry: (Ny, Nx, [Nz])
    % Structuring output for saving
    output_data.presence = final_result.presence;
    output_data.bz_layers = final_result.bz_layers;
    output_data.mesh = mesh;
    output_data.params = noise_params;
    output_data.is3D = mesh.is3D;

    if run_cfg.save_mesh
        writeALGFromPatternData(output_data, run_cfg.output_name);
    end

    if run_cfg.save_fig
        savePatternVisualisation(output_data, final_result, run_cfg);
    end
end
