function savePatternVisualisation(output_data, final_result, run_cfg)
% Saves a figure visualizing the fibrosis pattern (Supports 2D and 3D).

    % Extract Data
    % We pass the full 3D matrix (if 3D) to Python.
    presence = output_data.presence;

    output_name = run_cfg.output_name;
    fibrosis_type = run_cfg.fibrosis_type;
    density = run_cfg.target_density;
    seed = run_cfg.seed;

    % HEADLESS FIGURE SAVING (Via Python)
    % 1. Check if file already exists to avoid permission locks
    temp_mat_file = [output_name, '_temp.mat'];

    % 2. Save the FULL matrix to .mat
    % '-v7' ensures compatibility with SciPy
    save("-v7", temp_mat_file, "presence");

    % 3. Format Title
    if isfield(output_data, 'is3D') && output_data.is3D
        dim_tag = '3D';
    else
        dim_tag = '2D';
    end
    title_str = sprintf('%s (%s) | Density: %.2f | Seed: %d', fibrosis_type, dim_tag, density, seed);

    % 4. Call Python
    fprintf('Delegating %s visualization to Python...\n', dim_tag);
    cmd = sprintf("python save_fibrosis_image.py '%s' '%s' '%s'", temp_mat_file, output_name, title_str);

    [status, cmd_output] = system(cmd);

    if status ~= 0
        fprintf('Error calling Python bridge:\n%s\n', cmd_output);
    else
        fprintf('%s', cmd_output);
    end
end
