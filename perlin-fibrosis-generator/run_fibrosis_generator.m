function run_fibrosis_generator(fibrosis_type, density, seed, angle_deg_or_vec, dimension_mode, domain_dims, shape, core_dims, output_filename, save_mesh, save_figure)
%RUN_FIBROSIS_GENERATOR Entry point for batch execution via Python/CLI.
%   Dependencies: 'image' package must be installed in Octave.
%   (run 'pkg install -forge image' in Octave console once before using).
%
%   Inputs:
%       fibrosis_type    : String ('compact', 'interstitial', 'diffuse', 'patchy')
%       density          : Scalar or String (Target density, e.g., 0.4)
%       seed             : Scalar or String (RNG seed)
%       angle_deg_or_vec : Scalar, Vector, or String (Fiber angle in DEGREES or vector of angles in DEGREES for 3D)
%       dimension_mode   : String ('2D', '3D', 'custom')
%       domain_dims      : Vector [dx, Lx, Ly, (Lz)] (cm)
%       shape            : String (2D/3D: 'full', 2D: 'ellipse', 2D: 'rectangle', 3D: 'box', 3D: 'ellipsoid')
%       core_dims        : Vector [width, height, (depth)] (cm) - Ignored if shape='full'
%       output_filename  : String (Path/Filename without extension)
%       save_mesh        : Boolean (Whether to save mesh file)
%       save_figure      : Boolean (Whether to save figure file)

    %% 0. OCTAVE SETUP
    addpath(pwd);
    addpath(fullfile(pwd, 'modules'));
    addpath(fullfile(pwd, 'helpers'));

    % Load the image package explicitly.
    try pkg load image; catch error('Package "image" not found. Please install it running: pkg install -forge image'); end
    more off; % Forces command line output without pagination

    %% 1. INPUT SANITIZATION
    % Octave CLI calls might pass numbers as strings.
    if ischar(density) || isstring(density), density = str2double(density); end
    if ischar(seed) || isstring(seed), seed = str2double(seed); end

    angle_deg_or_vec = parseVectorInput(angle_deg_or_vec);
    domain_dims      = parseVectorInput(domain_dims);
    core_dims        = parseVectorInput(core_dims);

    % Normalize dimension mode string
    dimension_mode = upper(char(dimension_mode)); % '2D', '3D', 'CUSTOM'

    %% 2. FIBROSIS DATABASE DEFINITION
    % Order: [fibreness, fibre_sep, patchiness, feature_size, roughness, patch_size, alignment, orientation]
    switch lower(fibrosis_type)
        case 'compact',      base_params = [NaN, NaN, 0.44, 0.96, 0.59, 2.03, 2.47, 0]; mode = 'threshold';
        case 'interstitial', base_params = [0.30, 0.31, 0.32, 0.24, 0.96, 4.67, 1.89, 0]; mode = 'composition';
        case 'diffuse',      base_params = [NaN, NaN, 0.49, 0.07, 0.44, 1.22, 2.17, 0]; mode = 'composition';
        case 'patchy',       base_params = [0.38, 0.31, 0.42, 0.32, 0.78, 2.10, 2.50, 0]; mode = 'composition';
        otherwise, error('Unknown fibrosis type: %s', fibrosis_type);
    end

    %% 3. BUILD CONFIGURATION STRUCTS

    % 1. Domain Configuration
    % Expecting domain_dims = [dx, Lx, Ly, (Lz)]
    if length(domain_dims) < 3
        error('domain_dims must have at least 3 elements: [dx, Lx, Ly, (Lz)]');
    end

    domain_config.dx = domain_dims(1); % Grid Spacing (cm)
    domain_config.Lx = domain_dims(2); % Width (cm)
    domain_config.Ly = domain_dims(3); % Height (cm)
    if length(domain_dims) == 4
        domain_config.Lz = domain_dims(4); % Depth (cm)
    end

    if strcmp(dimension_mode, '3D')
        if domain_config.Lz <= 0
            warning('Dimension mode is 3D but Lz is 0. Forcing Lz=dx.');
            domain_config.Lz = domain_config.dx;
        end
        is3D = true;
    elseif strcmp(dimension_mode, '2D')
        domain_config.Lz = 0.0; % Enforce 0 for 2D
        is3D = false;
    elseif strcmp(dimension_mode, 'CUSTOM')
        % Future logic for patient specific meshes
        % domain_config.mesh_file = 'path/to/mesh.vtk';
        error('Custom mode not yet implemented.');
    else
        error('Invalid dimension_mode. Use "2D" or "3D".');
    end

    % 2. Fibrotic Core Geometry Configuration (Where the core is)
    geom_config.shape = shape;

    % For embedded shapes, define core dimensions
    if ~strcmpi(shape, 'full')
        % Expecting core_dims = [width, height, (depth)]
        if length(core_dims) < 2
            error('core_dims must have at least 2 elements: [width, height, (depth)]');
        end

        geom_config.fc_width  = core_dims(1); % Fibrotic Core Width (cm)
        geom_config.fc_height = core_dims(2); % Fibrotic Core Height (cm)
        if length(core_dims) == 3
            geom_config.fc_depth  = core_dims(3); % Fibrotic Core Depth (cm)
        end

        % Check if fibrotic core dimensions for embedded shapes fit within domain
        if geom_config.fc_width > domain_config.Lx || geom_config.fc_height > domain_config.Ly
            error('Fibrotic core dimensions exceed domain size.');
        end
        if is3D && geom_config.fc_depth > domain_config.Lz
            error('Fibrotic core depth exceeds domain depth.');
        end

        % Auto-center logic (Always center in the domain)
        geom_config.fc_center_x = domain_config.Lx / 2;
        geom_config.fc_center_y = domain_config.Ly / 2;
        geom_config.fc_center_z = domain_config.Lz / 2;

        % Check if geometry shape is compatible with dimension mode
        if is3D && any(strcmpi(shape, {'rectangle', 'ellipse'}))
            error('Shape "%s" incompatible with 3D. Use "box" or "ellipsoid".', shape);
        elseif ~is3D && any(strcmpi(shape, {'box', 'ellipsoid'}))
            error('Shape "%s" incompatible with 2D. Use "rectangle" or "ellipse".', shape);
        end
    end

    % 3. Border Zone Configuration
    % Logic: 'full' = No BZ
    % 'rectangle'/'box'/'ellipse'/'ellipsoid' = Yes BZ
    if strcmpi(shape, 'full')
        bz_config.active = false;
        bz_config.max_layers = 0;
    else
        bz_config.active = true;
        bz_config.threshold_rate = 1.0;  % Falloff rate (1=linear)
        bz_config.metric = 'chessboard'; % 'chessboard' (8-conn) or 'euclidean' (smooth)

        % Number of layers (10% of smallest core dimension)
        if is3D
            bz_config.bz_thickness = 0.1 * min([geom_config.fc_width, geom_config.fc_height, geom_config.fc_depth]);
            bz_config.max_layers = ceil(bz_config.bz_thickness / domain_config.dx);
        else
            bz_config.bz_thickness = 0.1 * min([geom_config.fc_width, geom_config.fc_height]);
            bz_config.max_layers = ceil(bz_config.bz_thickness / domain_config.dx);
        end
    end

    % 4. Noise Parameters (The texture)
    noise_params = vectorToParams(base_params, angle_deg_or_vec, is3D);
    noise_params.variable_direction = false; % Enable fiber direction variability
    noise_params.seed = seed; % RNG Seed

    % 5. Run Configuration
    run_config.fibrosis_type = fibrosis_type;
    run_config.mode = mode;
    run_config.target_density = density;
    run_config.seed = seed;
    run_config.output_name = output_filename;
    run_config.save_mesh = save_mesh;
    run_config.save_fig = save_figure;

    % Check if the output directory exists. if not, create it
    [output_dir, ~, ~] = fileparts(output_filename);
    if ~isempty(output_dir)
        [~, ~, ~] = mkdir(output_dir);
    end

    % === EXECUTE WORKFLOW ===
    try
        generateFibrosisWorkflow(domain_config, geom_config, bz_config, noise_params, run_config);
    catch ME
        error('Fibrosis generation failed: %s', ME.message);
        error('Traceback: %s', ME.stack(1).name);
    end
end

% --- Helpers ---
function val = parseVectorInput(input_val)
    % Handles: Scalar, Vector, or String representation "[1 2 3]"
    if ischar(input_val) || isstring(input_val)
        % Remove brackets and convert to numbers
        clean_str = strrep(strrep(input_val, '[', ''), ']', '');
        val = str2num(clean_str);
    else
        val = input_val;
    end
end

function params_struct = vectorToParams(vec, angle_deg_or_vec, is3D)
    % Maps input vector to struct. Handles 2D (8 params) and 3D (11 params)

    if is3D && length(vec) == 8
        % Convert 2D params to 3D equivalent
        vec = convert_params_2D_to_3D(vec);
    end

    % Common parameters
    params_struct.fibreness     = vec(1);
    params_struct.fibre_sep_y   = vec(2);
    params_struct.patchiness    = vec(3);
    params_struct.feature_size  = vec(4);
    params_struct.roughness     = vec(5);
    params_struct.patch_size    = vec(6);
    params_struct.alignment_y   = vec(7);

    % Check if is 3D
    if is3D
        % 3D Specifics
        params_struct.fibre_sep_z   = vec(8);
        params_struct.alignment_z   = vec(9);
        if length(angle_deg_or_vec) == 2
            params_struct.phi           = deg2rad(angle_deg_or_vec(1));
            params_struct.theta         = deg2rad(angle_deg_or_vec(2));
        elseif isscalar(angle_deg_or_vec)
            % If only one angle is provided, phi=theta=angle
            params_struct.phi           = deg2rad(angle_deg_or_vec);
            params_struct.theta         = deg2rad(angle_deg_or_vec);
        else
            error('For 3D, angle_deg_or_vec must be a scalar or a 2-element vector.');
        end

        params_struct.orientation = 0; % Not used in 3D logic
    else
        % 2D Defaults for compatibility
        params_struct.fibre_sep_z   = 1;
        params_struct.alignment_z   = 1;
        params_struct.phi           = 0;
        params_struct.theta         = 0;

        % In 2D, the external 'angle_deg' controls orientation
        params_struct.orientation = deg2rad(angle_deg_or_vec);
    end

    % Alias for compatibility with 2D code
    params_struct.fibre_sep = params_struct.fibre_sep_y;
    params_struct.alignment = params_struct.alignment_y;
end

function params3D = convert_params_2D_to_3D(params2D)
    % Convert parameters as needed in moving from 2D to 3D.
    % 2D: [fibreness, fibre_sep, patchiness, feature_size, roughness, patch_size, alignment, orientation]
    % 3D: [fibreness, fibre_sep_y, patchiness, feature_size, roughness, patch_size, alignment_y, fibre_sep_z, alignment_z, phi, theta]

    % Base parameters
    fibreness     = params2D(1);
    fibre_sep_y   = params2D(2);
    patchiness    = params2D(3);
    feature_size  = params2D(4);
    roughness     = params2D(5);
    patch_size    = params2D(6);
    alignment_y   = params2D(7);

    % Specify the modifier to anisotropy in the third dimension
    K = 1; % no modification

    % Convert parameters as needed in moving from 2D to 3D
    fibre_sep_z = fibre_sep_y;
    feature_size = feature_size .* alignment_y.^(-1/6) * K^(-1/3);
    alignment_z = alignment_y * K;

    phi = 0;    % Default angle in radians
    theta = 0;  % Default angle in radians

    % Convert from the 2D learned parameters to a 3D equivalent
    params3D = [fibreness, fibre_sep_y, patchiness, feature_size, roughness, patch_size, alignment_y, fibre_sep_z, alignment_z, phi, theta];
end
