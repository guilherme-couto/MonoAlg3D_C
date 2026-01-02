function geometry = computeFibrosisGeometry(mesh, geom_cfg, bz_cfg)
% Defines WHERE the fibrosis core and border zone layers are located
% Supports both 2D and 3D geometries based on mesh.Nz

    Nx = mesh.Nx;
    Ny = mesh.Ny;
    Nz = mesh.Nz;
    is3D = mesh.is3D;

    % Initialize Core Mask with appropriate dimensions
    if is3D
        core_mask = false(Ny, Nx, Nz);
    else
        core_mask = false(Ny, Nx);
    end

    % Dimensions of the fibrotic core in PIXELS for embedded shapes
    % Common pixel conversion for X and Y (used in almost all shapes)
    if isfield(geom_cfg, 'fc_width'), w_px = geom_cfg.fc_width / mesh.dx; end
    if isfield(geom_cfg, 'fc_height'), h_px = geom_cfg.fc_height / mesh.dx; end
    if isfield(geom_cfg, 'fc_center_x'), cx_px = geom_cfg.fc_center_x / mesh.dx; end
    if isfield(geom_cfg, 'fc_center_y'), cy_px = geom_cfg.fc_center_y / mesh.dx; end

    % --- DEFINE CORE SHAPE ---
    switch lower(geom_cfg.shape)

        case 'full'
            % Works for both 2D and 3D
            if is3D
                core_mask = true(Ny, Nx, Nz);
            else
                core_mask = true(Ny, Nx);
            end

        % === 2D SHAPES ===
        case 'rectangle'
            if is3D, error('Shape "rectangle" is for 2D meshes only. Use "box" for 3D.'); end

            [X, Y] = meshgrid(1:Nx, 1:Ny);
            x_min = cx_px - w_px/2; x_max = cx_px + w_px/2;
            y_min = cy_px - h_px/2; y_max = cy_px + h_px/2;

            core_mask = (X >= x_min & X <= x_max & Y >= y_min & Y <= y_max);

        case 'ellipse'
            if is3D, error('Shape "ellipse" is for 2D meshes only. Use "ellipsoid" for 3D.'); end

            [X, Y] = meshgrid(1:Nx, 1:Ny);
            a = w_px/2;
            b = h_px/2;

            % Normalized ellipse equation: ((x-cx)/a)^2 + ((y-cy)/b)^2 <= 1
            dist_sq = ((X - cx_px).^2) / a^2 + ((Y - cy_px).^2) / b^2;
            core_mask = (dist_sq <= 1);

        % === 3D SHAPES ===
        case 'box'
            if ~is3D, error('Shape "box" requires a 3D mesh (Lz > 0).'); end

            % Z-Dimension handling
            d_px = getDepthPx(geom_cfg, mesh, Nz);
            cz_px = getCenterZPx(geom_cfg, mesh, Nz);

            [X, Y, Z] = meshgrid(1:Nx, 1:Ny, 1:Nz);

            x_min = cx_px - w_px/2; x_max = cx_px + w_px/2;
            y_min = cy_px - h_px/2; y_max = cy_px + h_px/2;
            z_min = cz_px - d_px/2; z_max = cz_px + d_px/2;

            core_mask = (X >= x_min & X <= x_max & ...
                         Y >= y_min & Y <= y_max & ...
                         Z >= z_min & Z <= z_max);

        case 'ellipsoid'
            if ~is3D, error('Shape "ellipsoid" requires a 3D mesh (Lz > 0).'); end

            % Z-Dimension handling
            d_px = getDepthPx(geom_cfg, mesh, Nz);
            cz_px = getCenterZPx(geom_cfg, mesh, Nz);

            [X, Y, Z] = meshgrid(1:Nx, 1:Ny, 1:Nz);

            a = w_px/2;
            b = h_px/2;
            c = d_px/2;
            if c == 0, c = 1e-6; end % Prevent division by zero

            % Normalized ellipsoid equation: ((x-cx)/a)^2 + ((y-cy)/b)^2 + ((z-cz)/c)^2 <= 1
            dist_sq = ((X - cx_px).^2) / a^2 + ...
                      ((Y - cy_px).^2) / b^2 + ...
                      ((Z - cz_px).^2) / c^2;
            core_mask = (dist_sq <= 1);

        case 'stain'
            error('Stain shape logic needs migration.');

        otherwise
            error('Unknown shape: %s', geom_cfg.shape);
    end

    % --- BORDER ZONE CALCULATION (Dimension Agnostic) ---
    % Using Distance Transform
    if bz_cfg.active
        % bwdist calculates distance to nearest NON-ZERO pixel.
        % We want distance FROM the core (where core is 1).
        % So we compute distance on the INVERSE of the core.
        % 'chessboard' = 8-connectivity.
        % 'euclidean' = smooth circles.
        dist_map = bwdist(core_mask, bz_cfg.metric);

        % Convert distance to integer layers
        % Distance 0 = Inside Core
        % Distance 1 = Immediate neighbor (Layer 1)
        layer_map = double(dist_map);

        % Cap at max layers
        if isfield(bz_cfg, 'max_layers') && ~isinf(bz_cfg.max_layers)
            mask_valid_bz = (layer_map > 0) & (layer_map <= bz_cfg.max_layers);
        else
            mask_valid_bz = (layer_map > 0);
        end

        % Create Final Maps
        % Convention: 0 = Core, 1..N = BZ, -1 = Healthy
        final_bz_map = -ones(size(core_mask));
        final_bz_map(core_mask) = 0;
        final_bz_map(mask_valid_bz) = layer_map(mask_valid_bz);

        total_layers = max(final_bz_map(:));
    else
        % No BZ logic
        final_bz_map = -ones(size(core_mask));
        final_bz_map(core_mask) = 0;
        total_layers = 0;
    end

    % --- OUTPUT STRUCT ---
    geometry.core_mask = core_mask;       % Logical Matrix
    geometry.core_indices = find(core_mask);
    geometry.num_core_points = length(geometry.core_indices);
    geometry.bz_map = final_bz_map;       % Matrix (-1, 0, 1..N)
    geometry.total_layers = total_layers;
    geometry.is3D = is3D;
end

% --- Helper Functions for 3D Params ---
function d_px = getDepthPx(geom_cfg, mesh, Nz)
    if isfield(geom_cfg, 'fc_depth')
        d_px = geom_cfg.fc_depth / mesh.dx;
    else
        d_px = Nz; % Default to full depth
    end
end

function cz_px = getCenterZPx(geom_cfg, mesh, Nz)
    if isfield(geom_cfg, 'fc_center_z')
        cz_px = geom_cfg.fc_center_z / mesh.dx;
    else
        cz_px = Nz / 2; % Default to center
    end
end
