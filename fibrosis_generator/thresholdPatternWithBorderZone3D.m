function [presence_with_border_zone, fc_density, bz_layers, fibrotic_core, total_layers, threshold_map] = thresholdPatternWithBorderZone3D(field, density, mesh, params)
% thresholdPatternWithBorderZone - Applies a spatially varying threshold to a noise field.
%
% This function defines a central "Fibrotic Core" (FC) and a surrounding
% "Border Zone" (BZ). The fibrosis density is met within the FC, and a
% progressively increasing threshold is applied to the layers of the BZ
% to create a natural falloff effect.
%
% INPUTS:
%   field   - The 1D noise field [1 x (Nx*Ny*Nz)] with values typically in [0, 1].
%   density - The target fibrosis density to be achieved within the Fibrotic Core.
%   mesh    - A struct with mesh dimensions (Nx, Ny, Nz) and spacing (dx).
%   params  - A struct containing control parameters:
%     .core_shape          - Shape of the FC: 'cuboid', 'ellipsoid', 'elliptic_cylinder' or 'stain'/'stainoid'.
%     --- Shape Control (use one method) ---
%     .core_volume_fraction  - (Method 1) Desired volume of the FC as a fraction [0, 1].
%     .core_width            - (Method 2) Physical width of the FC (e.g., in cm).
%     .core_height           - (Method 2) Physical height of the FC (e.g., in cm).
%     .core_depth            - (Method 2) Physical depth of the FC (e.g., in cm).
%     --- Border Zone Control ---
%     .max_layers            - (Optional) Max number of layers in the BZ. Default: fill domain.
%     .threshold_rate        - (Optional) Exponent for threshold increase in BZ. 1=linear, 2=quadratic. Default: 1.
%     --- Blob Shape Specific ---
%     .stain_seed             - (Optional) Seed for 'stain' shape generation. Default: 1.
%     .stain_scale            - (Optional) Scale for 'stain' shape noise. Default: 2.
%     --- Visualization ---
%     .save_figure           - (Optional) Boolean, true to save a figure of the output. Default: true.
%     .figure_filename       - (Optional) String, filename for the saved figure. Default: 'debug_threshold.png'.
%
% OUTPUTS:
%   presence_with_border_zone - The final binary fibrosis pattern [1 x (Nx*Ny*Nz)].
%   fc_density                - The actual fibrosis density achieved within the FC.
%   bz_layers                 - An integer map where each element is its BZ layer number (FC=0, unassigned=-1).
%   fibrotic_core             - A binary mask indicating the FC (1) and BZ (0).
%   total_layers              - The total number of BZ layers generated.
%   threshold_map             - A map of the final threshold value applied to each element.

% --- Parameter Validation and Defaults ---
if ~isfield(params, 'max_layers'), params.max_layers = Inf; end
if ~isfield(params, 'threshold_rate'), params.threshold_rate = 1; end
if ~isfield(params, 'stain_seed'), params.stain_seed = 1; end
if ~isfield(params, 'stain_scale'), params.stain_scale = 2; end
if ~isfield(params, 'save_figure'), params.save_figure = true; end
if ~isfield(params, 'figure_filename'), params.figure_filename = 'debug_threshold.png'; end

Nx = mesh.Nx;
Ny = mesh.Ny;
Nz = mesh.Nz;

% =========================================================================
% STEP 1: Determine Fibrotic Core (FC) and Border Zone (BZ)
% =========================================================================
fibrotic_core = createFibroticCore3D(mesh, params);
fc_indices = find(fibrotic_core == 1);
num_fc_points = length(fc_indices);
if num_fc_points == 0
    error('Fibrotic Core generation resulted in an empty set. Try adjusting core_volume_fraction or stain parameters.');
end

% =========================================================================
% STEP 2: Discretize the BZ into layers
% =========================================================================
bz_layers = -ones(Ny, Nx, Nz); % Initialize with -1 (unassigned)
bz_layers(fc_indices) = 0; % FC is layer 0

current_layer = 1;
points_assigned_in_last_step = true;

% Propagate layers outwards from the core
while points_assigned_in_last_step && current_layer <= params.max_layers
    points_assigned_in_last_step = false;

    % Get indices of the previous layer's frontier
    frontier_indices = find(bz_layers == current_layer - 1);

    % Find unassigned neighbors of the frontier
    for idx = frontier_indices'
        [r, c, p] = ind2sub([Ny, Nx, Nz], idx); % r - row, c - column, p - page

        % Check 26-connectivity neighbors (including diagonals)
        neighbors = [r-1, c,   p;    % Top
                     r+1, c,   p;    % Bottom
                     r,   c-1, p;    % Left
                     r,   c+1, p;    % Right
                     r,   c,   p+1;  % Front
                     r,   c,   p-1;  % Back
                     r-1, c-1, p;    % Top-Left
                     r-1, c+1, p;    % Top-Right
                     r+1, c-1, p;    % Bottom-Left
                     r+1, c+1, p;    % Bottom-Right
                     r-1, c,   p+1;  % Front-Top
                     r+1, c,   p+1;  % Front-Bottom
                     r,   c-1, p+1;  % Front-Left
                     r,   c+1, p+1;  % Front-Right
                     r-1, c,   p-1;  % Back-Top
                     r+1, c,   p-1;  % Back-Bottom
                     r,   c-1, p-1;  % Back-Left
                     r,   c+1, p-1;  % Back-Right
                     r-1, c-1, p+1;  % Front-Top-Left
                     r-1, c+1, p+1;  % Front-Top-Right
                     r+1, c-1, p+1;  % Front-Bottom-Left
                     r+1, c+1, p+1;  % Front-Bottom-Right
                     r-1, c-1, p-1;  % Back-Top-Left
                     r-1, c+1, p-1;  % Back-Top-Right
                     r+1, c-1, p-1;  % Back-Bottom-Left
                     r+1, c+1, p-1]; % Back-Bottom-Right


        for n = 1:26
            nr = neighbors(n, 1);
            nc = neighbors(n, 2);
            np = neighbors(n, 3);

            % Check if neighbor is within bounds
            if nr >= 1 && nr <= Ny && nc >= 1 && nc <= Nx && np >= 1 && np <= Nz
                if bz_layers(nr, nc, np) == -1 % If unassigned
                    bz_layers(nr, nc, np) = current_layer;
                    points_assigned_in_last_step = true;
                end
            end
        end
    end

    if points_assigned_in_last_step
        current_layer = current_layer + 1;
    end
end

total_layers = max(bz_layers(:));
bz_layers = reshape(bz_layers, 1, []); % Reshape to 1D array for output

% =========================================================================
% STEP 3: Determine the threshold for the FC (threshold_fc)
% =========================================================================
% Bisection method to find the threshold that yields the target density *within the FC*
max_iters = 40;
tol = 1e-4;
window_threshold = [max(field); min(field)];
iters = 0;
err = 1;

while iters < max_iters && err > tol
    check_threshold = mean(window_threshold);

    % Calculate density only on points within the Fibrotic Core
    check_density = sum(field(fc_indices) >= check_threshold) / num_fc_points;

    if check_density < density
        window_threshold(1) = check_threshold;
    else
        window_threshold(2) = check_threshold;
    end

    err = abs(density - check_density);
    iters = iters + 1;
end
fc_density = check_density; % Actual density achieved in the FC
threshold_fc = mean(window_threshold);

% =========================================================================
% STEP 4: Define the threshold for each BZ layer
% =========================================================================
threshold_map = ones(1, Nx * Ny * Nz) + 2; % Initialize with a high value
threshold_map(fc_indices) = threshold_fc;

if total_layers > 0
    % Calculate the threshold for each layer, increasing up to >= 1.0
    for i = 1:total_layers
        layer_indices = find(bz_layers == i);

        % Calculate progression factor (0 to 1) based on layer number and rate
        progression = (i / total_layers) ^ params.threshold_rate;

        % Interpolate threshold between threshold_fc and 1.0
        threshold_layer = threshold_fc + (1.0 - threshold_fc) * progression;

        threshold_map(layer_indices) = threshold_layer;
    end
end
% Any points remaining at -1 (if max_layers was reached) will keep the high
% initial threshold, effectively preventing fibrosis there.

% =========================================================================
% STEP 5: Apply the threshold_map to the field
% =========================================================================
presence_with_border_zone = double(field >= threshold_map);

% Proportion of fibrosis inside the FC over total fibrosis
fibrosis_count_inside_fc = sum(fibrotic_core(:) & presence_with_border_zone(:));
fibrosis_count_outside_fc = sum(~fibrotic_core(:) & presence_with_border_zone(:));
fc_proportion = fibrosis_count_inside_fc / (fibrosis_count_inside_fc + fibrosis_count_outside_fc);
presence_with_categories = presence_with_border_zone + (fibrotic_core & presence_with_border_zone); % 0 = Healthy, 1 = FC Fibrosis, 2 = BZ Fibrosis


% =========================================================================
% STEP 6: (Optional) Visualize and Save Figure
% =========================================================================
% TODO
% if params.save_figure
%     % === Debug plots ===
%     figure('Visible', 'off');

%     subplot(3,2,1);
%     imagesc(flipud(reshape(fibrotic_core', Nx, Ny)')); axis image off; colorbar;
%     title(sprintf('Fibrotic Core (%s a=%.2g%%)', params.core_shape, sum(fibrotic_core(:)) * 100 / (Nx*Ny)));

%     subplot(3,2,2);
%     imagesc(flipud(reshape(bz_layers', Nx, Ny)')); axis image off; colorbar;
%     title(sprintf('Border Zone Layers (total=%d)', total_layers));

%     subplot(3,2,3);
%     imagesc(flipud(reshape(field', Nx, Ny)')); axis image off; colorbar;
%     title(sprintf('Noise Field'));

%     subplot(3,2,4);
%     imagesc(flipud(reshape(threshold_map', Nx, Ny)')); axis image off; colorbar; caxis([threshold_fc max(field)]);
%     title(sprintf('Threshold Map (rate=%.1f)', params.threshold_rate));

%     subplot(3,2,5);
%     imagesc(flipud(reshape(presence_with_categories', Nx, Ny)')); axis image off; colorbar;
%     title(sprintf('FC and BZ (FCd=%.4g%% | FCp=%.2g%%)', fc_density*100, fc_proportion*100));

%     subplot(3,2,6);
%     imagesc(flipud(reshape(presence_with_border_zone', Nx, Ny)')); axis image off; colorbar;
%     title(sprintf('Presence'));

%     print(params.figure_filename, '-dpng', '-r300');
% end

end

% --- Helper Function to Create the 3D Fibrotic Core ---
function fibrotic_core = createFibroticCore3D(mesh, params)
% This helper function creates the 3D binary mask for the Fibrotic Core (FC).

Nx = mesh.Nx;
Ny = mesh.Ny;
Nz = mesh.Nz;
dx = mesh.dx; % Physical size of one element
Lx = Nx * dx; % Total physical width
Ly = Ny * dx; % Total physical height
Lz = Nz * dx; % Total physical depth

fibrotic_core = zeros(Ny, Nx, Nz);

% --- Determine if using direct dimensions or volume fraction ---
use_direct_dims = isfield(params, 'core_width') && isfield(params, 'core_height') && isfield(params, 'core_depth');

switch lower(params.core_shape)
    case {'cuboid', 'ellipsoid', 'elliptic_cylinder'}
        w_pixels = 0; h_pixels = 0; d_pixels = 0;

        if use_direct_dims
            % --- METHOD 2: Use direct physical dimensions ---
            if params.core_width > Lx || params.core_height > Ly || params.core_depth > Lz
                error('Requested core dimensions (%.2f x %.2f x %.2f) exceed domain limits (%.2f x %.2f x %.2f).', ...
                      params.core_width, params.core_height, params.core_depth, Lx, Ly, Lz);
            end
            % CONVENTION: width->x, height->y, depth->z
            w_pixels = params.core_width / dx;  % x-dimension (cols)
            h_pixels = params.core_height / dx; % y-dimension (rows)
            d_pixels = params.core_depth / dx;  % z-dimension (pages)
        else
            % --- METHOD 1: Use volume fraction (fallback) ---
            if ~isfield(params, 'core_volume_fraction')
                error("For 'cuboid' or 'ellipsoid', you must provide either core dimensions (core_width, core_height) or core_volume_fraction.");
            end
            target_volume_pixels = params.core_volume_fraction * (Nx * Ny * Nz);
            aspect_ratio_yx = Ny / Nx;
            aspect_ratio_zx = Nz / Nx;

            if strcmpi(params.core_shape, 'cuboid')
                % For cuboid: V_cuboid = w*h*d = w * (w*ar_yx) * (w*ar_zx) = w^3 * ar_yx * ar_zx
                w_cuboid = (target_volume_pixels / (aspect_ratio_yx * aspect_ratio_zx))^(1/3);
                w_pixels = w_cuboid;
            elseif strcmpi(params.core_shape, 'ellipsoid')
                % For ellipsoid: V_ellipsoid = (pi/6) * V_cuboid => V_cuboid = V_ellipsoid * (6/pi)
                effective_cuboid_volume = target_volume_pixels * 6 / pi;
                w_pixels = (effective_cuboid_volume / (aspect_ratio_yx * aspect_ratio_zx))^(1/3);
            elseif strcmpi(params.core_shape, 'elliptic_cylinder')
                % For elliptic_cylinder: V_cylinder = Area_base * depth = (pi/4)*w*h*d
                % V = (pi/4) * w * (w*ar_yx) * (w*ar_zx) = (pi/4) * w^3 * ar_yx * ar_zx
                effective_cuboid_volume = target_volume_pixels * 4 / pi;
                w_pixels = (effective_cuboid_volume / (aspect_ratio_yx * aspect_ratio_zx))^(1/3);
            end

            h_pixels = w_pixels * aspect_ratio_yx;
            d_pixels = w_pixels * aspect_ratio_zx;
        end

        if strcmpi(params.core_shape, 'cuboid')
            x1 = floor((Nx - w_pixels) / 2) + 1;
            x2 = ceil((Nx + w_pixels) / 2);
            y1 = floor((Ny - h_pixels) / 2) + 1;
            y2 = ceil((Ny + h_pixels) / 2);
            z1 = floor((Nz - d_pixels) / 2) + 1;
            z2 = ceil((Nz + d_pixels) / 2);
            fibrotic_core(y1:y2, x1:x2, z1:z2) = 1;
        elseif strcmpi(params.core_shape, 'ellipsoid')
            [X, Y, Z] = meshgrid(1:Nx, 1:Ny, 1:Nz);
            cx = Nx / 2;
            cy = Ny / 2;
            cz = Nz / 2;
            a = w_pixels / 2; % semi-axis a
            b = h_pixels / 2; % semi-axis b
            c = d_pixels / 2; % semi-axis c
            ellipsoid_eq = ((X - cx).^2 / a^2) + ((Y - cy).^2 / b^2) + ((Z - cz).^2 / c^2);
            fibrotic_core(ellipsoid_eq <= 1) = 1;
        elseif strcmpi(params.core_shape, 'elliptic_cylinder')
            % 1. Create the 2D elliptical base
            [X, Y] = meshgrid(1:Nx, 1:Ny);
            cx = Nx / 2; cy = Ny / 2;
            a = w_pixels / 2; b = h_pixels / 2;
            ellipse_base = zeros(Ny, Nx);
            if a > 0 && b > 0
                ellipse_eq = ((X - cx).^2 / a^2) + ((Y - cy).^2 / b^2);
                ellipse_base(ellipse_eq <= 1) = 1;
            end

            % 2. Extrude the base along the Z-axis
            z1 = floor((Nz - d_pixels) / 2) + 1;
            z2 = floor(z1 + d_pixels - 1);

            % Replicate the base for each slice in the depth range
            for z_idx = z1:z2
                fibrotic_core(:, :, z_idx) = ellipse_base;
            end
        end

    case {'stain', 'stainoid'}
        if ~isfield(params, 'core_volume_fraction')
            error("For 'stain' shape, you must provide 'core_volume_fraction'.");
        end
        target_volume_fraction = params.core_volume_fraction;

        [X_norm, Y_norm, Z_norm] = meshgrid(linspace(-1, 1, Nx), linspace(-1, 1, Ny), linspace(-1, 1, Nz));
        points = [X_norm(:)'; Y_norm(:)'; Z_norm(:)'];
        [perm_table, ~] = generateTables3D(params.stain_seed);
        stain_noise = Octave3D(points * params.stain_scale, 1, 0.5, perm_table, []);
        stain_noise = reshape(stain_noise, Ny, Nx, Nz);

        [X_grad, Y_grad, Z_grad] = meshgrid(linspace(-Nx/2, Nx/2, Nx), linspace(-Ny/2, Ny/2, Ny), linspace(-Nz/2, Nz/2, Nz));
        gradient_radius = sqrt((X_grad/(Nx/2)).^2 + (Y_grad/(Ny/2)).^2 + (Z_grad/(Nz/2)).^2);
        containment_gradient = max(0, 1 - gradient_radius.^2);
        stain_noise = stain_noise .* containment_gradient;

        window_thresh = [max(stain_noise(:)), min(stain_noise(:))];
        iters = 0; max_iters = 30; tol = 0.01;

        while iters < max_iters
            check_thresh = mean(window_thresh);
            binary_mask = stain_noise >= check_thresh;
            current_volume_fraction = sum(binary_mask(:)) / (Nx * Ny * Nz);

            err = abs(current_volume_fraction - target_volume_fraction);
            if err < tol, break; end

            if current_volume_fraction > target_volume_fraction
                window_thresh(2) = check_thresh;
            else
                window_thresh(1) = check_thresh;
            end
            iters = iters + 1;
        end
        fibrotic_core = stain_noise >= mean(window_thresh);

    otherwise
        error("Unknown core_shape specified. Use 'cuboid', 'ellipsoid', or 'stainoid'.");
end

% Reshape to 1D array to match 'field' format
fibrotic_core = reshape(fibrotic_core, 1, []);
end
