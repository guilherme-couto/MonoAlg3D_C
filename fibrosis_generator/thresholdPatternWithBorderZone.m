function [presence_with_border_zone, fc_density, bz_layers, fibrotic_core, total_layers, threshold_map] = thresholdPatternWithBorderZone(field, density, mesh, params)
% thresholdPatternWithBorderZone - Applies a spatially varying threshold to a noise field.
%
% This function defines a central "Fibrotic Core" (FC) and a surrounding
% "Border Zone" (BZ). The fibrosis density is met within the FC, and a
% progressively increasing threshold is applied to the layers of the BZ
% to create a natural falloff effect.
%
% INPUTS:
%   field   - The 1D noise field [1 x (Nx*Ny)] with values typically in [0, 1].
%   density - The target fibrosis density to be achieved within the Fibrotic Core.
%   mesh    - A struct with mesh dimensions (Nx, Ny) and spacing (dx).
%   params  - A struct containing control parameters:
%     .core_shape          - Shape of the FC: 'rectangle', 'ellipse', or 'stain'.
%     --- Shape Control (use one method) ---
%     .core_area_fraction  - (Method 1) Desired area of the FC as a fraction [0, 1].
%     .core_width          - (Method 2) Physical width of the FC (e.g., in cm).
%     .core_height         - (Method 2) Physical height of the FC (e.g., in cm).
%     --- Border Zone Control ---
%     .max_layers          - (Optional) Max number of layers in the BZ. Default: fill domain.
%     .threshold_rate      - (Optional) Exponent for threshold increase in BZ. 1=linear, 2=quadratic. Default: 1.
%     --- Stain Shape Specific ---
%     .stain_seed           - (Optional) Seed for 'stain' shape generation. Default: 1.
%     .stain_scale          - (Optional) Scale for 'stain' shape noise. Default: 2.
%     --- Visualization ---
%     .save_figure         - (Optional) Boolean, true to save a figure of the output. Default: true.
%     .figure_filename     - (Optional) String, filename for the saved figure. Default: 'debug_threshold.png'.
%
% OUTPUTS:
%   presence_with_border_zone - The final binary fibrosis pattern [1 x (Nx*Ny)].
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

% =========================================================================
% STEP 1: Determine Fibrotic Core (FC) and Border Zone (BZ)
% =========================================================================
fibrotic_core = createFibroticCore(mesh, params);
fc_indices = find(fibrotic_core == 1);
num_fc_points = length(fc_indices);
if num_fc_points == 0
    error('Fibrotic Core generation resulted in an empty set. Try adjusting core_area_fraction or stain parameters.');
end

% =========================================================================
% STEP 2: Discretize the BZ into layers
% =========================================================================
bz_layers = -ones(Ny, Nx); % Initialize with -1 (unassigned)
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
        [r, c] = ind2sub([Ny, Nx], idx);

        % Check 8-connectivity neighbors (including diagonals)
        neighbors = [r-1, c;    % Top
                     r+1, c;    % Bottom
                     r,   c-1;  % Left
                     r,   c+1;  % Right
                     r-1, c-1;  % Top-Left
                     r-1, c+1;  % Top-Right
                     r+1, c-1;  % Bottom-Left
                     r+1, c+1]; % Bottom-Right

        for n = 1:8
            nr = neighbors(n, 1);
            nc = neighbors(n, 2);

            % Check if neighbor is within bounds
            if nr >= 1 && nr <= Ny && nc >= 1 && nc <= Nx
                if bz_layers(nr, nc) == -1 % If unassigned
                    bz_layers(nr, nc) = current_layer;
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
threshold_map = ones(1, Nx * Ny) + 2; % Initialize with a high value
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
if params.save_figure
    % === Debug plots ===
    figure('Visible', 'off');

    subplot(3,2,1);
    imagesc(flipud(reshape(fibrotic_core', Nx, Ny)')); axis image off; colorbar;
    title(sprintf('Fibrotic Core (%s a=%.2g%%)', params.core_shape, sum(fibrotic_core(:)) * 100 / (Nx*Ny)));

    subplot(3,2,2);
    imagesc(flipud(reshape(bz_layers', Nx, Ny)')); axis image off; colorbar;
    title(sprintf('Border Zone Layers (total=%d)', total_layers));

    subplot(3,2,3);
    imagesc(flipud(reshape(field', Nx, Ny)')); axis image off; colorbar;
    title(sprintf('Noise Field'));

    subplot(3,2,4);
    imagesc(flipud(reshape(threshold_map', Nx, Ny)')); axis image off; colorbar; caxis([threshold_fc max(field)]);
    title(sprintf('Threshold Map (rate=%.1f)', params.threshold_rate));

    subplot(3,2,5);
    imagesc(flipud(reshape(presence_with_categories', Nx, Ny)')); axis image off; colorbar;
    title(sprintf('FC and BZ (FCd=%.4g%% | FCp=%.2g%%)', fc_density*100, fc_proportion*100));

    subplot(3,2,6);
    imagesc(flipud(reshape(presence_with_border_zone', Nx, Ny)')); axis image off; colorbar;
    title(sprintf('Presence'));

    print(params.figure_filename, '-dpng', '-r300');
end

end

% --- Helper Function to Create the Fibrotic Core ---
function fibrotic_core = createFibroticCore(mesh, params)
% This helper function creates the binary mask for the Fibrotic Core (FC).

Nx = mesh.Nx;
Ny = mesh.Ny;
dx = mesh.dx; % Physical size of one element
Lx = Nx * dx; % Total physical width
Ly = Ny * dx; % Total physical height

fibrotic_core = zeros(Ny, Nx);

% --- Determine if using direct dimensions or area fraction ---
use_direct_dims = isfield(params, 'core_width') && isfield(params, 'core_height');

switch lower(params.core_shape)
    case {'rectangle', 'ellipse'}
        w_pixels = 0; h_pixels = 0;

        if use_direct_dims
            % --- METHOD 2: Use direct physical dimensions ---
            if params.core_width > Lx || params.core_height > Ly
                error('Requested core dimensions (%.2f x %.2f) exceed domain limits (%.2f x %.2f).', ...
                      params.core_width, params.core_height, Lx, Ly);
            end
            h_pixels = params.core_width / dx;  % Vertical size (rows) from core_width
            w_pixels = params.core_height / dx; % Horizontal size (cols) from core_height
        else
            % --- METHOD 1: Use area fraction (fallback) ---
            if ~isfield(params, 'core_area_fraction')
                error("For 'rectangle' or 'ellipse', you must provide either core dimensions (core_width, core_height) or core_area_fraction.");
            end
            target_area_pixels = params.core_area_fraction * (Nx * Ny);
            aspect_ratio = Nx / Ny;

            if strcmpi(params.core_shape, 'rectangle')
                % Area = w * h
                h_pixels = sqrt(target_area_pixels / aspect_ratio);
                w_pixels = aspect_ratio * h_pixels;
            else % Ellipse
                % To get the target area, the underlying rectangle must be 4/pi times larger.
                effective_rect_area = target_area_pixels * 4 / pi;
                h_pixels = sqrt(effective_rect_area / aspect_ratio);
                w_pixels = aspect_ratio * h_pixels;
            end
        end

        if strcmpi(params.core_shape, 'rectangle')
            x1 = floor((Nx - w_pixels) / 2) + 1;
            x2 = ceil((Nx + w_pixels) / 2);
            y1 = floor((Ny - h_pixels) / 2) + 1;
            y2 = ceil((Ny + h_pixels) / 2);
            fibrotic_core(y1:y2, x1:x2) = 1;
        else % Ellipse
            [X, Y] = meshgrid(1:Nx, 1:Ny);
            cx = Nx / 2;
            cy = Ny / 2;
            a = w_pixels / 2; % semi-axis a
            b = h_pixels / 2; % semi-axis b
            ellipse_eq = ((X - cx).^2 / a^2) + ((Y - cy).^2 / b^2);
            fibrotic_core(ellipse_eq <= 1) = 1;
        end

    case 'stain'
        if ~isfield(params, 'core_area_fraction')
            error("For 'stain' shape, you must provide 'core_area_fraction'.");
        end
        target_area_fraction = params.core_area_fraction;

        [X_norm, Y_norm] = meshgrid(linspace(-1, 1, Nx), linspace(-1, 1, Ny));
        points = [X_norm(:)'; Y_norm(:)'];
        [perm_table, ~] = generateTables(params.stain_seed);
        stain_noise = Octave2D(points * params.stain_scale, 1, 0.5, perm_table, []);
        stain_noise = reshape(stain_noise, Ny, Nx);

        [X_grad, Y_grad] = meshgrid(linspace(-Nx/2, Nx/2, Nx), linspace(-Ny/2, Ny/2, Ny));
        gradient_radius = sqrt((X_grad/(Nx/2)).^2 + (Y_grad/(Ny/2)).^2);
        containment_gradient = max(0, 1 - gradient_radius.^2);
        stain_noise = stain_noise .* containment_gradient;

        window_thresh = [max(stain_noise(:)), min(stain_noise(:))];
        iters = 0; max_iters = 30; tol = 0.01;

        while iters < max_iters
            check_thresh = mean(window_thresh);
            binary_mask = stain_noise >= check_thresh;
            binary_mask = imfill(binary_mask, 'holes');
            current_area_fraction = sum(binary_mask(:)) / (Nx * Ny);

            err = abs(current_area_fraction - target_area_fraction);
            if err < tol, break; end

            if current_area_fraction > target_area_fraction
                window_thresh(2) = check_thresh;
            else
                window_thresh(1) = check_thresh;
            end
            iters = iters + 1;
        end
        fibrotic_core = stain_noise >= mean(window_thresh);
        fibrotic_core = imfill(fibrotic_core, 'holes');
        current_area_fraction = sum(fibrotic_core(:)) / (Nx * Ny);

    otherwise
        error("Unknown core_shape specified. Use 'rectangle', 'ellipse', or 'stain'.");
end

% Reshape to 1D array to match 'field' format
fibrotic_core = reshape(fibrotic_core, 1, []);
end
