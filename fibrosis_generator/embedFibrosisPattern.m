function embedFibrosisPattern(params, density, seed_num, creation_mode, output_filename, save_mesh, save_figure)
%EMBEDFIBROSISPATTERN Embeds a fibrosis pattern into the center of a healthy tissue domain.
%
%   embedFibrosisPattern(params, density, seed_num, creation_mode, output_filename, save_mesh, save_figure)
%
%   Inputs:
%       params          - 1x8 vector, parameters for fibrosis generator
%       density         - scalar, desired fibrosis density in central region
%       seed_num        - scalar, seed for RNG
%       creation_mode   - string, mode of creation ('threshold' or 'composition')
%       output_filename - string, output file name (without extension)
%       save_mesh       - boolean, whether to save the mesh
%       save_figure     - boolean, whether to save the figure
%
%   This function creates a larger healthy tissue mesh, generates a smaller
%   fibrosis pattern, embeds the fibrosis in the center of the healthy
%   domain, and saves the final pattern to ALG format.
%
%   The spatial resolution (dx) and dimensions of both healthy tissue and
%   fibrotic region are fixed inside the function for consistency.

%% === DOMAIN CONFIGURATION ===
% Healthy domain size (cm)
Lx_healthy = 4.0;   % Width
Ly_healthy = 4.0;   % Height

% Fibrotic region size (cm)
Lx_fibro = 3.4;     % Width
Ly_fibro = 3.4;     % Height

% Spatial resolution (cm)
dx = 1/100; % 0.01 cm = 100 µm

%% === CREATE HEALTHY DOMAIN MESH ===
Nx_healthy = ceil(Lx_healthy / dx);
Ny_healthy = ceil(Ly_healthy / dx);

% Start with all healthy tissue (presence: 0 = healthy; bz_layers: -1 = healthy)
presence_healthy = zeros(Ny_healthy, Nx_healthy);
bz_layers_healthy = -ones(Ny_healthy, Nx_healthy);

%% === GENERATE FIBROSIS PATTERN ===
Nx_fibro = ceil(Lx_fibro / dx);
Ny_fibro = ceil(Ly_fibro / dx);
fibro_mesh = buildMesh(Nx_fibro, Ny_fibro, dx);

border_zone = true;
variable_direction = false;

% Generate binary fibrosis pattern
if strcmpi(creation_mode, 'threshold')
    fibro_pattern_data = generateOnePatternThreshold(params, density, seed_num, fibro_mesh, border_zone, variable_direction);
elseif strcmpi(creation_mode, 'composition')
    fibro_pattern_data = generateOnePatternComposition(params, density, seed_num, 0.005, fibro_mesh, border_zone, variable_direction);
else
    error('Invalid creation mode. Choose either "threshold" or "composition".');
end

%% === EMBED FIBROSIS INTO HEALTHY DOMAIN ===
start_x = floor((Nx_healthy - Nx_fibro) / 2) + 1;
start_y = floor((Ny_healthy - Ny_fibro) / 2) + 1;

presence_healthy(start_y:start_y+Ny_fibro-1, start_x:start_x+Nx_fibro-1) = fibro_pattern_data.presence;
bz_layers_healthy(start_y:start_y+Ny_fibro-1, start_x:start_x+Nx_fibro-1) = fibro_pattern_data.bz_layers;

%% === BUILD EMBEDDED PATTERN DATA ===
embedded_pattern_data.presence          = presence_healthy;
embedded_pattern_data.fc_density        = fibro_pattern_data.fc_density;
embedded_pattern_data.bz_layers         = bz_layers_healthy;
embedded_pattern_data.fiber_orientation = fibro_pattern_data.fiber_orientation;
embedded_pattern_data.mesh              = buildMesh(Nx_healthy, Ny_healthy, dx);
embedded_pattern_data.seed_num          = seed_num;

%% === SAVE FINAL PATTERN ===
if save_mesh
    writeALGFromPatternData(embedded_pattern_data, output_filename);
end

if save_figure

    % Save pattern as a figure
    figure('Visible', 'off');
    imagesc(embedded_pattern_data.presence);
    fibroclr = [[0.1, 0.5, 0.8]; [0.9, 0.5, 0.1]]; % Blue and orange
    colormap(fibroclr);
    axis equal off;
    title(sprintf('Embedded Fibrosis Pattern (Seed: %d) - Fibrosis Density %.3f', seed_num, embedded_pattern_data.fc_density));
    print([output_filename, '.png'], '-dpng', '-r300');

    save_debug_figure = true;

    if save_debug_figure

        % Proportion of fibrosis inside the FC over total fibrosis
        fibrosis_count_inside_fc = sum(fibro_pattern_data.fc(:) & fibro_pattern_data.presence(:));
        fibrosis_count_outside_fc = sum(~fibro_pattern_data.fc(:) & fibro_pattern_data.presence(:));
        fc_proportion = fibrosis_count_inside_fc / (fibrosis_count_inside_fc + fibrosis_count_outside_fc);
        fibro_fc = reshape(fibro_pattern_data.fc', Nx_fibro, Ny_fibro)';
        presence_with_categories = fibro_pattern_data.presence + (fibro_fc & fibro_pattern_data.presence); % 0 = Healthy, 1 = FC Fibrosis, 2 = BZ Fibrosis

        figure('Visible', 'off');

        % --- Figure Margins ---
        margin_top    = 0.08;
        margin_bottom = 0.05;
        margin_left   = 0.05;
        margin_right  = 0.05;

        % --- Proportions ---
        top_row_height    = 0.25;  % Height of the top row
        bottom_row_height = 0.55;  % Height of the bottom row
        v_spacing         = 0.07;  % Vertical space between rows
        h_spacing         = 0.04;  % Horizontal space between the 3 plots on top

        % --- Automatic layout calculations ---
        bottom_row_y = margin_bottom;
        top_row_y = margin_bottom + bottom_row_height + v_spacing;
        top_plot_width = (1 - margin_left - margin_right - 2*h_spacing) / 3;
        bottom_plot_width = 1 - margin_left - margin_right;

        % --- Top Row (smaller plots) ---
        axes('Position', [margin_left, top_row_y, top_plot_width, top_row_height]);
        imagesc(fibro_fc); axis image off; colormap(gca, 'default');
        title(['Fibrotic Core', sprintf('\n(a=%.2g%%)', sum(fibro_pattern_data.fc(:)) * 100 / (Nx_fibro*Ny_fibro))], 'FontSize', 8);

        axes('Position', [margin_left + top_plot_width + h_spacing, top_row_y, top_plot_width, top_row_height]);
        imagesc(flipud(fibro_pattern_data.bz_layers)); axis image off; colormap(gca, 'default');
        title(['Border Zone Layers', sprintf('\n(total=%d)', max(fibro_pattern_data.bz_layers(:)))], 'FontSize', 8);

        axes('Position', [margin_left + 2*top_plot_width + 2*h_spacing, top_row_y, top_plot_width, top_row_height]);
        imagesc(reshape(presence_with_categories', Nx_fibro, Ny_fibro)'); axis image off; colormap(gca, 'default');
        title(['FC and BZ', sprintf('\n(FCd=%.4g%% | FCp=%.2g%%)', fibro_pattern_data.fc_density*100, fc_proportion*100)], 'FontSize', 8);

        % --- Bottom Row (Main Plot) ---
        ax4 = axes('Position', [margin_left, bottom_row_y, bottom_plot_width, bottom_row_height]);
        imagesc(presence_healthy); axis image off;
        colormap(ax4, fibroclr);
        title('Embedded Fibrosis', 'FontSize', 10);

        % --- Figure configurations ---
        set(gcf, 'Position', [10, 10, 400, 600]);
        print([output_filename, '_debug.png'], '-dpng', '-r300');

    end

end

end
