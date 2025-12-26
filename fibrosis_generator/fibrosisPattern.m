function fibrosisPattern(params, density, seed_num, creation_mode, output_filename, save_mesh, save_figure)
%FIBROSISPATTERN creates a fibrosis pattern and saves it to a file.
%
%   fibrosisPattern(params, density, seed_num, creation_mode, output_filename, save_mesh, save_figure)
%
%   Inputs:
%       params          - 1x8 vector, parameters for fibrosis generator
%       density         - scalar, desired fibrosis density
%       seed_num        - scalar, seed for RNG
%       creation_mode   - string, mode of creation ('threshold' or 'composition')
%       output_filename - string, output file name (without extension)
%       save_mesh       - boolean, whether to save the mesh
%       save_figure     - boolean, whether to save the figure
%
%   This function generates a fibrosis pattern domain, and saves it to ALG format.
%
%   The spatial resolution (dx) and dimensions of the fibrotic region are fixed inside the function for consistency.

%% === DOMAIN CONFIGURATION ===
% Fibrotic region size (cm)
Lx_fibro = 4.0;     % Width
Ly_fibro = 4.0;     % Height

% Spatial resolution (cm)
dx = 1/100; % 0.01 cm = 100 µm

%% === GENERATE FIBROSIS PATTERN ===
Nx_fibro = ceil(Lx_fibro / dx);
Ny_fibro = ceil(Ly_fibro / dx);
fibro_mesh = buildMesh(Nx_fibro, Ny_fibro, dx);

border_zone = false;
variable_direction = false;

% Generate binary fibrosis pattern
if strcmpi(creation_mode, 'threshold')
    fibro_pattern_data = generateOnePatternThreshold(params, density, seed_num, fibro_mesh, border_zone, variable_direction);
elseif strcmpi(creation_mode, 'composition')
    fibro_pattern_data = generateOnePatternComposition(params, density, seed_num, 0.005, fibro_mesh, border_zone, variable_direction);
else
    error('Invalid creation mode. Choose either "threshold" or "composition".');
end

%% === SAVE FINAL PATTERN ===
if save_mesh
    writeALGFromPatternData(fibro_pattern_data, output_filename);
end

if save_figure

    % Save pattern as a figure
    figure('Visible', 'off');
    imagesc(fibro_pattern_data.presence);
    fibroclr = [[0.1, 0.5, 0.8]; [0.9, 0.5, 0.1]]; % Blue and orange
    colormap(fibroclr);
    axis equal off;
    title(sprintf('Fibrosis Pattern (Seed: %d) - Fibrosis Density %.3f', seed_num, fibro_pattern_data.fc_density));
    print([output_filename, '.png'], '-dpng', '-r300');

    save_fields = true;

    if save_fields

        figure('Visible', 'off');

        % LAYOUT CONFIGURATION
        margin_top    = 0.08;
        margin_bottom = 0.05;
        margin_left   = 0.05;
        margin_right  = 0.05;

        % --- Heights ---
        top_row_height    = 0.35;  % Height of the top row (3 plots)
        bottom_row_height = 0.45;  % Height of the bottom row (2 plots)

        % --- Spacings ---
        v_spacing = 0.07;   % Vertical space between rows
        h_spacing = 0.001;  % Horizontal space between plots

        % --- Automatic layout calculations ---
        % Widths
        top_plot_width    = (1 - margin_left - margin_right - 2*h_spacing) / 3; % 3 plots, 2 spaces
        bottom_plot_width = (1 - margin_left - margin_right - 1*h_spacing) / 2; % 2 plots, 1 space

        % Vertical Positions (calculated from bottom to top)
        bottom_row_y = margin_bottom;
        top_row_y    = bottom_row_y + bottom_row_height + v_spacing;

        % PLOTS
        % --- Line 1: Components ---
        ax1 = axes('Position', [margin_left, top_row_y, top_plot_width, top_row_height]);
        imagesc(fibro_pattern_data.Ob); axis image off; colormap(ax1, 'default');
        title(['Base Noise (Nb)'], 'FontSize', 10);

        ax2 = axes('Position', [margin_left + top_plot_width + h_spacing, top_row_y, top_plot_width, top_row_height]);
        imagesc(fibro_pattern_data.Od); axis image off; colormap(ax2, 'default');
        title(['Density Variation (Nd)'], 'FontSize', 10);

        ax3 = axes('Position', [margin_left + 2*top_plot_width + 2*h_spacing, top_row_y, top_plot_width, top_row_height]);
        if nnz(fibro_pattern_data.F) > 0
            imagesc(fibro_pattern_data.F, 'Parent', ax3);
            colormap(ax3, 'default');
        end

        axis(ax3, 'image', 'off');
        title(ax3, ['Fiber Selection (F)'], 'FontSize', 10);

        % --- Line 2: Results ---
        ax4 = axes('Position', [margin_left, bottom_row_y, bottom_plot_width, bottom_row_height]);
        imagesc(fibro_pattern_data.noise); axis image off;
        colormap(ax4, 'gray');
        title('Combined Noise Field', 'FontSize', 12);

        ax5 = axes('Position', [margin_left + bottom_plot_width + h_spacing, bottom_row_y, bottom_plot_width, bottom_row_height]);
        imagesc(fibro_pattern_data.presence); axis image off;
        colormap(ax5, fibroclr);
        title('Resulting Pattern', 'FontSize', 12);

        % --- Figure configurations ---
        set(gcf, 'Position', [10, 10, 800, 600]);
        print([output_filename, '_fields.png'], '-dpng', '-r300');

    end
end

end
