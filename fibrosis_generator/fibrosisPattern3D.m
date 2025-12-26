function fibrosisPattern3D(params, density, seed_num, creation_mode, output_filename, save_mesh, save_figure)
%FIBROSISPATTERN3D creates a 3D fibrosis pattern and saves it to a file.
%
%   fibrosisPattern3D(params, density, seed_num, creation_mode, output_filename, save_mesh, save_figure)
%
%   Inputs:
%       params          - 1x9 vector, parameters for fibrosis generator
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
Lx_fibro = 2.0;     % Width
Ly_fibro = 2.0;     % Height
Lz_fibro = 0.6;     % Depth

% Spatial resolution (cm)
dx = 1/100; % 0.01 cm = 100 µm

%% === GENERATE FIBROSIS PATTERN ===
Nx_fibro = ceil(Lx_fibro / dx);
Ny_fibro = ceil(Ly_fibro / dx);
Nz_fibro = ceil(Lz_fibro / dx);

fibro_mesh = buildMesh3D(Nx_fibro, Ny_fibro, Nz_fibro, dx);

border_zone = false;
variable_direction = false;

% Generate binary fibrosis pattern
if strcmpi(creation_mode, 'threshold')
    fibro_pattern_data3D = generateOnePatternThreshold3D(params, density, seed_num, fibro_mesh, border_zone, variable_direction);
elseif strcmpi(creation_mode, 'composition')
    fibro_pattern_data3D = generateOnePatternComposition3D(params, density, seed_num, 0.005, fibro_mesh, border_zone, variable_direction);
else
    error('Invalid creation mode. Choose either "threshold" or "composition".');
end

%% === SAVE FINAL PATTERN ===
% writeALGFromPatternData(fibro_pattern_data, output_filename);
mesh3D = fibro_pattern_data3D.mesh;
presence3D = fibro_pattern_data3D.presence;
fc_density = fibro_pattern_data3D.fc_density;

%% === SAVE FINAL PATTERN ===
if save_mesh
    writeALGFromPatternData3D(fibro_pattern_data3D, output_filename);
end

% disp(['Pattern density: ', num2str(fibro_pattern_data3D.fc_density)]);

end
