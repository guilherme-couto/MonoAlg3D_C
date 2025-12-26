function embedFibrosisPattern3D(params, density, seed_num, creation_mode, output_filename, save_mesh, save_figure)
%EMBEDFIBROSISPATTERN3D Embeds a 3D fibrosis pattern into the center of a healthy tissue domain.
%
%   embedFibrosisPattern3D(params, density, seed_num, creation_mode, output_filename, save_mesh, save_figure)
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
Lx_healthy = 1.0;   % Width
Ly_healthy = 1.0;   % Height
Lz_healthy = 1.0;   % Depth

% Fibrotic region size (cm)
Lx_fibro = 0.85;     % Width
Ly_fibro = 0.85;     % Height
Lz_fibro = 0.85;     % Depth

% Spatial resolution (cm)
dx = 1/100; % 0.01 cm = 100 µm

%% === CREATE HEALTHY DOMAIN MESH ===
Nx_healthy = ceil(Lx_healthy / dx);
Ny_healthy = ceil(Ly_healthy / dx);
Nz_healthy = ceil(Lz_healthy / dx);

% Start with all healthy tissue (presence: 0 = healthy; bz_layers: -1 = healthy)
presence_healthy = zeros(Ny_healthy, Nx_healthy, Nz_healthy);
bz_layers_healthy = -ones(Ny_healthy, Nx_healthy, Nz_healthy);

%% === GENERATE FIBROSIS PATTERN ===
Nx_fibro = ceil(Lx_fibro / dx);
Ny_fibro = ceil(Ly_fibro / dx);
Nz_fibro = ceil(Lz_fibro / dx);
fibro_mesh = buildMesh3D(Nx_fibro, Ny_fibro, Nz_fibro, dx);

border_zone = true;
variable_direction = false;

% Generate binary fibrosis pattern
if strcmpi(creation_mode, 'threshold')
    fibro_pattern_data3D = generateOnePatternThreshold3D(params, density, seed_num, fibro_mesh, border_zone, variable_direction);
elseif strcmpi(creation_mode, 'composition')
    fibro_pattern_data3D = generateOnePatternComposition3D(params, density, seed_num, 0.005, fibro_mesh, border_zone, variable_direction);
else
    error('Invalid creation mode. Choose either "threshold" or "composition".');
end

%% === EMBED FIBROSIS INTO HEALTHY DOMAIN ===
start_x = floor((Nx_healthy - Nx_fibro) / 2) + 1;
start_y = floor((Ny_healthy - Ny_fibro) / 2) + 1;
start_z = floor((Nz_healthy - Nz_fibro) / 2) + 1;

presence_healthy(start_y:start_y+Ny_fibro-1, start_x:start_x+Nx_fibro-1, start_z:start_z+Nz_fibro-1) = fibro_pattern_data3D.presence;
bz_layers_healthy(start_y:start_y+Ny_fibro-1, start_x:start_x+Nx_fibro-1, start_z:start_z+Nz_fibro-1) = fibro_pattern_data3D.bz_layers;

%% === BUILD EMBEDDED PATTERN DATA ===
embedded_pattern_data3D.presence           = presence_healthy;
embedded_pattern_data3D.fc_density         = fibro_pattern_data3D.fc_density;
embedded_pattern_data3D.bz_layers          = bz_layers_healthy;
embedded_pattern_data3D.fiber_orientations = fibro_pattern_data3D.fiber_orientations;
embedded_pattern_data3D.mesh               = buildMesh3D(Nx_healthy, Ny_healthy, Nz_healthy, dx);
embedded_pattern_data3D.seed_num           = seed_num;

%% === SAVE FINAL PATTERN ===
if save_mesh
    writeALGFromPatternData3D(embedded_pattern_data3D, output_filename);
end

end
