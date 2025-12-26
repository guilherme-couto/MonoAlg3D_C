function writeALGFromPatternData(pattern_data, filename)
% Saves the fibrosis pattern to a .alg file compatible with the finite volume solver.
% The file uses scientific notation and comma as separator.
%
% Inputs:
%   pattern_data - pattern structure with fields:
%      - presence           : the binary presence pattern (2D matrix of 0s and 1s)
%      - fc_density         : the density of fibrosis in the core region
%      - bz_layers          : border zone layer map (-1 = Healthy, 0 = Fibrotic Core, >0 = Border Zone Layer)
%      - fiber_orientation  : the orientation of fibers (in radians)
%      - mesh               : the mesh structure used for the pattern
%      - seed_num           : the seed number used for random generation
%
    [Ny, Nx] = size(pattern_data.presence);   % Matrix: Ny rows (y), Nx columns (x)
    dx_um = pattern_data.mesh.dx * 1e4;         % Convert voxel size from cm to µm

    % Fiber direction (unit vector)
    fiber_x = cos(pattern_data.fiber_orientation);
    fiber_y = sin(pattern_data.fiber_orientation);
    fiber_z = 0;

    % Z-center is fixed for all volumes (2D slice)
    z_center = dx_um / 2;

    % Get number of layers to normalize bz_multiplier
    num_bz_layers = max(pattern_data.bz_layers(:)) + 1;  % +1 to include healthy layer

    % Open file for writing
    fid = fopen([filename, '.alg'], 'wt');

    for j = 1:Ny     % y-direction (rows)
        for i = 1:Nx % x-direction (columns)
            x_center  = (i - 0.5) * dx_um;
            y_center  = ((Ny - j + 0.5)) * dx_um;
            fibro_tag = pattern_data.presence(j, i);  % row = j, col = i
            bz_layer  = pattern_data.bz_layers(j, i);  % row = j, col = i
            if bz_layer == -1
                bz_multiplier = 1;  % Healthy layer
            else
                bz_multiplier = bz_layer / num_bz_layers;
            end

            % Write the data in the required format
            % x_center, y_center, z_center, dx/2, dy/2, dz/2, fibro_tag, fiber_x, fiber_y, fiber_z, bz_layer, bz_multiplier
            fprintf(fid, '%d,%d,%d,%d,%d,%d,%d,%.5g,%.5g,%.5g,%d,%.5g\n', ...
                x_center, y_center, z_center, ...
                dx_um*0.5, dx_um*0.5, dx_um*0.5, ...
                fibro_tag, ...
                fiber_x, fiber_y, fiber_z, ...
                bz_layer, bz_multiplier);
        end
    end

    fclose(fid);
end
