function writeALGFromPatternData3D(pattern_data, filename)
% Saves the fibrosis pattern to a .alg file compatible with the finite volume solver.
% The file uses scientific notation and comma as separator.
%
% Inputs:
%   pattern_data - pattern structure with fields:
%      - presence            : the binary presence pattern (3D matrix of 0s and 1s)
%      - fc_density          : the density of fibrosis in the core region
%      - bz_layers           : border zone layer map (-1 = Healthy, 0 = Fibrotic Core, >0 = Border Zone Layer)
%      - fiber_orientations  : [phi, theta], the orientation of fibers (in radians)
%      - mesh                : the mesh structure used for the pattern
%      - seed_num            : the seed number used for random generation
%

    % Get the dimensions directly from the new mesh struct fields
    Nx = pattern_data.mesh.Nx;
    Ny = pattern_data.mesh.Ny;
    Nz = pattern_data.mesh.Nz;

    presence = reshape(pattern_data.presence, Ny, Nx, Nz);
    bz_layers = reshape(pattern_data.bz_layers, Ny, Nx, Nz);

    dx_um = pattern_data.mesh.dx * 1e4;         % Convert voxel size from cm to µm

    % Fiber direction (unit vector)
    % phi -> rotation around z-axis, theta -> rotation around the new y-axis
    phi = pattern_data.fiber_orientations(1);
    theta = pattern_data.fiber_orientations(2);
    fiber_x = cos(phi) * cos(theta);
    fiber_y = -sin(phi) * cos(theta);
    fiber_z = sin(theta);

    % Get number of layers to normalize bz_multiplier
    num_bz_layers = max(bz_layers(:)) + 1;  % +1 to include healthy layer

    % Open file for writing
    fid = fopen([filename, '.alg'], 'wt');

    for j = 1:Ny         % y-direction (rows)
        for i = 1:Nx     % x-direction (columns)
            for k = 1:Nz % z-direction (pages)
                x_center  = (i - 0.5) * dx_um;
                % y_center  = ((Ny - j + 0.5)) * dx_um;
                y_center  = (j - 0.5) * dx_um;
                z_center  = (k - 0.5) * dx_um;

                fibro_tag = presence(i, j, k);  % row = j, col = i, page = k
                bz_layer  = bz_layers(i, j, k);  % row = j, col = i, page = k

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
    end

    fclose(fid);
end
