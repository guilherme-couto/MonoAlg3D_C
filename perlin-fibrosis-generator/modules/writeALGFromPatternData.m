function writeALGFromPatternData(pattern_data, filename)
% Saves the fibrosis pattern to a .alg file.
% Assumes Input Data is (Ny, Nx, [Nz]) where Row 1 is Y=0 (Bottom).

    Nx = pattern_data.mesh.Nx;
    Ny = pattern_data.mesh.Ny;
    Nz = pattern_data.mesh.Nz;
    is3D = pattern_data.mesh.is3D;

    dx_um = pattern_data.mesh.dx * 1e4; % Convert cm to µm

    % Fiber direction setup
    % When fiber varies per point, this will need to be updated
    fiber_x = cos(pattern_data.params.orientation);
    fiber_y = sin(pattern_data.params.orientation);
    fiber_z = 0;
    if is3D
        fiber_x = cos(pattern_data.params.phi) * cos(pattern_data.params.theta);
        fiber_y = -sin(pattern_data.params.phi) * cos(pattern_data.params.theta);
        fiber_z = sin(pattern_data.params.theta);
    end

    num_bz_layers = max(pattern_data.bz_layers(:)) + 1; % +1 to include healthy layer

    fid = fopen([filename, '.alg'], 'wt');
    if fid == -1, error('Could not create output file'); end

    % === LOOP STRUCTURE ===
    % Standard Order: Z -> Y -> X
    for k = 1:Nz
        % Z Coordinate (Bottom-Up)
        z_center = (k - 0.5) * dx_um;

        for j = 1:Ny
            % Y Coordinate (Bottom-Up)
            % Row 1 (j=1) is Y=0. Row Ny is Y=Max.
            y_center = (j - 0.5) * dx_um;

            for i = 1:Nx
                % X Coordinate (Left-Right)
                x_center = (i - 0.5) * dx_um;

                % Extract Data directly (Ny, Nx, Nz) implies (Row, Col, Slice)
                % Access (j, i, k)
                if is3D
                    fibro_tag = pattern_data.presence(j, i, k);
                    bz_layer  = pattern_data.bz_layers(j, i, k);
                else
                    fibro_tag = pattern_data.presence(j, i);
                    bz_layer  = pattern_data.bz_layers(j, i);
                end

                % BZ Logic
                if bz_layer == -1
                    bz_multiplier = 1; % healthy layer
                else
                    bz_multiplier = bz_layer / num_bz_layers;
                end

                % Write
                fprintf(fid, '%d,%d,%d,%d,%d,%d,%d,%.5g,%.5g,%.5g,%d,%.5g\n', ...
                    x_center, y_center, z_center, ...
                    dx_um*0.5, dx_um*0.5, dx_um*0.5, ...
                    fibro_tag, fiber_x, fiber_y, fiber_z, ...
                    bz_layer, bz_multiplier);
            end
        end
    end
    fclose(fid);
end
