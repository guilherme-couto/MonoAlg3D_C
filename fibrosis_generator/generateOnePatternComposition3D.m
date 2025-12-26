function pattern_data3D = generateOnePatternComposition3D(params, density, seed_num, tolerance, mesh, border_zone, variable_direction)
% This function takes the provided set of parameter values, and creates the
% representation of a fibrosis pattern with the desired density.
%
% INPUTS:
%
% params - the parameter values for the generator to use (1x9 vector)
% density - density of fibrosis in the pattern to be generated
% seed_num - seed for the random number generator
% tolerance - tolerance for density matching
% mesh - mesh to specify size of pattern
% border_zone - (bool) whether to include a border zone in the pattern
% variable_direction - (bool) whether to allow variable fiber direction
%
% PARAMETER INFORMATION:
%
% 1 - FIBRENESS: The extent to which patterns exhibit long, thin fibres
%     aligned in consistent directions
%     ::: If set to NaN, a pattern without fibres will be created :::
% 2 - FIBRE SEPARATION: The average spacing between fibres (in units
%     matching those used in input mesh
% 3 - PATCHINESS: The extent of inconsistency in pattern density (high
%     patchiness will produce distinct regions of higher and lower density)
% 4 - FEATURE SIZE: The overall size of obstacle features in the mesh (in
%     units matching the mesh)
% 5 - ROUGHNESS: The roughness of feature edges (values from [0,1], may
%     cause issues for values of precisely 0 or 1)
% 6 - PATCH SIZE: The size of regions over which density varies (with
%     extent according to PATCHINESS)
% 7 - FIBRE ALIGNMENT: The extent to which non-fibre features are aligned
%     to the fibre direction (i.e. extent of feature anisotropy)
% 8 - PHI: An angle (in radians) defining the rotation of fibres
%     around the z-axis
% 9 - THETA: An angle (in radians) defining the rotation of fibres
%     around the y-axis
%
% OUTPUTS:
%
% pattern_data3D - the generated pattern structure with fields:
%   - presence            : the binary presence pattern (3D matrix of 0s and 1s)
%   - fc                  : the fibrotic core region (3D matrix of 0s and 1s)
%   - fc_density          : the density of fibrosis in the core region
%   - bz_layers           : border zone layer map (-1 = Healthy, 0 = Fibrotic Core, >0 = Border Zone Layer)
%   - fiber_orientations  : [phi, theta], the orientation of fibers (in radians)
%   - mesh                : the mesh structure used for the pattern
%   - seed_num            : the seed number used for random generation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the modifier to anisotropy in the third dimension
K = 1;  % no modification

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate permutation and offset tables from seed
[permute_table, offset_table1, offset_table2] = generateTables3D(seed_num);

% Convert parameters as needed in moving from 2D to 3D
feature_size3 = params(4) .* params(7).^(-1/6) * K^(-1/3);
anisotropy_y = params(7);
anisotropy_z = params(7) * K;

% Set angles (phi = around z-axis; theta = around y-axis)
phi = params(8);
theta = params(9);

% Convert from the 2D learned parameters to a 3D equivalent
params3D = [params(1), params(2), params(2), params(3), feature_size3, params(5), params(6), anisotropy_y, anisotropy_z, phi, theta];

% Calculate the threshold for the current density
threshold = getCompositionThreshold(density, tolerance);

% Use the 3D version of the generator to get a binary pattern
[presence, fc_density, bz_layers, fibrotic_core] = createFibroPatternComplete3D(mesh, threshold, params3D, permute_table, offset_table1, offset_table2, border_zone, variable_direction);

% TEMPORARY: Create a fibrotic_core the same size as the presence pattern to measure density
% fibrotic_core = ones(size(presence));
% fc_density = getPatternDensity(presence, fibrotic_core);

% Get the density of the pattern generated
actual_density = fc_density;

% Initialise counters
iterations = 0;
tries_counter = 0;
max_tries = 10;

% While the density is not the desired one, keep generating patterns and summing them
while abs(actual_density - density) > tolerance

    iterations = iterations + 1;

    % Generate a new pattern
    [permute_table, offset_table1, offset_table2] = generateTables3D(seed_num + iterations);
    abs_diff = abs(actual_density - density);
    threshold = getCompositionThreshold(abs_diff, tolerance);
    aux_presence = createFibroPatternComplete3D(mesh, threshold, params3D, permute_table, offset_table1, offset_table2, border_zone, variable_direction);

    % Sum the patterns
    aux_presence = aux_presence | presence;

    % Get the density of the sum of the patterns
    aux_density = getPatternDensity(aux_presence, fibrotic_core);

    % Update the tries counter
    tries_counter = tries_counter + 1;

    % Only update the presence if the density obtained is less than the desired density
    if aux_density < density + tolerance
        presence = aux_presence;
        actual_density = aux_density;
        tries_counter = 0;
    end

    % If the number of tries is greater than the maximum number of tries, restart the presence
    if tries_counter > max_tries
        % Generate a new pattern
        threshold = getCompositionThreshold(density, tolerance);
        [presence, fc_density] = createFibroPatternComplete3D(mesh, threshold, params3D, permute_table, offset_table1, offset_table2, border_zone, variable_direction);
        % actual_density = getPatternDensity(presence, fibrotic_core);
        actual_density = fc_density;
        tries_counter = 0;
    end
end

% Store all pattern data in a single struct
pattern_data3D.presence            = presence;
pattern_data3D.fc                  = fibrotic_core;
pattern_data3D.fc_density          = actual_density;
pattern_data3D.bz_layers           = bz_layers;
pattern_data3D.fiber_orientations  = [phi, theta];
pattern_data3D.mesh                = mesh;
pattern_data3D.seed_num            = seed_num;

end
