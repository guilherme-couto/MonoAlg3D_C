function pattern_data = generateOnePatternThreshold(params, density, seed_num, mesh, border_zone, variable_direction)
% This function takes the provided set of parameter values, and creates the
% representation of a fibrosis pattern with the desired density.
%
% INPUTS:
%
% params - the parameter values for the generator to use (1x8 vector)
% density - density of fibrosis in the pattern to be generated
% seed_num - seed for the random number generator
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
% 8 - DIRECTION: An angle (in radians) defining the orientation of fibres
%     and/or feature anisotropy
%
% OUTPUTS:
%
% pattern_data - the generated pattern structure with fields:
%   - presence           : the binary presence pattern (2D matrix of 0s and 1s)
%   - fc                 : the fibrotic core region (2D matrix of 0s and 1s)
%   - fc_density         : the density of fibrosis in the core region
%   - bz_layers          : border zone layer map (-1 = Healthy, 0 = Fibrotic Core, >0 = Border Zone Layer)
%   - fiber_orientation  : the orientation of fibers (in radians)
%   - mesh               : the mesh structure used for the pattern
%   - seed_num           : the seed number used for random generation

% mesh = buildMesh(250, 400, 1/136); % Original values from paper

% Generate permutation and offset tables from seed
[permute_table, offset_table] = generateTables(seed_num);

% Generate the fibrosis pattern
[presence, fc_density, bz_layers, fibrotic_core, Ob, Od, F, noise] = createFibroPatternComplete(mesh, density, params, permute_table, offset_table, border_zone, variable_direction);

% Store all pattern data in a single struct
pattern_data.presence          = presence;
pattern_data.fc                = fibrotic_core;
pattern_data.fc_density        = fc_density;
pattern_data.bz_layers         = bz_layers;
pattern_data.fiber_orientation = params(8);
pattern_data.mesh              = mesh;
pattern_data.seed_num          = seed_num;

pattern_data.Ob = Ob;
pattern_data.Od = Od;
pattern_data.F = F;
pattern_data.noise = noise;

end
