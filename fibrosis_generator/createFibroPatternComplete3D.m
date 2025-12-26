function [presence, fc_density, bz_layers, fibrotic_core, O_b, O_d, F] = createFibroPatternComplete3D(mesh, density, params, Ps, offsets1, offsets2, border_zone, variable_direction)
% Generates a binary fibrosis pattern (3D)
% with optional border zone and variable fibre direction
%
% INPUTS:
%   mesh               : struct with mesh.points, Nx, Ny, Nz
%   density            : desired fibrosis density [0,1]
%   params             : vector with parameters:
%                        [ fibreness, fibre_sep1, fibre_sep2, patchiness, feature_size, roughness,
%                          patch_size, fibre_alignment_y, fibre_alignment_z, phi, theta ]
%   Ps                 : Perlin noise permutation tables
%   offsets1            : octave noise offsets
%   offsets2            : octave noise offsets
%   border_zone        : (bool) whether to apply border zone
%   variable_direction : (bool) whether to apply variable fibre direction
%
% OUTPUTS:
%   presence       : binary pattern (1 = fibrotic)
%   fc_density     : actual fibrosis density within the core
%   bz_layers      : border zone layer map
%   fibrotic_core  : binary mask indicating the FC (1) and BZ (0)
%   O_b            : base noise field
%   O_d            : large-scale density noise field
%   F              : fibre-aligned sinusoidal field (empty if unstructured)

% === PARAMETER DECOMPOSITION ===
params = num2cell(params);
[fibreness, fibre_sep1, fibre_sep2, patchiness, feature_size, roughness, patch_size, fibre_alignment_y, fibre_alignment_z, phi, theta] = deal(params{:});

% Non-fibre cases are properly set to have zero fibreness
if isnan(fibreness)
    fibreness  = 0;
    fibre_sep1 = 1;
    fibre_sep2 = 1;
end

% === MESH DATA EXTRACTION ===
points = mesh.points;
Nx = mesh.Nx;
Ny = mesh.Ny;
Nz = mesh.Nz;

if variable_direction

    % TODO
    a = 0;

else
    % === ROTATED POINTS ===
    % Create rotated points using the 3D rotation matrix
    R = [ cos(phi) * cos(theta), sin(phi), -cos(phi) * sin(theta);
         -cos(theta) * sin(phi), cos(phi),  sin(phi) * sin(theta);
          sin(theta)           , 0       ,  cos(theta) ];

    R_points = R * points'; % 3 x N matrix
end

% === PRECOMPUTE PERMUTATION TABLES ===
for k = 1:size(Ps,1)
    Ps2(k,:) = Ps(k, Ps(k,:)+1);
    Ps3(k,:) = Ps(k, Ps2(k,:)+1);
    Ps4(k,:) = Ps(k, Ps3(k,:)+1);
end

% === STEP 1: MAIN FIBROSIS NOISE FIELD (O_b) ===
P_f_points = [ R_points(1,:) / fibre_alignment_y^(1/3) / fibre_alignment_z^(1/3);
               R_points(2,:) * fibre_alignment_y^(2/3) / fibre_alignment_z^(1/3);
               R_points(3,:) / fibre_alignment_y^(1/3) * fibre_alignment_z^(2/3) ];

O_b = Octave3D(P_f_points / feature_size, 4, roughness, Ps, offsets1);

% === STEP 2: FIBRE-SELECTING FIELD (STRUCTURED PATTERN) ===
n_fibres_similarity   = 4;
wiggle_feature_length = 4;
phasefield_strength   = 5;

phasefield_points = [ R_points(1,:) / wiggle_feature_length;
                      R_points(2,:) / (n_fibres_similarity * fibre_sep1);
                      R_points(3,:) / (n_fibres_similarity * fibre_sep2) ];

phasefield1 = Octave3D(phasefield_points, 4, 0.5, Ps2, offsets1);
phasefield2 = Octave3D(phasefield_points, 4, 0.5, Ps3, offsets2);

F = 0.5 + 0.5 * cos(2 * pi * (R_points(2,:) / fibre_sep1 + phasefield_strength * (phasefield1 - 0.5))) .* cos(2 * pi * (R_points(3,:) / fibre_sep2 + phasefield_strength * (phasefield2 - 0.5)));
F = F.^15;

% Fibre-modulated noise
noise = (1 - fibreness + fibreness * F) .* O_b;

% === STEP 3: LARGE-SCALE DENSITY VARIATION FIELD ===
% Use base points (non-rotated) stretched by patch size for density variation field
O_d = Octave3D(points' / patch_size, 4, 0.5, Ps4, offsets1);

% === STEP 4: FINAL COMBINATION ===
noise = noise + patchiness * O_d;

% === STEP 5: BINARY THRESHOLDING ===
if border_zone

    % Define the parameters for the new thresholding function in a struct
    threshold_params.core_shape = 'ellipsoid';
    threshold_params.threshold_rate = 1.0;
    threshold_params.max_layers = 20;
    threshold_params.core_width = 2.0;
    threshold_params.core_height = 2.0;
    threshold_params.core_depth = 0.6;

    % The 'density' variable passed to createFibroPatternComplete is now used as the target density for the Fibrotic Core.
    [presence, fc_density, bz_layers, fibrotic_core] = thresholdPatternWithBorderZone3D(noise, density, mesh, threshold_params);

else

    presence = thresholdPattern(noise, density);
    fibrotic_core = ones(size(presence)); % All points are considered part of the Fibrotic Core
    fc_density = getPatternDensity(presence, fibrotic_core);  % Actual density achieved in the core
    bz_layers = zeros(Ny, Nx, Nz);  % Without Border Zone layers, everything is Fibrotic Core
    bz_layers = reshape(bz_layers, 1, []);  % Reshape to 1D array

end

% === STEP 6: RESHAPE TO MATRIX FORMAT ===
presence  = reshape(presence', Nx, Ny, Nz);
bz_layers = reshape(bz_layers', Nx, Ny, Nz);
O_b       = reshape(O_b', Nx, Ny, Nz);
O_d       = reshape(O_d', Nx, Ny, Nz);

if ~isnan(fibreness)

    F = reshape(F', Nx, Ny, Nz);

end

end


function local_points = projectToRotatedBasis(points, beta)
%   Projects a set of 3D points into a local rotated basis
%
%   INPUTS:
%       points : N x 3 matrix of 3D coordinates
%       beta   : N x 1 vector of rotation angles (in radians) for each point
%
%   OUTPUT:
%       local_points : 2 x N matrix with points projected into local bases
%
%   Each point is projected onto its local fibre-aligned coordinate system
%   defined by the unit vectors:
%       e1 = [cos(beta); sin(beta)]  (local x, aligned with fibre)
%       e2 = [-sin(beta); cos(beta)] (local y, transverse to fibre)

    x = points(:,1);
    y = points(:,2);
    z = points(:,3);

    % Local basis vectors
    e1x = cos(beta);
    e1y = sin(beta);
    e2x = -sin(beta);
    e2y = cos(beta);

    % Projection onto local basis
    % Rx​ = ⟨(xi​,yi​),e1​⟩ = xi​*cos(βi​) + yi​*sin(βi​)
    % Ry = ⟨(xi,yi),e2⟩ = −xi*sin⁡(βi) + yi*cos⁡(βi)
    Rx = e1x .* x + e1y .* y;  % projection onto e1 (fibre direction)
    Ry = e2x .* x + e2y .* y;  % projection onto e2 (transverse direction)

    local_points = [Rx'; Ry'; z'];  % Output as 3 x N
end
