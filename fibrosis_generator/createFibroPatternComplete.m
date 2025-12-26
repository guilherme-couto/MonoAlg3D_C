function [presence, fc_density, bz_layers, fibrotic_core, O_b, O_d, F, noise] = createFibroPatternComplete(mesh, density, params, Ps, offsets, border_zone, variable_direction)
% Generates a binary fibrosis pattern
% with optional border zone and variable fibre direction
%
% INPUTS:
%   mesh               : struct with mesh.points, Nx, Ny
%   density            : desired fibrosis density [0,1]
%   params             : vector with parameters:
%                        [ fibreness, fibre_sep, patchiness, feature_size, roughness,
%                          patch_size, alignment_ratio, fibre_orientation ]
%   Ps                 : Perlin noise permutation tables
%   offsets            : octave noise offsets
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
%   noise          : combined noise field used for pattern generation

% === PARAMETER DECOMPOSITION ===
params = num2cell(params);
[fibreness, fibre_sep, patchiness, feature_size, roughness, patch_size, fibre_alignment, fibre_orientation] = deal(params{:});

% === MESH DATA EXTRACTION ===
points = mesh.points;
Nx = mesh.Nx;
Ny = mesh.Ny;
Ly = max(points(:,2));           % Total height of domain (transmural direction)

if variable_direction

    % === TRANSMURAL DEPTH CALCULATION ===
    y_coords = points(:,2);
    d = y_coords / Ly;               % Normalized transmural position [0,1]

    % === TRANSMURAL ROTATION PARAMETERS ===
    beta_endo = deg2rad(40);     % Endocardium: +80°
    beta_epi  = deg2rad(-50);    % Epicardium: -60°

    % === LOCAL FIBRE ANGLE FIELD ===
    beta = beta_endo * (1 - d) + beta_epi * d;

    % === PROJECT POINTS INTO LOCAL FIBRE BASIS ===
    R_points = projectToRotatedBasis(points, beta);  % 2 x N matrix

else

    % === ROTATED POINTS ===
    R_points = [ [ cos(fibre_orientation), sin(fibre_orientation)];
                 [-sin(fibre_orientation), cos(fibre_orientation)] ] * points';

end

% === PRECOMPUTE PERMUTATION TABLES ===
for k = 1:size(Ps,1)
    Ps2(k,:) = Ps(k, Ps(k,:)+1);
    Ps3(k,:) = Ps(k, Ps2(k,:)+1);
end

% === STEP 1: MAIN FIBROSIS NOISE FIELD (O_b) ===
P_f_points = [ R_points(1,:) / sqrt(fibre_alignment);
               R_points(2,:) * sqrt(fibre_alignment) ];

O_b = Octave2D(P_f_points / feature_size, 4, roughness, Ps, offsets);

% === STEP 2: FIBRE-SELECTING FIELD (STRUCTURED PATTERN) ===
use_structured = ~(isnan(fibreness) || isnan(fibre_sep));
if use_structured

    n_fibres_similarity   = 4;
    wiggle_feature_length = 4;
    phasefield_strength   = 5;

    phasefield_points = [ R_points(1,:) / wiggle_feature_length;
                          R_points(2,:) / (n_fibres_similarity * fibre_sep) ];

    phasefield = Octave2D(phasefield_points, 4, 0.5, Ps2, offsets);

    F = 0.5 + 0.5 * cos(2*pi * (R_points(2,:) / fibre_sep + phasefield_strength * (phasefield - 0.5)));
    F = F.^15;

    % Fibre-modulated noise
    noise = (1 - fibreness + fibreness * F) .* O_b;

else

    % Unstructured pattern: no F
    F = [];
    F = zeros(Ny, Nx);
    noise = O_b;  % no fibre modulation

end

% === STEP 3: LARGE-SCALE DENSITY VARIATION FIELD ===
O_d = Octave2D(points' / patch_size, 4, 0.5, Ps3, offsets);

% === STEP 4: FINAL COMBINATION ===
noise = noise + patchiness * O_d;

% === STEP 5: BINARY THRESHOLDING ===
if border_zone

    % Define the parameters for the new thresholding function in a struct
    threshold_params.core_shape = 'ellipse';
    threshold_params.threshold_rate = 1.0;
    threshold_params.max_layers = 28;
    threshold_params.core_width = 2.8;
    threshold_params.core_height = 2.8;
    % threshold_params.core_area_fraction = 0.53;
    % threshold_params.stain_seed = Ps(1,1);

    % The 'density' variable passed to createFibroPatternComplete is now used as the target density for the Fibrotic Core.
    [presence, fc_density, bz_layers, fibrotic_core] = thresholdPatternWithBorderZone(noise, density, mesh, threshold_params);

else

    presence = thresholdPattern(noise, density);
    fibrotic_core = ones(size(presence)); % All points are considered part of the Fibrotic Core
    fc_density = getPatternDensity(presence, fibrotic_core);  % Actual density achieved in the core
    bz_layers = zeros(Ny, Nx);  % Without Border Zone layers, everything is Fibrotic Core
    bz_layers = reshape(bz_layers, 1, []);  % Reshape to 1D array

end

% === STEP 6: RESHAPE TO MATRIX FORMAT ===
presence  = reshape(presence', Nx, Ny)';
bz_layers = reshape(bz_layers', Nx, Ny)';
noise     = reshape(noise', Nx, Ny)';
O_b       = reshape(O_b', Nx, Ny)';
O_d       = reshape(O_d', Nx, Ny)';
if use_structured
    F = reshape(F', Nx, Ny)';
end

% === STEP 7: FLIP FOR IMAGE DISPLAY ===
presence  = flipud(presence);
bz_layers = flipud(bz_layers);
noise     = flipud(noise);
O_b       = flipud(O_b);
O_d       = flipud(O_d);
if use_structured
    F = flipud(F);
end

end

function local_points = projectToRotatedBasis(points, beta)
%   Projects a set of 2D points into a local rotated basis
%
%   INPUTS:
%       points : N x 2 matrix of 2D coordinates
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

    local_points = [Rx'; Ry'];  % Output as 2 x N
end
