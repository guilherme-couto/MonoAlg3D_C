function combined_noise = computeNoiseField(mesh, params, seed)
% Generates a normalized floating point noise field [0, 1].
% Does NOT know about fibrosis, just texture.
% Supports 2D and 3D generation based on mesh dimensions.

    % Check Dimensionality
    is3D = mesh.is3D;

    % === COMMON PARAMETER EXTRACTION ===
    fibreness    = params.fibreness;
    patchiness   = params.patchiness;
    feature_size = params.feature_size;
    roughness    = params.roughness;
    patch_size   = params.patch_size;

    variable_direction = params.variable_direction;

    % Generate Tables
    if nargin < 3, seed = params.seed; end
    [perm_table, offsets1, offsets2] = generateTables(seed, is3D);

    % Precompute perm tables for structured noise
    perm_table2 = zeros(size(perm_table), 'int32');
    perm_table3 = zeros(size(perm_table), 'int32');
    perm_table4 = zeros(size(perm_table), 'int32');

    for k = 1:size(perm_table,1)
        perm_table2(k,:) = perm_table(k, perm_table(k,:)+1);
        perm_table3(k,:) = perm_table(k, perm_table2(k,:)+1);
        if is3D
            perm_table4(k,:) = perm_table(k, perm_table3(k,:)+1);
        end
    end

    % =====================================================================
    % BRANCH: 3D GENERATION
    % =====================================================================
    if is3D
        % Extract 3D Specific Params
        fibre_sep_y = params.fibre_sep_y;
        fibre_sep_z = params.fibre_sep_z;
        align_y     = params.alignment_y;
        align_z     = params.alignment_z;
        phi         = params.phi;
        theta       = params.theta;

        % Handle Non-Fibrous Case
        if isnan(fibreness)
            fibreness = 0; fibre_sep_y = 1; fibre_sep_z = 1;
        end

        % Coordinate Rotation (Global Direction - No variable direction yet)
        points = mesh.points; % Nx3 Matrix

        % Rotation Matrix
        R = [ cos(phi) * cos(theta), sin(phi), -cos(phi) * sin(theta);
             -cos(theta) * sin(phi), cos(phi),  sin(phi) * sin(theta);
              sin(theta)           , 0       ,  cos(theta) ];

        P_rot = R * points'; % 3xN

        % STEP 1: MAIN NOISE FIELD (O_b)
        % Anisotropic scaling
        P_f = [ P_rot(1,:) / align_y^(1/3) / align_z^(1/3);
                P_rot(2,:) * align_y^(2/3) / align_z^(1/3);
                P_rot(3,:) / align_y^(1/3) * align_z^(2/3) ];

        O_b = Octave3D(P_f / feature_size, 4, roughness, perm_table, offsets1);

        % STEP 2: FIBRE-SELECTING FIELD (STRUCTURED PATTERN)
        n_fibres_similarity = 4; wiggle_feature_length = 4; phasefield_strength = 5;

        phase_pts = [ P_rot(1,:) / wiggle_feature_length;
                    P_rot(2,:) / (n_fibres_similarity * fibre_sep_y);
                    P_rot(3,:) / (n_fibres_similarity * fibre_sep_z) ];

        phasefield1 = Octave3D(phase_pts, 4, 0.5, perm_table2, offsets1);
        phasefield2 = Octave3D(phase_pts, 4, 0.5, perm_table3, offsets2); % Uses second offset table

        term_y = cos(2*pi * (P_rot(2,:) / fibre_sep_y + phasefield_strength * (phasefield1 - 0.5)));
        term_z = cos(2*pi * (P_rot(3,:) / fibre_sep_z + phasefield_strength * (phasefield2 - 0.5)));

        F = 0.5 + 0.5 * term_y .* term_z;
        F = F.^15;

        combined_noise = (1 - fibreness + fibreness * F) .* O_b;

        % STEP 3: DENSITY VARIATION FIELD (O_d)
        % Use original points (non-rotated) stretched by patch size for density variation field
        O_d = Octave3D(points' / patch_size, 4, 0.5, perm_table4, offsets1);

        % Final Combination
        combined_noise = combined_noise + patchiness * O_d;

    % =====================================================================
    % BRANCH: 2D GENERATION
    % =====================================================================
    else
        % Extract 2D Params
        fibre_sep   = params.fibre_sep;
        align       = params.alignment;
        orientation = params.orientation;
        var_dir     = params.variable_direction;

        % Handle Non-Fibrous Case
        if isnan(fibreness)
            fibreness = 0; fibre_sep = 1;
        end

        points = mesh.points; % Nx2
        Ly = mesh.Ny * mesh.dx;

        % Rotation Logic
        if var_dir
            y_coords = points(:,2);
            depth = y_coords / Ly; % Normalized transmural position [0,1]
            beta_endo = deg2rad(40);
            beta_epi = deg2rad(-50);
            beta = beta_endo * (1 - depth) + beta_epi * depth;
            R_points = projectToRotatedBasis2D(points, beta); % 2xN
        else
            R_points = [ [ cos(orientation), sin(orientation)];
                         [-sin(orientation), cos(orientation)] ] * points'; % 2xN
        end

        % STEP 1:  MAIN FIBROSIS NOISE FIELD (O_b)
        P_f = [ R_points(1,:) / sqrt(align);
                R_points(2,:) * sqrt(align) ];

        O_b = Octave2D(P_f / feature_size, 4, roughness, perm_table, offsets1);

        % STEP 2: FIBRE-SELECTING FIELD (STRUCTURED PATTERN) (F)
        n_fibres_similarity = 4; wiggle_feature_length = 4; phasefield_strength = 5;
        phase_pts = [ R_points(1,:) / wiggle_feature_length;
                      R_points(2,:) / (n_fibres_similarity * fibre_sep) ];

        phasefield = Octave2D(phase_pts, 4, 0.5, perm_table2, offsets1);

        F = 0.5 + 0.5 * cos(2*pi * (R_points(2,:) / fibre_sep + phasefield_strength * (phasefield - 0.5)));
        F = F.^15;

        combined_noise = (1 - fibreness + fibreness * F) .* O_b;

        % STEP 3: DENSITY VARIATION FIELD (O_d)
        % Use original points (non-rotated) stretched by patch size for density variation field
        O_d = Octave2D(points' / patch_size, 4, 0.5, perm_table3, offsets1);

        combined_noise = combined_noise + patchiness * O_d;
    end

    % Ensure no values > 1 or < 0 (Safety clipping)
    % combined_noise = max(0, min(1, combined_noise));
    % (Optional, depending on how Octave2D/3D behaves)
end

function local_points = projectToRotatedBasis2D(points, beta)
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
