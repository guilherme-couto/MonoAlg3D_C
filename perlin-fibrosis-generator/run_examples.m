clc;
fprintf("Initializing GNT Framework...\n\n");

% --- SCENARIO 1: 2D Diffuse Full ---
% Domain: dx=0.01, 4x4 cm
scn1.type = 'diffuse';
scn1.dim_mode = '2D';
scn1.shape = 'full';
scn1.angle = 45;
scn1.density = 0.25;
% Format: [dx, Lx, Ly]
scn1.domain = [0.01, 4.0, 4.0];
% Core dims ignored for 'full'
scn1.core   = [];
scn1.desc   = [scn1.type, ' ', scn1.shape, ' ', scn1.dim_mode];

% --- SCENARIO 2: 2D Interstitial Ellipse (Embedded) ---
% Domain: dx=0.01, 2x4 cm
scn2.type = 'interstitial';
scn2.dim_mode = '2D';
scn2.shape = 'ellipse';
scn2.angle = 60;
scn2.density = 0.1;
% Format: [dx, Lx, Ly]
scn2.domain = [0.01, 2.0, 4.0];
% Format: [width, height]
scn2.core   = [1.5, 3.0];
scn2.desc   = [scn2.type, ' ', scn2.shape, ' ', scn2.dim_mode];

% --- SCENARIO 3: 3D Patchy Full ---
% Domain: dx=0.01, 1x1x1 cm
scn3.type = 'patchy';
scn3.dim_mode = '3D';
scn3.shape = 'full';
scn3.angle = [0, 0]; % [Phi, Theta]
scn3.density = 0.2;
% Format: [dx, Lx, Ly, Lz]
scn3.domain = [0.01, 1.0, 1.0, 1.0];
% Core dims ignored for 'full'
scn3.core   = [];
scn3.desc   = [scn3.type, ' ', scn3.shape, ' ', scn3.dim_mode];

% --- SCENARIO 4: 3D Compact Box (Embedded) ---
% Domain: dx=0.01, 2x1x0.5 cm
scn4.type = 'compact';
scn4.dim_mode = '3D';
scn4.shape = 'box';
scn4.angle = [0, 0];
scn4.density = 0.3;
scn4.domain = [0.01, 2.0, 1.0, 0.5];
% Core: 1.5x0.8x0.3 cm centered box
scn4.core   = [1.5, 0.8, 0.3];
scn4.desc   = [scn4.type, ' ', scn4.shape, ' ', scn4.dim_mode];

test_scenarios = {scn1, scn2, scn3, scn4};
seed = 42;
save_mesh = true;
save_figure = true;

% LOAD OR INSTALL IMAGE PACKAGE DEPENDENCY
fprintf('Checking dependencies...\n');
try
    pkg load image;
    fprintf('  [OK] Image package loaded.\n\n');
catch
    fprintf('  [FAIL] Image package not found. Installing...\n');
    try
        pkg install -forge image;
        pkg load image;
        fprintf('  [OK] Image package installed and loaded.\n\n');
    catch
        error('Failed to install/load image package. Please check your Octave setup.');
    end
end

% --- EXECUTION LOOP ---
for i = 1:length(test_scenarios)
    s = test_scenarios{i};
    fname = sprintf('examples/%s_%s_%s', s.type, s.shape, s.dim_mode);

    fprintf('>>> Running Scenario %d: %s\n', i, s.desc);
    try
        % func(type, dens, seed, angle, dim_mode, domain_dims, shape, core_dims, out, mesh, fig)
        run_fibrosis_generator(s.type, s.density, seed, s.angle, s.dim_mode, ...
                               s.domain, s.shape, s.core, fname, save_mesh, save_figure);
        fprintf('  [OK] Success.\n');
    catch err
        fprintf(2, '  [FAIL] Error: %s\n', err.message);
        for k=1:length(err.stack), fprintf('      %s: %d\n', err.stack(k).name, err.stack(k).line); end
    end
    fprintf('\n');
end
