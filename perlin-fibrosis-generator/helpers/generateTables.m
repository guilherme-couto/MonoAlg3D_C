function [permute_table, offset_table1, offset_table2] = generateTables(seed, is3D)
% Generates permutation and offset tables. Supports both 2D and 3D.
%
% INPUTS:
%   seed - RNG seed
%   is3D - (Optional) Boolean. If true, generates Nx3 offsets. Default: false.

    if nargin < 2
        is3D = false;
    end

    % Set the seed
    rng(seed);

    % Number of octaves/frequencies
    N_freqs = 8;

    % Permutation table (Same for 2D/3D)
    permute_table = zeros(N_freqs, 256, 'int32');
    for j = 1:N_freqs
        permute_table(j,:) = int32(randperm(256) - 1);
    end

    % Offset dimensions
    if is3D
        num_cols = 3;
    else
        num_cols = 2;
    end

    % Offset tables
    offset_table1 = rand(N_freqs, num_cols) - 0.5;

    % Only generate second table if requested (used in 3D structure logic)
    if nargout > 2
        offset_table2 = rand(N_freqs, num_cols) - 0.5;
    else
        offset_table2 = [];
    end
end
