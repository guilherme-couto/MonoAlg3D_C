function [permute_table, offset_table] = generateTables(seed)
% This function generates the permutation and offset tables for the provided seed.
%
% INPUTS:
%
% seed - the seed to be used for the random number generator
%
% OUTPUTS:
%
% permute_table - the permutation table generated
% offset_table - the offset table generated

% Set the seed for the random number generator
rng(seed);

% Assume a decent safe number like eight for the number of offsets
N_freqs = 8;

% Permutation tables for this seed
for j = 1:N_freqs
    permute_table(j,:) = int32(randperm(256) - 1);
end

% Offset table for this seed
offset_table = rand(N_freqs, 2) - 0.5;

end
