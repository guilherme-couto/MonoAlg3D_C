function threshold = getCompositionThreshold(density, tolerance)
% This function returns the threshold to be used for the density of the pattern.
%
% INPUTS:
%
% density - the density of the pattern to be generated
% tolerance - the tolerance for density matching
%
% OUTPUTS:
%
% threshold - the threshold to be used for the density of the pattern

if density >= 0.1 - tolerance
    threshold = 0.1;
else
    threshold = density;
end

end
