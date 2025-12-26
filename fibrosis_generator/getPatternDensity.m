function [density, fibrosis_count_inside_fc] = getPatternDensity(presence, fibrotic_core)

    % Count how many elements are in the presence inside the fibrotic core
    fibrosis_count_inside_fc = sum(fibrotic_core(:) & presence(:));
    density = fibrosis_count_inside_fc/sum(fibrotic_core(:));

end
