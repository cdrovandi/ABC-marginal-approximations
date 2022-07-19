function ssx = summStat_pair_corr(x, extra_args)
% the pair correlation summary statistics used in Browning et al (2018)
% Input: x - a 4 by 1 cell array containing location information for time
%            point 0, 12, 24 and 36 hours.
%        PC_dr - a vector of distances
% Output: ssx - summary statistics

    PC_dr = extra_args.PC_dr;

    if nargin < 2
        PC_dr = 50;
    end
    npc = length(PC_dr);
    ssx = zeros(1, 3*(1+npc));
    for i = 2:4
        ssx((i-2)*(npc+1)+1) = size(x{i},1);
        ssx((i-2)*(npc+1)+2 : (i-1)*(npc+1)) = arrayfun(@(y) pair_corr(x{i},y), PC_dr);
    end

end