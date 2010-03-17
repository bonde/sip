function [nVals nLocs] = PruneZeroes(vals, locs)

    function C = CountPositive(S)
        C = 0;
        for i = 1:length(S)
            if sign(S(i)) > 0
                C = C + 1;
            end
        end
    end

nVals = zeros(CountPositive(vals), 1);
nLocs = zeros(size(nVals));

j = 1;
for i = 1:length(vals)
    if sign(vals(i)) == 1
        nVals(j) = vals(i);
        nLocs(j) = locs(i);
        j = j + 1;
    end
end
end
