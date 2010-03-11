function [vals locs] = EdgeDetect(S, thresshold)

Dx = ImDerivative(S, 'dx', 'f', 1);
D2x = ImDerivative(S, 'dx', 'f', 2);
[C locs] = ZeroCrossings(D2x);
vals = zeros(size(locs));

for j = 1:length(locs)
    %if sign(D2x(locs(j))) == -1 && Dx(locs(j)) >= thresshold
    if Dx(locs(j)) >= thresshold
        vals(j) = S(locs(j));
    end
end

[vals locs] = PruneZeroes(vals, locs);

end
