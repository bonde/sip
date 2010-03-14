function [vals locs] = EdgeDetect(S, threshold)

% Get first and second derivative
Dx = ImDerivative(S, 'dx', 'b', 1);
Dx2 = ImDerivative(S, 'dx', 'f', 2);

% Get zero crossings in second derivative
[C locs] = ZeroCrossings(Dx2);
vals = zeros(size(locs));

% Only chose values above a threshold
for j = 1:length(locs)
    if Dx(locs(j)) >= threshold
        vals(j) = S(locs(j));
    end
end

% Remove zero values in array
% This leaves only the edge(s) of interest
[vals locs] = PruneZeroes(vals, locs);

end
