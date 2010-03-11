function [C locs] = ZeroCrossings(S)
prevsign = 0;
C = 0;

for n = 1:length(S)
    thissign = sign(S(n));
    if (thissign * prevsign) == -1
        C = C + 1;
        locs(C) = n;
    end
    if thissign ~= 0
        prevsign = thissign;
    end
end

end
