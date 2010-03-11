function D = BackwardDiff(S)
[M N] = size(S);

S = Pad(S, 1, 0);
D = zeros(size(S));

for n = 2:size(S,2)
    D(n) = S(n) - S(n - 1);
end

% Depad
D = D(:, 2:N+1);

end
