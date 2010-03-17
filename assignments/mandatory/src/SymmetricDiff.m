function D = SymmetricDiff(S)
[M N] = size(S);

% Pad signal
S = Pad(S, 1, 1);

D = zeros(size(S));

for n = 2:(N)
    D(n) = (S(n + 1) - S(n - 1))/2;
end

D = D(:, 2:N+1);

end
