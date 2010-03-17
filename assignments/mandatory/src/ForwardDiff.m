function D = ForwardDiff(S)
[M N] = size(S);

S = Pad(S, 0, 1);
D = zeros(size(S));

for n = 1:size(S,2) - 1
    D(n) = S(n + 1) - S(n);
end

% Depad
D = D(:, 1:N);

end
