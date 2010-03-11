function G = Pad(S, f, b)
    % Pad a 1D-signal with n in both directions
    [M N] = size(S);
    G = zeros([M N+f+b]);
    G(:, f+1:N+f) = S(:, :);
end
