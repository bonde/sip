function G = AdaptiveFilter(S, factor, limit)

epsilon = 0.00000001;
padding = 10;

% Pad the original signal and apply a gaussian
S = ExpandSignal(S, padding);
GS = Gauss(S, 5/4);
N = length(S);

% Init target array
G = zeros(size(S));

% Get second derivative of blurred signal
Dx2 = ImDerivative(GS, 'dx', 's', 2);

% For each value we convolve with a different sigma
for n = 1:length(S)
    % Avoid division by zero
    val = max(abs(Dx2(n)), epsilon);

    % Decide if we should really blur this
    if val < limit
        % Find sigma
        sigma = 1/val*factor;

        % Convolve with original (padded) signal
        h = fspecial('gaussian', [1 7], sigma);
        tmp = conv(S, h, 'same');

        % Place value of interest in target array
        G(n) = tmp(n);
    else
        G(n) = S(n);
    end
end

% Remove the padding
G = G(:, padding+1:N-padding);

end
