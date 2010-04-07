function G = FilterImage(I, d0, band, method)

    function h = MakeHighFrequencyEmphasisFilter(I, d0, alpha, beta)
        [im in] = size(I);

        h = zeros(size(I));

        for u = 1:im
            for v = 1:in
                D = sqrt((u - im/2)^2 + (v - in/2)^2);
                if D > 0
                    h(u, v) = alpha + beta*(1/(1 + (D/d0)^(2)));
                elseif D == 0;
                    h(u, v) = 1;
                end
            end
        end
    end

    function h = MakeIdealFilter(I, d0)
        [im in] = size(I);

        h = zeros(size(I));

        for u = 1:im
            for v = 1:in
                D = sqrt((u - im/2)^2 + (v - in/2)^2);
                if D < d0
                    h(u, v) = 1;
                end
            end
        end
    end

    function h = MakeButterworthFilter(I, d0, n)
        [im in] = size(I);

        h = zeros(size(I));

        for u = 1:im
            for v = 1:in
                D = sqrt((u - im/2)^2 + (v - in/2)^2);
                h(u, v) = 1/(1 + (D/d0)^(2*n));
            end
        end
    end

% Create filter
if strcmp(method, 'ideal')
    h = MakeIdealFilter(I, d0);
elseif strcmp(method, 'butter')
    n = 2;
    h = MakeButterworthFilter(I, d0, n);
elseif strcmp(method, 'emphasis')
    alpha = 40;
    beta = 2;
    h = MakeHighFrequencyEmphasisFilter(I, d0, alpha, beta);
else
    % Fail is fail
    error('No such method');
end

if strcmp(band, 'low')
    % Do nothing
elseif strcmp(band, 'high')
    h = (1 - h);
else
    % Fail
    error('Select valid band')
end

% Inspect the filter
%figure, imshow(h, []);

% Save the filter to disk
%imwrite(h, ['../report/images/' method '_' band '_filter_' num2str(d0) '.png'], 'png');

% Transform image
F = fft2(I);
F = fftshift(F);
%ftshow(F);

% Pointwise multiplication with the filter
G = h.*F;

% Inspect result (Fourier)
%ftshow(G);

% Reverse transform to spatial
G = abs(ifft2(fftshift(G)));

end
