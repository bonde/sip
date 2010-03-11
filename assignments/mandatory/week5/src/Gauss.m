function G = Gauss(I, sigma)
    % Compute the gaussian with specified sigma
    % Convolution is done in the frequency domain
    F = fft2(I);

    h = fspecial('gaussian', size(I), sigma);
    h = fft2(h);

    G = abs(fftshift(ifft2(h.*F)));
end
