function assignment81( ~ )

% Import files from Jon Sporring
addpath('../jon/');

close all;

[I cmap] = imread('../../../images/starsDegraded.jpg');

% Adjust the image size to be used for scale.m
% Pixels are mirrored at the edges (this should really have it's own method)
[M N] = size(I);

if mod(size(I,1), 2)
    I = padarray(I, [1 0], 'replicate', 'post');
end
if mod(size(I,2), 2)
    I = padarray(I, [0 1], 'replicate', 'post');
end

numUnModified = CountStars(I);

fprintf('\n');
fprintf('No. of stars in unfiltered image = %d\n', numUnModified);
fprintf('\n');

% Inverse filtering with epsilon
epsilon = 0.2;
inv = InverseFilter(I, 2, epsilon, 0);
inv = uint8(inv);

numInvEpsilon = CountStars(inv);

fprintf('No. of stars in inverse filtered image (epsilon) = %d\n', numInvEpsilon);
fprintf('\n');

% Inverse filtering with Euler-Lagrange
fxx = real(ifft2(scale(fft2(I), 1.1, 0, 2)));
fyy = real(ifft2(scale(fft2(I), 1.1, 2, 0)));

fxy = double(I) - 2*(fxx+fyy);
fxy = uint8(fxy);

numInvEL = CountStars(fxy);

fprintf('No. of stars in inverse filtered image (Euler-Lagrange) = %d\n', numInvEL);
fprintf('\n');

% Inverse filtering with threshold
invt = uint8(InverseFilter(I, 2, 0.15, 1));

numInvThres = CountStars(invt);

fprintf('No. of stars in inverse filtered image (threshold) = %d\n', numInvThres);
fprintf('\n');

end
