function assignment81( ~ )

% Use external package from http://staff.science.uva.nl/~rein/nldiffusionweb/material.html
addpath('../jon/');
addpath('../diffusion/');

close all;

scale1 = 1;
scale2 = 2;
gridn = 4;

%[I cmap] = imread('../../../images/FINGERPRINT.jpg');
%I = imread('../../../images/lenna.tiff');
%I = imread('../../../images/R1.tiff');
%I = imread('../../../images/square.tiff');
%[I cmap] = imread('../../../images/apple.tiff');
[I cmap] = imread('../../../images/starsDegraded.jpg');

%I = imadjust(I);
%I = double(I);

% Adjust the image size to be used for scale.m
% Pixels are mirrored at the edges
[M N] = size(I);

if mod(size(I,1), 2)
    I = padarray(I, [1 0], 'replicate', 'post');
end
if mod(size(I,2), 2)
    I = padarray(I, [0 1], 'replicate', 'post');
end

numUnModified = CountStars(I, cmap);

fprintf('\n');
fprintf('No. of stars in unfiltered image = %d\n', numUnModified);
fprintf('\n');

% Inverse filtering with epsilon
epsilon = 0.2;
inv = InverseFilter(I, 2, epsilon, 0);
inv = uint8(inv);

numInvEpsilon = CountStars(inv, cmap);

fprintf('No. of stars in inverse filtered image (epsilon) = %d\n', numInvEpsilon);
fprintf('\n');

% Inverse filtering with Taylor series
fxx = real(ifft2(scale(fft2(I), 1.1, 0, 2)));
fyy = real(ifft2(scale(fft2(I), 1.1, 2, 0)));

fxy = double(I) - 2*(fxx+fyy);
fxy = uint8(fxy);

numInvTaylor = CountStars(fxy, cmap);

fprintf('No. of stars in inverse filtered image (Taylor) = %d\n', numInvTaylor);
fprintf('\n');

% Inverse filtering with threshold
invt = uint8(InverseFilter(I, 2, 0.15, 1));

numInvThres = CountStars(invt, cmap);

fprintf('No. of stars in inverse filtered image (threshold) = %d\n', numInvThres);
fprintf('\n');

end
