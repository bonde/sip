function assignment621( ~ )

% Use external package from http://staff.science.uva.nl/~rein/nldiffusionweb/material.html
addpath('../diffusion/');
addpath('../jon/');

close all;

scale1 = 3;
scale2 = 5;
gridn = 4;

[I cmap] = imread('../../../images/FINGERPRINT.jpg');
%I = imread('../../../images/lenna.tiff');
%I = imread('../../../images/square.tiff');
[I cmap] = imread('../../../images/apple.tiff');
I = imread('../../../images/R1.tiff');

[M N] = size(I);

if mod(size(I,1), 2)
    I = padarray(I, [1 0], 'replicate', 'post');
end
if mod(size(I,2), 2)
    I = padarray(I, [0 1], 'replicate', 'post');
end

I = double(I);
%[u v] = StructureTensor(I, scale1, scale2, gridn, 1);
[u v] = AltTensor(I, scale1, scale2, gridn, 0);

[p q] = GradientVector(I, scale1, gridn);

I = gD(I, scale1, 0, 0);

maxmeshscale = max(size(I,1), size(I,2));

[x y] = meshgrid(linspace(1, maxmeshscale, maxmeshscale/gridn));

u = interp2(u, x, y);
v = interp2(v, x, y);

p = interp2(p, x, y);
q = interp2(q, x, y);

figure(1);
imshow(I, []);
hold on;
contour(I, [100 100], 'r');

% Note that u and v are swapped (and p and q for that matter)
quiver(x-.5*u, y-.5*v, u, v, 0, '.b');

quiver(x-.5*p, y-.5*q, p, q, 'g');
hold off;

%figure, imshow(Cc, []);

end
