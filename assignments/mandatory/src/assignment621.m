function assignment621( ~ )

% Use external package from http://staff.science.uva.nl/~rein/nldiffusionweb/material.html
addpath('../diffusion/');

close all;

sigma1 = 4;
sigma2 = 4;
step = 4;

[I cmap] = imread('../../../images/FINGERPRINT.jpg');
%I = imread('../../../images/lenna.tiff');

I = double(I);
[u v] = Tensor(I, step, sigma1, sigma2);

[p q] = GradientVector(I, step, sigma1);

I = gD(I, sigma1, 0, 0);

% Guess we have to scale or something
[x y] = meshgrid(linspace(1, size(I,1), size(I,1)/step));

% Interpolate (why I really don't know)
u = interp2(u, x, y);
v = interp2(v, x, y);

p = interp2(p, x, y);
q = interp2(q, x, y);

figure(1);
imshow(I, []);
hold on;
contour(I, [100 100], 'r');

% Note that u and v are swapped (and p and q for that matter)
quiver(x-.5*v, y-.5*u, v, u);

quiver(x-.5*q, y-.5*p, q, p, 'g');
hold off;

end
