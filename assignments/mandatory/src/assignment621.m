function assignment621( ~ )

% Use external package from http://staff.science.uva.nl/~rein/nldiffusionweb/material.html
addpath('../diffusion/');

close all;

scale1 = 4;
scale2 = 4;
gridn = 4;

[I cmap] = imread('../../../images/FINGERPRINT.jpg');
%I = imread('../../../images/lenna.tiff');
%I = imread('../../../images/square.tiff');

I = double(I);
[u v Cc] = StructureTensor(I, scale1, scale2, gridn);

[p q] = GradientVector(I, scale1, gridn);

I = gD(I, scale1, 0, 0);

[x y] = meshgrid(linspace(1, size(I,1), size(I,1)/gridn));

u = interp2(u, x, y);
v = interp2(v, x, y);

p = interp2(p, x, y);
q = interp2(q, x, y);

figure(1);
imshow(I, []);
hold on;
contour(I, [100 100], 'r');

% Note that u and v are swapped (and p and q for that matter)
quiver(x-.5*u, y-.5*v, u, v);

quiver(x-.5*p, y-.5*q, p, q, 'g');
hold off;

%figure, imshow(Cc, []);

end
