function [u v] = Tensor(I, scale, sigma1, sigma2)

% Use external package from http://staff.science.uva.nl/~rein/nldiffusionweb/material.html
addpath('../diffusion/');

% Find the derivatives on the scale defined by sigma1
Dxx = gD(I, sigma1, 2, 0);
Dyy = gD(I, sigma1, 0, 2);
Dxy = gD(I, sigma1, 1, 1);

% Convolve the derivatives with a gaussian
a = gD(Dxx, sigma2, 0, 0);
b = gD(Dxy, sigma2, 0, 0);
c = gD(Dyy, sigma2, 0, 0);

% Get vectors u and v
u = 2*b./(c - a + sqrt(c.^2 - 2*a.*c + a.^2+4*b.^2));
v = 2*b./(-c + a + sqrt(c.^2 - 2*a.*c + a.^2+4*b.^2));

% Normalize u and v (according to scale)
normal = ones(size(u));
lu = sqrt(u.^2 + normal.^2);
lv = sqrt(v.^2 + normal.^2);
u = 0.75*scale*u./lu;
v = 0.75*scale*v./lv;

end
