function [u v] = Tensor(I, scale, varargin)

addpath('../diffusion/');

if nargin == 2
    sigma1 = 5;
    sigma2 = 4;
elseif nargin == 3
    sigma1 = varargin(1);
    sigma1 = sigma1{1};
    sigma2 = 4;
elseif nargin == 4
    sigma1 = varargin(1);
    sigma1 = sigma1{1};
    sigma2 = varargin(2);
    sigma2 = sigma2{1};
end

Dxx = gD(I, sigma1, 2, 0);
Dyy = gD(I, sigma1, 0, 2);
Dxy = gD(I, sigma1, 1, 1);

a = gD(Dxx, sigma2, 0, 0);
b = gD(Dxy, sigma2, 0, 0);
c = gD(Dyy, sigma2, 0, 0);

u = 2*b./(c - a + sqrt(c.^2 - 2*a.*c + a.^2+4*b.^2));
v = 2*b./(-c + a + sqrt(c.^2 - 2*a.*c + a.^2+4*b.^2));

normal = ones(size(u));
lu = sqrt(u.^2 + normal.^2);
lv = sqrt(v.^2 + normal.^2);

u = 0.75*scale*u./lu;
v = 0.75*scale*v./lv;

end
