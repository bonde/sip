function assignment621( ~ )

addpath('../diffusion/');

close all;

sigma1 = 3;
sigma2 = 4;
step = 4;

[I cmap] = imread('../../../images/FINGERPRINT.jpg');
%I = imread('../../../images/lenna.tiff');

I = double(I);
[u v] = Tensor(I, step, sigma1, sigma2);

I = gD(I, sigma1, 0, 0);

[x y] = meshgrid(linspace(1, size(u,1), size(u,1)/step));

u = interp2(u, x, y);
v = interp2(v, x, y);

figure(1);
imshow(I, []);
hold on;
contour(I, [100 100], 'r');
quiver(x-.5*v, y-.5*u, v, u);
hold off;

end
