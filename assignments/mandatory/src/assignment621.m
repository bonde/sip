function assignment621( ~ )

addpath('../diffusion/');

tensor = 0;
grad = 1;

close all;

sigma1 = 6;
sigma2 = 4;
step = 6;

[I cmap] = imread('../../../images/FINGERPRINT.jpg');
%I = imread('../../../images/lenna.tiff');

I = double(I);
if tensor
    [u v] = Tensor(I, step, sigma1, sigma2);
end

if grad
    p = gD(I, sigma1, 1, 0);
    q = gD(I, sigma1, 0, 1);

    % Normalize u and v (according to scale)
    normal = ones(size(p));
    lp = sqrt(p.^2 + normal.^2);
    lq = sqrt(q.^2 + normal.^2);
    u = 0.75*step*p./lp;
    v = 0.75*step*q./lq;
end

I = gD(I, sigma1, 0, 0);

%if grad
%    [p q] = gradient(I);
%end

[x y] = meshgrid(linspace(1, size(I,1), size(I,1)/step));

if tensor
    u = interp2(u, x, y);
    v = interp2(v, x, y);
end
if grad
    p = interp2(p, x, y);
    q = interp2(q, x, y);
end

figure(1);
imshow(I, []);
hold on;
contour(I, [110 110], 'r');
if tensor
    quiver(x-.5*v, y-.5*u, v, u);
end
if grad
    quiver(x-.5*p, y-.5*q, p, q, 'g');
end
hold off;

end
