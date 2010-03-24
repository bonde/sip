function assignment621( ~ )

% Use external package from http://staff.science.uva.nl/~rein/nldiffusionweb/material.html
addpath('../jon/');
addpath('../diffusion/');

close all;

scale1 = 3;
scale2 = 3;
gridn = 4;

%[I cmap] = imread('../../../images/FINGERPRINT.jpg');
%I = imread('../../../images/lenna.tiff');
I = imread('../../../images/square.tiff');
[I cmap] = imread('../../../images/apple.tiff');

% Adjust the image size to be used for scale.m
% Pixels are mirrored at the edges
[M N] = size(I);

if mod(size(I,1), 2)
    I = padarray(I, [1 0], 'replicate', 'post');
end
if mod(size(I,2), 2)
    I = padarray(I, [0 1], 'replicate', 'post');
end

E = edge(I, 'canny', 0.5, 1.4);
[H2 th rho] = hough(E);

%[u v] = StructureTensor(I, scale1, scale2, gridn);
[H T R] = SIPHoughTransform(E);

%E = gD(double(E), scale1, 0, 0);
%E = gD(double(I), scale1, 0, 0);

%[x y] = meshgrid(linspace(1, size(I,2), size(I,2)/gridn));

%u = interp2(u, x, y);
%v = interp2(v, x, y);

%p = interp2(p, x, y);
%q = interp2(q, x, y);

P  = houghpeaks(H, 5, 'threshold', ceil(0.3*max(H(:))));
%P  = houghpeaks(H, 23, 'threshold', 12);
lines = houghlines(E, T, R, P, 'FillGap', 5, 'MinLength', 7);

figure(1);
imshow(I, []);
hold on;
max_len = 0;
for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];
       plot(xy(:,1),xy(:,2), 'LineWidth', 2, 'Color', 'red');

       % Determine the endpoints of the longest line segment
       len = norm(lines(k).point1 - lines(k).point2);
       if ( len > max_len)
           max_len = len;
           xy_long = xy;
       end
end
% highlight the longest line segment
plot(xy_long(:,1),xy_long(:,2), 'LineWidth', 2, 'Color', 'green');
%contour(I, [100 100], 'r');

% Note that u and v are swapped (and p and q for that matter)
%quiver(x-.5*u, y-.5*v, u, v, 0, '.b');

%quiver(x-.5*p, y-.5*q, p, q, 'g');
hold off;

figure, imshow(H2, []);
figure, imshow(H, []);

end
