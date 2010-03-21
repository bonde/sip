function assignment611( ~ )

close all;

% Use external package from http://staff.science.uva.nl/~rein/nldiffusionweb/material.html
addpath('../diffusion/');

A = 260;
B = 1200;
T_w = 40;
T_h = T_w;
theta = 90/2;

gridn = 16;

dimensions = 256;

s = zeros(dimensions,dimensions);
r = zeros(size(s));

for w = 1:dimensions
    for h = 1:dimensions
        s(w,h) = A*sin(w*2*pi/T_w + 45 + h*2*pi/T_h + 90);
        r(w,h) = B*sin(w*2*pi/T_w + 45 - h*2*pi/T_h + 90);
    end
end

super = s + r;

[u v Cc] = StructureTensor(s, 4, 4, gridn);

figure(12);
imshow(Cc, []);

[x y] = meshgrid(linspace(1, size(s,1), size(s,1)/gridn));

u = interp2(u, x, y);
v = interp2(v, x, y);

figure(1);
imshow(r, []);

figure(2);
imshow(super, []);

figure(3);
imshow(s, []);
hold on;
quiver(x-.5*u, y-.5*v, u, v);
hold off;

[u v Cc] = StructureTensor(super, 4, 4, gridn);

[x y] = meshgrid(linspace(1, size(s,1), size(s,1)/gridn));

u = interp2(u, x, y);
v = interp2(v, x, y);

figure(4);
imshow(super, []);
hold on;
quiver(x-.5*u, y-.5*v, u, v);
hold off;

figure(13);
imshow(Cc, [])

end
