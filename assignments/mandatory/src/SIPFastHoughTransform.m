function [H theta rhorange] = SIPFastHoughTransform(I, E, scale1, scale2)

rhores = 1;
thetares = 1;
ntheta = 180/thetares;

[M N] = size(I);

% Taken from MATLAB docs
D = sqrt((M - 1)^2 + (N - 1)^2);
rdiagonal = rhores*(ceil(D)/rhores);
nrho = 2*(ceil(D)/rhores) + 1;
tdiagonal = 90;

%[u v] = StructureTensor(I, scale1, scale2, 1, 0);
[u v] = AltTensor(I, scale1, scale2, 1, 0);

H = zeros([nrho ntheta]);
t = size(H);

for m = 1:M
    for n = 1:N
        if E(m,n)
            U = u(m,n);
            V = v(m,n);
            theta = atan2(U, V)*thetares;

            rho = ceil(n*cos(theta) + m*sin(theta));
            theta = ceil(rad2deg(theta));

            H(rho + rdiagonal, theta + tdiagonal) = H(rho + rdiagonal, theta + tdiagonal) + 1;
        end
    end
end

tstep = pi/ntheta;
theta = -pi/2: tstep : pi/2 - tstep;
theta = rad2deg(theta);

rhorange = -rdiagonal: rhores : rdiagonal - rhores;

end
