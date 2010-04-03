function [H theta rhorange] = SIPHoughTransform(E)
% rho is equal to d from Jahne
% Set the resolution of rho and theta
% TODO: Make rho-resolution work
rhores = 1;
thetares = 0.25;
ntheta = 180/thetares;

[M N] = size(E);

% Taken from MATLAB docs about the Hough transform
D = sqrt((M - 1)^2 + (N - 1)^2);
diagonal = ceil(D/rhores);
nrho = 2*(ceil(D)/rhores) + 1

% Range of rho
rhorange = -diagonal: rhores : diagonal;

% Construct the Hough accumulator
H = zeros([nrho ntheta]);

% Initialize theta (in radians)
step = pi/ntheta;
theta = -pi/2 : step : pi/2 - step;

% For each pixel on an edge calculate rho (d) for every angle in theta
% and save the result in the Hough accumulator.
for m = 1:M
    for n = 1:N
        if E(m,n)
            for t = 1:ntheta
                rho = round(n*cos(theta(t)) + m*sin(theta(t)) );
                rho = round(rho*rhores);
                rho = rho + diagonal;
                H(rho, t) = H(rho, t) + 1;
            end
        end
    end
end

% Convert theta to degrees
theta = rad2deg(theta);

end
