function [H theta rhorange] = SIPHoughTransform(E)
% rho is equal to d from Jahne
rhores = 1;
ntheta = 180;

[M N] = size(E);

% Taken from MATLAB docs
D = sqrt((M - 1)^2 + (N - 1)^2);
diagonal = rhores*(ceil(D)/rhores);
nrho = 2*(ceil(D)/rhores) + 1;

% Range of rho
rhorange = -diagonal:diagonal;

H = zeros([nrho ntheta]);

step = pi/ntheta;
theta = -pi/2 : step : pi/2-step;

for m = 1:M
    for n = 1:N
        if E(m,n)
            for t = 1:ntheta
                rho = round(n*cos(theta(t)) + m*sin(theta(t)));
                H(rho + diagonal, t) = H(rho + diagonal, t) + 1;
            end
        end
    end
end

theta = rad2deg(theta);

end
