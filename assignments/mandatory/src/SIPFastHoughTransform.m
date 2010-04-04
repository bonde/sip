function [H thetarange rhorange] = SIPFastHoughTransform(I, E, scale1, scale2)

    % Naive helper method
    function i = GetIndex(L, V)
        % Very cruel special case hack
        % Please forgive me
        if V == 90
            V = -90;
        end

        for k = 1:length(L)
            if L(k) == V
                i = k;
                return
            end
        end
        i = NaN;
    end

thetares = 0.25;
ntheta = 180/thetares;
thetarange = -90 : thetares : 90 - thetares;
tdiagonal = 90/thetares;

[M N] = size(I);

% Taken from MATLAB docs
rhores = 1;
D = sqrt((M - 1)^2 + (N - 1)^2);
rdiagonal = rhores*(ceil(D)/rhores);
nrho = 2*(ceil(D)/rhores) + 1;
rhorange = -rdiagonal: rhores : rdiagonal;

% Structure tensor from Jon
[u v] = AltTensor(I, scale1, scale2, 1, 0);

H = zeros([nrho ntheta]);

for m = 1:M
    for n = 1:N
        % If this is an edge
        if E(m,n)
            % Lookup the structure tensor and calculate theta
            U = u(m,n);
            V = v(m,n);
            theta = atan2(V, U);

            % Calculate rho
            rho = round(n*cos(theta) + m*sin(theta)) + rdiagonal;

            % Convert theta and lookup the index (slow) :/
            theta = rad2deg(theta);
            theta = round(theta/thetares)*thetares;
            theta = GetIndex(thetarange, theta);

            % Increment the accumulator if we have a real number
            if ~( isnan(rho) ) && ~( isnan(theta) )
                H(rho, theta) = H(rho, theta) + 1;
            end
        end
    end
end

end
