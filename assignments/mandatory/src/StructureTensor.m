function [u v Cc] = StructureTensorTensor(I, scale1, scale2, gridn)

% Find the derivatives on the defined scale
Dxx = gD(I, scale1, 0, 1).*gD(I, scale1, 0, 1);
Dyy = gD(I, scale1, 1, 0).*gD(I, scale1, 1, 0);
Dxy = gD(I, scale1, 0, 1).*gD(I, scale1, 1, 0);

% Convolve the derivatives with a gaussian
a = gD(Dxx, scale2, 0, 0);
b = gD(Dxy, scale2, 0, 0);
c = gD(Dyy, scale2, 0, 0);

% Allocate
u = zeros(size(I));
v = zeros(size(I));
Cc = zeros(size(I));

for m = 1:size(I,1)
    for n = 1:size(I,2)
        % Construct the structure tensor
        S = [a(m,n) b(m,n); b(m,n) c(m,n)];

        % Find the eigenvectors and -values
        [e L] = eig(S);

        % Calculate the coherency
        if L(1,1) + L(2,2) > 0
            Cc(m,n) = ((L(1,1) - L(2,2))./(L(1,1) + L(2,2)))^2;
        else
            Cc(m,n) = 0;
        end

        % Check which eigenvalue is dominant and set u and v accordingly
        if L(1,1) >= L(2,2)
            u(m,n) = e(1,1);
            v(m,n) = e(2,1);
        elseif L(1,1) < L(2,2)
            u(m,n) = e(1,2);
            v(m,n) = e(2,2);
        end
    end
end

% Normalize u and v (according to grid)
normal = ones(size(u));
lu = sqrt(u.^2 + normal.^2);
lv = sqrt(v.^2 + normal.^2);
u = 0.75*gridn*u./lu;
v = 0.75*gridn*v./lv;

end
