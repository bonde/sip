function [u v] = StructureTensor(I, scale1, scale2, gridn)
% Awwwww, merge with Jon's code. Very ugly. Needs fix.

FI = fft2(I);

% Find the derivatives on the defined scale
%Dxx = gD(I, scale1, 0, 1).*gD(I, scale1, 0, 1);
%Dyy = gD(I, scale1, 1, 0).*gD(I, scale1, 1, 0);
%Dxy = gD(I, scale1, 0, 1).*gD(I, scale1, 1, 0);
Dy = real(ifft2(scale(FI,scale1,1,0)));
Dx = real(ifft2(scale(FI,scale1,0,1)));

% Allocate the structure tensor
%S = zeros([2 2 size(I)]);

% Convolve the derivatives with a gaussian
%S(1, 1, :, :) = gD(Dxx, scale2, 0, 0);
%S(1, 2, :, :) = gD(Dxy, scale2, 0, 0);
%S(2, 1, :, :) = S(1, 2, :, :);
%S(2, 2, :, :) = gD(Dyy, scale2, 0, 0);

% We produce the entrance of the structure tensor [[a,b],[b,c]]
%S(1, 1, :, :) = real(ifft2(scale(fft2(Dx.*Dx),scale2,0,0)));
%S(1, 2, :, :) = real(ifft2(scale(fft2(Dx.*Dy),scale2,0,0)));
%S(2, 1, :, :) = S(1, 2, :, :);
%S(2, 2, :, :) = real(ifft2(scale(fft2(Dy.*Dy),scale2,0,0)));
a = real(ifft2(scale(fft2(Dx.*Dx),scale2,0,0)));
b = real(ifft2(scale(fft2(Dx.*Dy),scale2,0,0)));
c = real(ifft2(scale(fft2(Dy.*Dy),scale2,0,0)));

% Allocate
%u = zeros(size(I));
%v = zeros(size(I));
%Cc = zeros(size(I));

%for m = 1:size(I,1)
%    for n = 1:size(I,2)
        % Construct the structure tensor
        %S = [a(m,n) b(m,n); b(m,n) c(m,n)];

        % Find the eigenvectors and -values
%        [e L] = eig(S(:,:,m,n));

        % Calculate the coherency
%        if L(1,1) + L(2,2) > 0
%            Cc(m,n) = ((L(1,1) - L(2,2))./(L(1,1) + L(2,2)))^2;
%        else
%            Cc(m,n) = 0;
%        end

        % Check which eigenvalue is dominant and set u and v accordingly
%        if L(1,1) > L(2,2)
%            u(m,n) = e(1,1);
%            v(m,n) = e(2,1);
%        elseif L(1,1) < L(2,2)
%            u(m,n) = e(1,2);
%            v(m,n) = e(2,2);
%        else
%            u(m,n) = 0;
%            v(m,n) = 0;
%        end
%    end
%end


d = c.^2-2*c.*a+a.^2+4*b.^2;
u = 2*b./(c-a+sqrt(d));
v = ones(size(u));

% Normalize u and v (according to grid)
l = sqrt(u.^2 + v.^2);
u = .75*gridn*u./l;
v = .75*gridn*v./l;

end
