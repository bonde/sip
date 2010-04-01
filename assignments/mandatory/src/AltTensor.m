function [u v] = AltTensor(I, scale1, scale2, gridn, orientation)

FI = fft2(I);

% Finally we calulcate the eigenvectors of the structure tensor
% We measure the derivatives
Ly = real(ifft2(scale(FI,scale1,1,0)));
Lx = real(ifft2(scale(FI,scale1,0,1)));
% We produce the entrance of the structure tensor [[a,b],[b,c]] 
a = real(ifft2(scale(fft2(Lx.*Lx),scale2,0,0)));
b = real(ifft2(scale(fft2(Lx.*Ly),scale2,0,0)));
c = real(ifft2(scale(fft2(Ly.*Ly),scale2,0,0)));
% We calulate the eigen values e and vector v.
d = c.^2-2*c.*a+a.^2+4*b.^2;
e = 0.5*(c+a+sqrt(d));

if orientation
    u = 2*b./(-c+a+sqrt(d));
    v = ones(size(u));
    %v = 2*b./(c-a+sqrt(d));
    %u = ones(size(v));
else
    u = 2*b./(c-a+sqrt(d));
    v = ones(size(u));
end

l = sqrt(u.^2+v.^2);
u = .75*gridn*u./l;
v = .75*gridn*v./l;

end
