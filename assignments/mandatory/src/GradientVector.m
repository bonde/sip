function [dx dy] = GradientVector(I, scale, sigma)
    dx = gD(I, sigma, 1, 0);
    dy = gD(I, sigma, 0, 1);
    %[dx dy] = gradient(I);

    % Normalize u and v (according to scale)
    normal = ones(size(dx));
    ldx = sqrt(dx.^2 + normal.^2);
    ldy = sqrt(dy.^2 + normal.^2);
    dx = 0.75*scale*dx./ldx;
    dy = 0.75*scale*dy./ldy;
end
