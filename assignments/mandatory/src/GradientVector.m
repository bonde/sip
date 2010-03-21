function [dx dy] = GradientVector(I, scale, gridn)
    dx = gD(I, scale, 0, 1);
    dy = gD(I, scale, 1, 0);

    % Normalize u and v (according to grid)
    normal = ones(size(dx));
    ldx = sqrt(dx.^2 + normal.^2);
    ldy = sqrt(dy.^2 + normal.^2);
    dx = 0.75*gridn*dx./ldx;
    dy = 0.75*gridn*dy./ldy;
end
