function D = ImDerivative(I, derivative, varargin)
% Methods for finding derivatives in images/signals
% using the definitions from JÄhne p. 340 and 345
% i.e.
% Approximations for first derivatives:
%       -D_x  = g(x1, x2) - g(x1 - delta, x2)/delta                 [1. -1]
%       +D_x  = g(x1 + delta, x2) - g(x1, x2)/delta                 [1 -1.]
%        D2_x = g(x1 + delta, x2) - g(x1 - delta, x2)/2*delta       1/2 [1 0 -1]
%
% Second derivative:
%       D^2_x = -D_x +D_x      =       [1. -1] \star [1 -1.]    =   -D_x(+D_x(g)) = [1 -2  1]
%
% The Laplacian is then
%       L = D^2_x + D^2_y
% with filter mask
%                   [ 1 ]      [ 0  1  0 ]
% L = [1 -2 1 ]  +  [-2 ]   =  [ 1 -4  1 ]
%                   [ 1 ]      [ 0  1  0 ]

    function DiffMethod = SetMethod(method)
        if strcmp(method, 'b')
            DiffMethod = @BackwardDiff;
            return
        elseif strcmp(method, 'f')
            DiffMethod = @ForwardDiff;
            return
        elseif strcmp(method, 's')
            DiffMethod = @SymmetricDiff;
            return
        else
            error('Specify a method.');
            return
        end
    end

if nargin == 2
    DiffMethod = SetMethod('s');
    order = 1;
elseif nargin == 3;
    method = varargin(1);
    method = varargin{1};
    DiffMethod = SetMethod(method);
    order = 1;
elseif nargin > 3
    method = varargin(1);
    method = varargin{1};
    DiffMethod = SetMethod(method);
    order = varargin(2);
    order = order{1};
end


[M N] = size(I);
D = zeros(size(I));

if strcmp(derivative, 'dx')
    if order == 1
        for m = 1:M
            D(m, :) = DiffMethod(I(m, :));
        end
    elseif order == 2
        for m = 1:M
            forw = ForwardDiff(I(m, :));
            D(m, :) = BackwardDiff(forw);
        end
    else
        error('Not implemented. Only 1st and 2nd order derivatives.');
    end
    return
elseif strcmp(derivative, 'dy')
    if order == 1
        for n = 1:N
            D(:, n) = transpose(DiffMethod(transpose(I(:, n))));
        end
    elseif order == 2
        for n = 1:N
            forw = ForwardDiff(transpose(I(:, n)));
            D(:, n) = transpose(BackwardDiff(forw));
        end
    else
        error('Not implemented. Only 1st and 2nd order derivatives.');
    end
    return
elseif strcmp(derivative, 'sum')
    % The sum of the derivatives
    Dx = ImDerivative(I, 'dx', varargin);
    Dy = ImDerivative(I, 'dy', varargin);
    D = abs(Dx) + abs(Dy);
    return
elseif strcmp(derivative, 'mag')
    % Magnitude of derivatives
    Dx = ImDerivative(I, 'dx', varargin);
    Dy = ImDerivative(I, 'dy', varargin);
    D = sqrt(Dx .* Dx + Dy .* Dy);
    return
else
    error('Define a derivative "dx", "dy", "sum" or "mag"');
    return
end

end
