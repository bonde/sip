function D = ImDerivative(I, derivative, varargin)

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
    % Order-derivative in the x-direction
    for i = 1:order
        for m = 1:M
            D(m, :) = DiffMethod(I(m, :));
        end
        I = zeros(size(I));
        I(:,:) = D(:,:);
    end
    return
elseif strcmp(derivative, 'dy')
    % Order-derivative in the y-direction
    for i = 1:order
        for n = 1:N
            D(:, n) = transpose(DiffMethod(transpose(I(:, n))));
        end
        I = zeros(size(I));
        I(:,:) = D(:,:);
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
