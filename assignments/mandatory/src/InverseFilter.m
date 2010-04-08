function Inv = InverseFilter(I, s, epsilon, threshold)

%% !! NOTE !!
% This is a copy/paste-hack from Jon Sporring's scale.m

% We use the Gaussian kernel that the program generates
% and deconvolve instead.

dr = 0;
dc = dr;

% First transform the input image
F = fft2(I);

% Compute the gaussian with specified sigma
if (s < 0)
  error('s must be larger than or equal 0');
else 
  if s == 0
    if (dr == 0) & (dc == 0)
      Is = I;
    else
      Is = zeros(size(I));
    end
  else
    rows = size(I,1);
    cols = size(I,2);

    G = zeros(rows,cols);
    if s == Inf
      G(1,1) = 1;
      Is = I.*G;
    else
      if((rows > 2 & rem(rows,2) ~= 0) | (cols > 2 & rem(cols,2) ~= 0))
	error('The image must have even side lengths');
      else
	% Calculate the Fourier transform of a gaussian fct.
	if (rows > 1) & (cols > 1) 
	  % 2 dimensional image
	  G(1:rows/2+1, 1:cols/2+1) = exp(-(repmat(([0:rows/2]'/rows).^2,[1,cols/2+1])+repmat(([0:cols/2]/cols).^2,[rows/2+1,1]))*(s*2*pi)^2/2);  
	  G(rows/2+1:rows, 1:cols/2+1) = flipud(G(2:rows/2+1, 1:cols/2+1));
	  G(1:rows/2+1, cols/2+1:cols) = fliplr(G(1:rows/2+1, 2:cols/2+1));
	  G(rows/2+1:rows, cols/2+1:cols) = fliplr(flipud(G(2:rows/2+1, 2:cols/2+1)));
	else
	  % 1 dimensional image
	  [val,ind] = max([rows,cols]);
	  G(1:val/2+1) = exp(-([0:val/2]'/val).^2*(s*2*pi)^2/2);
	  G(val/2+1:val) = flipdim(G(2:val/2+1),ind);
	end
	
	% Calculate the Differentiation matrix
	j = sqrt(-1);
	if (rows > 1) & (cols > 1)
	  x = [0:rows/2-1,-rows/2:-1]/rows;
	  y = [0:cols/2-1,-cols/2:-1]/cols;
	  DG = (x.^dr)'*(y.^dc)*(j*2*pi)^(dr+dc);
	else
	  if rows > 1
	    x = [0:rows/2-1,-rows/2:-1]/rows;
	    DG = (j*2*pi*x').^dr;
	  else
	    y = [0:cols/2-1,-cols/2:-1]/cols;
	    DG = (j*2*pi*y).^dc;
	  end
	end

    % Start hacking
    % Should we use a threshold or not...
    if threshold
        [M N] = size(I);
        for m = 1:M
            for n = 1:N
                if G(m,n) > epsilon
                    F(m,n) = F(m,n)/G(m,n);
                else
                    F(m,n) = 0;
                end
            end
        end
        F = F.*DG;
    else
        F = F./((epsilon + G).*DG);
    end
      end
    end
  end

Inv = real(ifft2(F));

end
