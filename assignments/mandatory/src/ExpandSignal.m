function G = ExpandSignal(S, n)

[M N] = size(S);

G = Pad(S, n, n);

G(1, 1:n) = S(1,1);
G(1, n+N+1:N+2*n) = S(1,N);

end
