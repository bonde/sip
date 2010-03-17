function S = GenerateSignal( ~ )

randn('seed', 31);

N = 64;
S = [zeros(1,20),ones(1,3),zeros(1,2),ones(1,5)];
S = [S,zeros(1,N-length(S))] + 0.1*randn(1,N);
S = cumsum(S);

end
