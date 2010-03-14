function assignment522 ( ~ )

close all;

S = SignalWeek5();
S = ExpandSignal(S, 10);
GS = Gauss(S, 0.4);
S = S(11:74);
GS = GS(11:74);
A = AdaptiveFilter(S, 0.2, 0.27);

hold on;
range = [1:64] - 64;
plot(range, S, 'color', [1 0.4 0.7]);
plot(range, A);
%plot(range, S, 'b');
%plot(range, GS,'g');
xlabel('Depth in meters')
ylabel('Temperature in Degree Celsius')
hold off;

end
