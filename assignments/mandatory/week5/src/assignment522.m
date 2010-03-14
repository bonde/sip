function assignment522 ( ~ )

close all;

S = SignalWeek5();
A = AdaptiveFilter(S, 0.2);

hold on;
range = [1:64] - 64;
%plot(range, S, 'b');
plot(range, A, 'm');
xlabel('Depth in meters')
ylabel('Temperature in Degree Celsius')
hold off;

end
