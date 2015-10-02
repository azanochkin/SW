main
%
h.method.sensefncname = 'S4';
h0 = original(h);
%
h.rule.name = 'Sense';
h.rule.delta = 0.1*h0.result.sense;
h1 = implicit(h);
%
subplot(2,2,1);
plotSW(h0);
hold on;
plotSW(h1,60,[1 0.6 0]);
%
plotrates(h1);
subplot(2,2,3)
bar([h0.result.grad,h1.result.grad])
title(['clear sense = ',num2str(h1.result.sense),' relev sense = ', num2str(h1.result.sense/h0.result.sense)])