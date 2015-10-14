h0 = SW(data,day,1:30);
%
sen = 0.1*h0.result.sense;
h1 = SW(data,day,1:30,'method','implicit','rule','Sense','delta',sen);
%
subplot(2,2,1);
plotSW(h0);
hold on;
plotSW(h1,'time',60,'color',[1 0.6 0]);
%
plotrates(h1);
subplot(2,2,3)
bar([[h0.result.grad],h1.result.grad])
title(['clear sense = ',num2str(h1.result.sense),' relev sense = ', num2str(h1.result.sense/h0.result.sense)])