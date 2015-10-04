u = 1:30;
h0 = SW(data,day,u);
%
h1 = SW(data,day,u,'method','implicit','rule','Sense','delta',0.003);%0.1*h0.result.sense);
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