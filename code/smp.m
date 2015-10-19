day = 2359;%2359
color = [0 0 0];
color(rem(day,3)+1) = 0.9;
%
h0 = SW(data,day,1:20);
subplot(1,2,1);
plotSW(h0,'color',color);
lg = legend('location','southeast');
set(lg,'fontsize',8);
title 'Original'
hold on
%
sen = 0.15*h0.result.sense;%0.003
h1 = SW(data,day,u,'method','implicit','rule','Sense','delta',sen,'norm','volatility');
%h1 = SW(data,day,1:20,'method','implicit','rule','l2','delta',sqrt(20),'norm','volatility');
subplot(1,2,2);
plotSW(h1,'color',color);
lg = legend('location','southeast');
set(lg,'fontsize',8);
title 'Regularized, 20 years'
hold on
%%
day = 2359;%2359
color = [0 0 0];
color(rem(day,3)+1) = 0.9;
%
sen = 0.003;
h1 = SW(data,day,1:30,'method','implicit','rule','Sense','delta',sen,'norm','volatility');
plotSW(h1,'color',color);
lg = legend('location','southeast');
set(lg,'fontsize',8);
title 'Regularized, 30 years'
hold on