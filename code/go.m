% script require day and data to be set. Very usefull for using
% figure('keypressfcn','clf;day=day+1;go')
h0 = SW(data,day,1:20);
%
%sen = 0.15;%0.15*h0.result.sense;%0.003
%h1 = SW(data,day,1:20,'method','implicit','rule','Sense','delta',sen);
h1 = SW(data,day,1:20,'method','implicit','rule','l2','delta',sqrt(20),'norm','volatility');
%h1 = SW(data,day,1:20,'method','cauchy','rule','linf','delta',1,'norm','spread');
%
%subplot(2,2,1);
subplot(1,2,1);
plotSW(h0,'color',[0 0 0]);
hold on;
plotSW(h1,'time',60,'color',[1 0.4 0]);
legend('location','southeast')
title(h1.data.date)
%
%subplot(2,2,3)
subplot(1,2,2);
width1 = 0.8;
bar(1e-2*h1.result.grad/h1.result.annuity,width1,'FaceColor',[0.2,0.2,0.5],....
                     'EdgeColor','none','displayname','regularized');
hold on
width2 = width1/3;
bar(1e-2*h0.result.grad/h0.result.annuity,width2,'FaceColor',[0,0.7,0.7],...
                     'EdgeColor',[0,0.7,0.7],'displayname','original');
legend('location','southwest')
ylabel 'Sensitivity, %'
xlabel 'Term, years'
axis([-1 21 -1.25 1.25])
title(['Regul.sens = ',num2str(100*h1.result.sense,'%6.2f'),'%, or ', num2str(100*h1.result.sense/h0.result.sense,'%6.2f'),'% of SW sens.'])
%plotrates(h1);