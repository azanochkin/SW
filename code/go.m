% script require day and data to be set. Very usefull for using
% figure('keypressfcn','clf;day=day+1;go')
h0 = SW(data,day,1:30);
%
sen = 0.03*h0.result.sense;%0.003
%h1 = SW(data,day,1:30,'method','implicit','rule','Sense','delta',sen,'norm','volatility');
%h1 = SW(data,day,1:20,'method','implicit','rule','l2','delta',sqrt(20),'norm','volatility');
h1 = SW(data,day,1:30,'method','cauchy','rule','linf','delta',1,'norm','spread');
%%
h0 = SW(data,day,1:30,'functional','new','nsubiter',5,'alpha',0.01,'UFR',0.01);
h2 = SW(data,day,1:30,'alpha',0.01,'UFR',0.01);
h1 = SW(data,day,1:30,'method','iterative','rule','l2','delta',sqrt(20),...
    'functional','new','nsubiter',5,'alpha',0.01,'UFR',0.01);
%%
subplot(4,2,[1,3,5]);
%subplot(1,2,1);
plotSW(h0,'time',100,'color',[0 0 0]);
ylim('manual');
hold on;
plotSW(h1,'time',100,'color',[1 0.4 0]);
plotSW(h2,'time',100,'color',[0.4 0.4 1]);
%
%plotSW(SW(data,day,1:30,'method','implicit','rule','Sense','delta',0.15*h0.result.sense,'norm','volatility'),'time',60,'color',[0 0.9 0.2]);
legend('location','southeast')
title(h1.data.date)
%%
subplot(4,2,7)
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
axis([-1 31 -0.15 0.15])
title(['Regul.sens = ',num2str(100*h1.result.sense,'%6.2f'),'%, or ', num2str(100*h1.result.sense/h0.result.sense,'%6.2f'),'% of SW sens.'])
%%
subplot(1,2,2);
plotrates(h1);