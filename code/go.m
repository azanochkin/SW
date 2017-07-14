% script require day and data to be set. Very usefull for using
% figure('keypressfcn','clf;day=day+1;go')
%sen = 0.03*h0.result.sense;%0.003
%h1 = SW(data,day,1:30,'method','implicit','rule','Sense','delta',sen,'norm','volatility');
%h1 = SW(data,day,1:20,'method','implicit','rule','l2','delta',sqrt(20),'norm','volatility');
%h1 = SW(data,day,1:30,'method','cauchy','rule','linf','delta',1,'norm','spread');
%%
clf
prof = [1:100; ones(1,100)]';
T = 60;
u = 1:20;
alpha = 0.1;
UFR = log(1 + 0.042);
kernel = 'SW';
h0 = SW(data,day,u,...
    'functional','new','nsubiter',5,...
    'alpha',alpha,'UFR',UFR,'kernel',kernel,...
    'profile',prof);
h1 = SW(data,day,u,...
    'method','cauchy','rule','linf','delta',1,...
    'functional','new','nsubiter',3,'norm','spread',...
    'alpha',alpha,'UFR',UFR,'kernel',kernel,...
    'profile',prof);
%
subplot(4,2,[1,3,5]);
%subplot(1,2,1);
plotSW(h0,'time',T,'color',[0 0 0]);
ylim('manual');
hold on;
plotSW(h1,'time',T,'color',[1 0.4 0]);
%
legend('location','southeast')
title(h1.data.date)
%
subplot(4,2,7)
width1 = 0.8;
vec1 = 1e-2*h0.result.grad/h0.result.annuity;
b(1) = bar(abs(vec1),width1,'FaceColor',[0.2,0.2,0.5],'EdgeColor','none',...
    'displayname',[h0.method.name,' ',h0.method.functional]);
hold on
width2 = width1/3;
vec2 = 1e-2*h1.result.grad/h1.result.annuity;
b(2) = bar(abs(vec2),width2,'FaceColor',[0,0.7,0.7],'EdgeColor',[0,0.7,0.7],...
    'displayname',[h1.method.name,' ',h1.method.functional]);
plot(abs(vec1),'v','MarkerEdgeColor',[0.2,0.2,0.5],...
    'MarkerFaceColor','w','MarkerIndices',find(vec1<0))
plot(abs(vec1),'^','MarkerEdgeColor',[0.2,0.2,0.5],...
    'MarkerFaceColor','w','MarkerIndices',find(vec1>=0))
plot(abs(vec2),'v','MarkerEdgeColor',[0,0.7,0.7],...
    'MarkerFaceColor','w','MarkerIndices',find(vec2<0))
plot(abs(vec2),'^','MarkerEdgeColor',[0,0.7,0.7],...
    'MarkerFaceColor','w','MarkerIndices',find(vec2>=0))
set(gca,'YScale','log')
legend(b,'location','northwest')
ylabel 'Sensitivity, %'
xlabel 'Term, years'
xlim([0 31])
%axis([-1 31 -10.25 10.25])
title(sprintf('Regul.sens = %6.3f%%, or %6.3f%% of SW sens.',...
    100*h1.result.sense,100*h1.result.sense/h0.result.sense))
%
subplot(1,2,2);
plotrates(h1);