% script require day and data to be set. Very usefull for using
% figure('keypressfcn','clf;day=day+1;go')
a = 0.1;
h0 = SW(data,day,1:30,'alpha',a);
%
sen = 0.1*h0.result.sense;%0.003
h1 = SW(data,day,1:30,'method','implicit','rule','Sense','delta',sen,'alpha',a);
%
subplot(2,2,1);
plotSW(h0);
hold on;
plotSW(h1,'time',60,'color',[1 0.6 0]);
%
plotrates(h1);
subplot(2,2,3)
width1 = 0.8;
bar(h1.result.grad,width1,'FaceColor',[0.2,0.2,0.5],....
                     'EdgeColor','none');
hold on
width2 = width1/3;
bar(h0.result.grad,width2,'FaceColor',[0,0.7,0.7],...
                     'EdgeColor',[0,0.7,0.7]);

title(['clear sense = ',num2str(h1.result.sense),' relev sense = ', num2str(h1.result.sense/h0.result.sense)])