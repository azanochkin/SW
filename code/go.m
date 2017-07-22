%% Load data
load ./data/data_update_may2015.mat
%% Set parameters
clf
dt='14.07.2016';
T = 60;
alpha = 0.01;
tau = 10e-4;
convpnt = 60;
UFR = 0.0365;
kernel = 'SW';
prof = [1:100; ((1:100))]';
%
h0 = SW(data,dt,...
    'functional','new','nsubiter',5,...
    'alpha',alpha,'tau',tau,'convpnt',convpnt,...
    'UFR',UFR,'kernel',kernel,...
    'profile',prof);
h1 = SW(data,dt,...
    'method','cauchy','rule','linf','delta',1,...
    'functional','new','nsubiter',3,'norm','spread',...
    'alpha',h0.method.alpha,'tau',tau,'convpnt',convpnt,...
    'UFR',UFR,'kernel',kernel,...
    'profile',prof);
%
subplot(4,2,[1,3,5]);
%subplot(1,2,1);
plotSW(h0,'time',T,'color',[0 0 0]);
ylim('manual');
hold on;
plotSW(h1,'time',T,'color',[1 0.4 0]);
area([0 T],1e2*(UFR+tau)*[1 1],1e2*(UFR-tau),...
    'EdgeColor','none','facealpha',0.4,'ShowBaseLine','off',...
    'displayname','convergence interval')
%
legend('location','southeast')
title(h1.data.date)
%
subplot(4,2,7)
width1 = 0.8;
vec1 = 1e-2*h0.result.grad/h0.result.annuity;
b(1) = bar(h0.data.tenor,vec1,width1,'FaceColor',[0 0 0],'EdgeColor','none',...
    'displayname',[h0.method.name,' ',h0.method.functional]);
hold on
width2 = width1/3;
vec2 = 1e-2*h1.result.grad/h1.result.annuity;
b(2) = bar(h1.data.tenor,vec2,width2,'FaceColor',[1 0.4 0],'EdgeColor',[0.7 0.7 0 ],...
    'displayname',[h1.method.name,' ',h1.method.functional]);
%set(gca,'YScale','log')
legend('location','northwest')
ylabel 'Sensitivity, %'
xlabel 'Term, years'
%xlim([0 31])
% axis([-1 31 -10.25 10.25])
% title(sprintf('Regul.sens = %6.3f%%, or %6.3f%% of SW sens.',...
%     100*h1.result.sense,100*h1.result.sense/h0.result.sense))
%
subplot(1,2,2);
plotrates(h1);