%day = 670;
u = 1:30;
h0 = SW(data,day,1:30);
sen = 0.1*h0.result.sense;%0.003
h1 = SW(data,day,1:30,'method','implicit','rule','Sense','delta',sen,'norm','volatility');
%
[~,~,p,U,D,Q0,q0,H,ann_vec,ann_add] = getInitData(h1);
DltSq = h1.method.DeltaSq;
tld = diag(U'*D*(1+H*Q0*(Q0'*H*Q0\(p-q0))));
lambda = 200;
Dlt = lambda*tld*DltSq*tld;
QHQ = Q0'*H*Q0;
%
y = p-q0;
k = length(y);
T1 = inv(QHQ);
T2 = inv(Dlt);
m = 3e4;
alpha = k/2; beta1 = 1.5;
beta2 = 0.8;
arrtau1 = gamrnd(alpha,1/beta1,[m 1]);
arrtau2 = gamrnd(alpha,1/beta2,[m 1]);
%
for j=1:10
    beta1 = ones(1,m);
    beta2 = ones(1,m);
    alpha = 0 + k/2;
    parfor i = 1:m
        tau1 = arrtau1(i);
        tau2 = arrtau2(i);
        x1 = (tau1*T1+tau2*T2)\(tau2*T2*y) + chol(tau1*T1+tau2*T2)\randn(size(y));
        x2 = y-x1;
        beta1(i) = 0.1 + 0.5 * (x1'*T1*x1);
        beta2(i) = 0.1 + 0.5 * (x2'*T2*x2);   
    end
    ind = randi(m,[1,m]);
    arrtau1 = gamrnd(alpha,1./beta1(ind));
    ind = randi(m,[1,m]);
    arrtau2 = gamrnd(alpha,1./beta2(ind));
    fprintf('%3d: mean(tau1) = %4.2f; mean(tau2) = %4.2f \n',j,mean(arrtau1),mean(arrtau2))
end
% generate x1 sample
n = 2e5;
ind = randi(m,[1,n]);
arrtau1 = gamrnd(alpha,1./beta1(ind));
ind = randi(m,[1,n]);
arrtau2 = gamrnd(alpha,1./beta2(ind));
x1 = zeros(30,n);
parfor i = 1:n
    tau1 = arrtau1(i);
    tau2 = arrtau2(i);
    x1(:,i) = (tau1*T1+tau2*T2)\(tau2*T2*y) + chol(tau1*T1+tau2*T2)\randn(size(y));
    %x1(:,i) = (tau1*T1+tau2*T2)\(tau2*T2*y);
end
x2 = repmat(y,1,n) - x1;
%
h2 = h1;
h3 = h1;
h2.result.xi = Q0*(QHQ\mean(x1,2));
h3.result.xi = Q0*(QHQ\median(x1,2));
plotSW([h0 h1 h2 h3])
%% parametr density function
x = linspace(0,5e1,200);
[X,B1] = meshgrid(x,beta1);
pd1 = gampdf(X,alpha,1./B1);
%
[X,B2] = meshgrid(x,beta2);
pd2 = gampdf(X,alpha,1./B2);
hold on
plot(x,sum(pd1,1)/m,'r',x,sum(pd2,1)/m);
%% hist x2
N = 50;
rngs = linspace(-0.9e-2,0.9e-2,N);
cnts = zeros(30,N);
for j = 1:30
    cnts(j,:) = histc(x2(j,:),rngs);
end
%bar3(rngs,cnts',0.7)
[X,Y] = meshgrid(1:30,rngs);
plot3(X,Y,cnts','linewidth',2)
rotate3d on
%%
% h0 = SW(data,day,1:30);
% sen = 0.15*h0.result.sense;%0.003
% h1 = SW(data,day,1:30,'method','implicit','rule','Sense','delta',sen,'norm','volatility');
% plotSW(h0,'color',[0 0 0]);
% hold on;
% plotSW(h1,'time',60,'color',[1 0.4 0]);
% legend('location','southeast')
% title(h1.data.date)
% %% 
% tau1 = mean(arrtau1);
% tau2 = mean(arrtau2);
% x1 = (tau1*T1+tau2*T2)\(tau2*T2*y);
% h2 = h1;
% h2.result.xi = Q0*(QHQ\x1);
% plotSW([h1 h2])
