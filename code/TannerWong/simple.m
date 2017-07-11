m = 5e5;
y = 10;
t1 = 1;
tau1 = gamrnd(1,1,[m 1]);
t2 = 1;
tau2 = gamrnd(1,1,[m 1]);
rng('shuffle');
for j = 1:20
    x1 = (tau2*y)./(tau1+tau2) + randn(m,1)./sqrt(tau1+tau2);
    x2 = y-x1;
    %x2 = (1./(tau1*t1+tau2*t2)).*(tau1*t1*y) + randn(1,m)./sqrt(tau1*t1+tau2*t2);
    %x1 = y-x2;
    %
    alpha = 0.5+0.5;
    beta1 = 0.5+(x1.*t1.*x1)/2;
    ind = randi(m,[m,1]);
    tau1 = t1*gamrnd(alpha,1./beta1(ind));
    %
    beta2 = 0.5+(x2.*t2.*x2)/2;
    ind = randi(m,[m,1]);
    tau2 = t2*gamrnd(alpha,1./beta2(ind));
    %
    fprintf('%3d: mean(x1) = %4.2f; mean(x2) = %4.2f \n',j,mean(x1),mean(x2))
    %hist([x1 x2],linspace(-20,30,50))
    %ginput();
end 
%% x1,x2 pdf 
nGrid = 200;
x = linspace(-20,30,nGrid);
pd1 = zeros(1,nGrid);
pd2 = zeros(1,nGrid);
%
sigmasq = tau1+tau2;
mu1 = (tau2*y)./sigmasq;
for i=1:nGrid
    t = x(i);
    pd1(i) = sum(sqrt(sigmasq).*exp(-(sigmasq/2).*(t - mu1).^2))/m;
    pd2(i) = sum(sqrt(sigmasq).*exp(-(sigmasq/2).*(t - (y-mu1)).^2))/m;
end
pd1 = pd1/sqrt(2*pi);
pd2 = pd2/sqrt(2*pi);
%loglog(0.01./x,sum(pd1,1)/m,'r',x,sum(pd2,1)/m,'b');
plot(x,pd1,'b',x,pd2,'r');
%% theta pdf
nGrid = 200;
x = linspace(0,0.5,nGrid);
pd1 = zeros(m,nGrid);
pd2 = zeros(m,nGrid);
%
G = gamma(alpha);
xdeg = (x.^(alpha-1));
for i=1:m
    pd1(i,:) = (beta1(i)^alpha)*xdeg.*exp(-beta1(i)*x)/G;
    pd2(i,:) = (beta2(i)^alpha)*xdeg.*exp(-beta2(i)*x)/G;
end
%loglog(0.01./x,sum(pd1,1)/m,'r',x,sum(pd2,1)/m,'b');
plot(x,sum(pd1,1)/m,'b',x,sum(pd2,1)/m,'r');
%% integration
%y = 0.01;
N = 30;
res = zeros(1,N);
yGrid = exp(linspace(log(1e-2),log(1e2),N));
for j = 1:N
    y = yGrid(j);
    n = 4e4;
    ind = randi(m,[n,1]);
    tau1 = gamrnd(alpha*ones(n,1),1./beta1(ind));
    x2 = y*tau1;
    %
    integ = 0;
    G = gamma(alpha);
    xdeg = (x2.^(alpha-1));
    for i=1:m
        integ = integ + sum((beta2(i)^alpha)*xdeg.*exp(-beta2(i)*x2)/G)/n;
    end
    integ = (integ*sqrt(1+(1/y)^2))/m;
    fprintf('y = %6.3f ; pdf(y) = %4.2e \n',y,integ)
    res(j) = integ;
end
plot(yGrid,res)
%%
n = 1e7;
ind = randi(m,[n,1]);
tau1 = gamrnd(alpha*ones(n,1),1./beta1(ind));
ind = randi(m,[n,1]);
tau2 = gamrnd(alpha*ones(n,1),1./beta2(ind));
res = tau1./tau2;
hist(res(res<10),1000)
%%
hist([tau1 ,tau2],y(1:3:end));
%hist([t1./tau1 ;t2./tau2]',50);
%hist(tau1./tau2)

