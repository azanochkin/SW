m = 3000;
y = 10;
tau1 = gamrnd(1,1,[1 m]);
t1 = 1;
tau2 = gamrnd(1,1,[1 m]);
t2 = 1;
rng('shuffle');
for j = 1:300
    x1 = (1./(tau1*t1+tau2*t2)).*(tau2*t2*y) + randn(1,m)./sqrt(tau1*t1+tau2*t2);
    x2 = y-x1;
    %x2 = (1./(tau1*t1+tau2*t2)).*(tau1*t1*y) + randn(1,m)./sqrt(tau1*t1+tau2*t2);
    %x1 = y-x2;
    alpha = 1+1/2;
    beta1 = 1+0.5*(t1*x1.*tau1.*x1);
    tau1 = sum(gamrnd(alpha*ones(m),repmat(1./beta1,m,1)),2)'/m;
    beta2 = 1+0.5*(t2*x2.*tau2.*x2);
    tau2 = sum(gamrnd(alpha*ones(m),repmat(1./beta2,m,1)),2)'/m;
    fprintf('%3d: x1/x2 = %4.2f; tau1/tau2 = %4.2f \n',j,mean(x1./x2),mean(tau1./tau2))
end    
%%
hist([tau1 ;tau2]',20);
%hist([t1./tau1 ;t2./tau2]',50);
%hist(tau1./tau2)

