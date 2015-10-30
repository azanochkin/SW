%day = 670;
u = 1:30;
h2 = SW(data,day,u,'method','implicit','rule','Sense','delta',0.003);
h1 = SW(data,day,u);
d = exp(-h1.method.ufr*u');
H = Wilson_Heart(0.1,u,u);
d1 = diag(d)*H*h1.result.xi;
d2 = diag(d)*H*h2.result.xi;
%
[ m,n,p,U,D,Q0,q0,H,ann_vec,ann_add] = getInitData( h1);
W = diag(d)*H*diag(d);
DltSq = h1.method.DeltaSq;
tld = diag(U'*D*(1+H*Q0*(Q0'*H*Q0\(p-q0))));
lambda = 5e2;
Dlt = lambda*tld*DltSq*tld;
QHQ = Q0'*H*Q0;
CY = QHQ*((QHQ+Dlt)\(p-q0));
CZ = Dlt*((QHQ+Dlt)\(p-q0));
% disp(CY'*inv(QHQ)*CY);
% disp((p-q0)'*inv(QHQ)*(p-q0));
%
y = p-q0;
k = length(y);
T1 = inv(QHQ);
T2 = inv(Dlt);
m = 2000;
alpha1 = 1; beta1 = 1;
alpha2 = 1; beta2 = 1;
arrtau1 = gamrnd(alpha1,1/beta1,[m 1]);
arrtau2 = gamrnd(alpha2,1/beta2,[m 1]);
%
for j=1:30
    alpha1 = ones(1,m); beta1 = ones(1,m);
    alpha2 = ones(1,m); beta2 = ones(1,m);
    for i = 1:m
        tau1 = arrtau1(i);
        tau2 = arrtau2(i);
        x1 = (tau1*T1+tau2*T2)\(tau2*T2*y) + chol(tau1*T1+tau2*T2)\randn(size(y));
        x2 = y-x1;
        alpha1(i) = 1+k/2;
        alpha2(i) = 1+k/2;
        beta1(i) = 1 + 0.5 * (tau1*x1'*T1*x1);
        beta2(i) = 1 + 0.5 * (tau2*x2'*T2*x2);
    end
    arrtau1 = sum(gamrnd(repmat(alpha1,m,1),repmat(1./beta1,m,1)),2)'/m;
    arrtau2 = sum(gamrnd(repmat(alpha2,m,1),repmat(1./beta2,m,1)),2)'/m;
    fprintf('%3d: tau1/tau2 = %4.2f \n',j,mode(arrtau1./arrtau2))
end
%%
tau1 = mode(arrtau1)
tau2 = mode(arrtau2)
x1 = (tau1*T1+tau2*T2)\(tau2*T2*y);
h0 = h1;
h1.result.xi = Q0*(QHQ\x1);
plotSW([h0 h1])
