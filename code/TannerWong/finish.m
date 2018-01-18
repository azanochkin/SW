%% Load data
load ./data/data_update_may2015.mat
%% Set parameters
dt='04.07.2016';
alpha = 0.01;
UFR = 0.0365;
kernel = 'SW';
mask = true(26,1);
mask(21:end) = false;
%
h = SW(data,dt,'method','Tikhonov','functional','new','mask',mask,...
    'norm','simple','alpha',alpha,'UFR',UFR,'kernel',kernel);
%% Compute
[m,n,p,U,D,Q0,q0,H] = getInitData(h);
S = h.method.DeltaSq;
pTol = 1e-8;
M = 5e2;
N = 10e3;
%%
alpha_x = 0.5 * m;
alpha_r = 0.5 * n;
beta_x = 0.5 * m * ones(1,N);
beta_r = 0.5 * n * ones(1,N);
%%
norm_x = nan(M,N);
norm_r = nan(M,N);
norm_y = nan(M,N);
DX = nan(m,N);
DR = nan(n,N);
for l = 1:50
    disp(l)
    if fasle
        x = linspace(0,10,300);
        [X,BX] = meshgrid(x,beta_x);
        pdx = gampdf(X,alpha_x,1./BX);    
        [X,BR] = meshgrid(x,beta_r);
        pdr = gampdf(X,alpha_r,1./BR);
        plot(x,sum(pdx,1)/N,'r',x,sum(pdr,1)/N);
        ylim([0, 1.8])
    end
    %
    tau_x = gamrnd(alpha_x,1./beta_x);
    tau_r = gamrnd(alpha_r,1./beta_r);
    dx_init = mean(DX,2);
    dr_init = (p - Q0'*exp(dx_init))./(U'*D'*exp(dx_init));
    %
    for i = 1:N
        lambda = tau_x(i)/tau_r(i);
        dx = zeros(m,1);
        dr = zeros(n,1);
        flag = true;
        while flag
            ex = exp(dx);
            Q = diag(ex)*(Q0 + D*U*diag(dr));
            P = diag(U'*D*ex);
            flag = max(abs(p'-sum(Q))) > pTol;
            if flag
                beta = (Q'*H*Q + lambda*P'*S*P)\(p - Q'*(1 - dx) + P' * dr);
                dr = lambda * S * P * beta;
                dx = H * Q * beta;
            end
        end
        % A = 0.5*(inv(H)-diag(dx)+Q*inv(P*S*P)*Q') + Q*inv(P)*diag(beta)*U'*D*ex;
        % A = inv(A+A');
        A = H - H*Q*((Q'*H*Q+lambda*P'*S*P)\(Q'*H));
        A = 0.5 * (A + A') / tau_x(i);
        L = chol(A,'lower');

        Y = randn(m,M);
        X = repmat(dx,1,M) + L*Y;
        R = (p - Q0'*exp(X))./(U'*D'*exp(X));
        norm_x(:,i) = dot(X,H\X);
        norm_r(:,i) = dot(R,S\R);
        norm_y(:,i) = dot(Y,Y);
        
        DX(:,i) = dx;%%%!!! prob != 1
        DR(:,i) = mean(R,2);%%!!! plob != 1
    end
    norm_x = norm_x .* tau_x;
    norm_r = norm_r .* tau_r;
    pr = exp(-0.5*(norm_x + norm_r - norm_y));
%     pr = M * pr./sum(pr);
    cumpr = cumsum(pr);
    ind = randi(N,[1,N]);
    rnd = rand([N,1]);
    for k = 1:N
        i = ind(k);
        j = find(cumpr(:,i) >= rnd(k)*cumpr(end,i),1,'first');
        beta_x(k) = 0 + 0.5 * norm_x(j,i);
        beta_r(k) = 0 + 0.5 * norm_r(j,i);
    end
end
%%
T = 60;
h1 = h;
x = mean(DX,2);
h1.result.xi = H\x;
h1.result.r = h.method.r0 + (p - Q0'*exp(x))./(U'*D'*exp(x));
%
subplot(1,2,1);
plotSW(h,'time',T,'color',[0 0 0]);
ylim('manual');
hold on;
plotSW(h1,'time',T,'color',[1 0.4 0]);
subplot(1,2,2);
plotrates(h1);