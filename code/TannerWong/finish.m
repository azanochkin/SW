%% Load data
load ./data/data_update_may2015.mat
data.Date = datetime(data.Date,'InputFormat','dd.MM.yyyy');
%% Set parameters
Ndt = 39;
% dt = data.Date(end-Ndt+1:end);
dt = data.Date(round(linspace(1,length(data.Date),Ndt)));
alpha = 0.01;
UFR = 0.0365;
kernel = 'SW';
mask = true(26,1);
mask(21:end) = false;
%
for it = 1:Ndt
    arrh(it) = SW(data,dt(it),'functional','new','mask',mask,...
        'norm','simple','alpha',alpha,'UFR',UFR,'kernel',kernel);
    if any(xor(mask,arrh(it).data.mask))
        warning(sprintf('in %i date mask not full',it))
    end
end
%
pTol = 1e-8;
M = 1e2;
N = 1e2;
m = 20;
n = 20;
%% Compute
coef = 4e-2;
alpha_x = coef * 0.5 * m * Ndt;
alpha_r = coef * 0.5 * Ndt;
beta_x = coef * 0.5 * m * Ndt * ones(1,N);
beta_r = coef * 0.5 * Ndt * ones(n,N);
%%
for l = 1:2000
    tau_x = gamrnd(alpha_x,1./beta_x);
    tau_r = gamrnd(alpha_r,1./beta_r);
    %
    norm_x = zeros(M,N);
    norm_r = zeros(n,M,N);
    norm_y = zeros(M,N);
    %
    for it = 1:Ndt
        [m,n,p,U,D,Q0,~,H] = getInitData(arrh(it));
        S = 1e-1*arrh(it).method.DeltaSq;
        for i = 1:N
            lS = tau_x(i)*S./tau_r(:,i);
            dx = zeros(m,1);
            dr = zeros(n,1);
            flag = true;
            while flag
                ex = exp(dx);
                Q = diag(ex)*(Q0 + D*U*diag(dr));
                P = diag(U'*D*ex);
                flag = max(abs(p'-sum(Q))) > pTol;
                if flag
                    beta = (Q'*H*Q + P'*lS*P)\(p - Q'*(1 - dx) + P' * dr);
                    dr = lS * P * beta;
                    dx = H * Q * beta;
                end
            end
            % A = 0.5*(inv(H)-diag(dx)+Q*inv(P*S*P)*Q') + Q*inv(P)*diag(beta)*U'*D*ex;
            % A = inv(A+A');
            A = H - H*Q*((Q'*H*Q+P'*lS*P)\(Q'*H));
            A = 0.5 * (A + A') / tau_x(i);
            L = chol(A,'lower');

            Y = randn(m,M);
            X = dx + L*Y;
            expX = exp(X);
            R = (p - Q0'*expX)./((U'*D)*expX);
            norm_x(:,i) = norm_x(:,i) + dot(X,H\X)';
            norm_r(:,:,i) = norm_r(:,:,i) + R.*(S\R);
            norm_y(:,i) = norm_y(:,i) + dot(Y,Y)';

    %         DX(:,i) = dx;%%%!!! prob != 1
    %         DR(:,i) = mean(R,2);%%!!! plob != 1
        end
    end
    nrt = squeeze(sum(norm_r.* repmat(permute(tau_r,[1 3 2]),1,M)));
    pr = norm_x.* tau_x + nrt - norm_y;
    pr = exp(-0.5*(pr - mean(pr)));
%     pr = M * pr./sum(pr);
    cumpr = cumsum(pr);
    ind = randi(N,[1,N]);
    rnd = rand([N,1]);
    for k = 1:N
        i = ind(k);
        j = find(cumpr(:,i) >= rnd(k)*cumpr(end,i),1,'first');
        beta_x(k) = 0 + 0.5 * norm_x(j,i);
        beta_r(:,k) = 0 + 0.5 * norm_r(:,j,i);
    end
    alpha_x = 0.5 * m * Ndt;
    alpha_r = 0.5 * Ndt;
    %%
    fprintf('%3i: sigma^2_x = %6.3f; sigma^2_r = %6.3f\n',l,mean(beta_x)/alpha_x,mean(beta_r(:))/alpha_r)
    if false
        cla
        x = linspace(0,3,300);
%         [X,BX] = meshgrid(x,beta_x);
%         pdx = gampdf(X,alpha_x,1./BX);    
%         plot(x,sum(pdx,1)/N,'r-.','DisplayName','tau_x');
        hold on
        for ii = 1:n
            [X,BR] = meshgrid(x,beta_r(ii,:));
            pdr = gampdf(X,alpha_r,1./BR);
            fr = sum(pdr,1)/N;
            [~,maxind] = max(fr);
            plot(x,fr,'DisplayName',sprintf('\\tau_{%i}',ii));
            text(x(maxind),fr(maxind)+0.01,sprintf('\\downarrow \\tau_{%i}',ii))
        end
        legend show
        hold off
        pause(0.1)
    end
end
%%
% T = 60;
% h1 = h;
% x = mean(DX,2);
% h1.result.xi = H\x;
% h1.result.r = h.method.r0 + (p - Q0'*exp(x))./(U'*D'*exp(x));
% %
% subplot(1,2,1);
% plotSW(h,'time',T,'color',[0 0 0]);
% ylim('manual');
% hold on;
% plotSW(h1,'time',T,'color',[1 0.4 0]);
% subplot(1,2,2);
% plotrates(h1);