function h = iterative_new(h)
    function [dr,xi,dxidr,dxidp] = fun(lambda)
        Sigma = lambda*DltSq;
        xi = zeros(m,1);
        dxidr = zeros(m,n);
        dxidp = zeros(m,n);
        for i = 1:nIter
            dxi = zeros(m,1);
            dr = zeros(n,1);
            for j = 1:nsubiter
                eHxi = exp(H*(xi+dxi));
                Q = diag(eHxi)*(Q0 + D*U*diag(dr));
                P = diag(U'*D*eHxi);
                beta = (Q'*H*Q+P'*Sigma*P)\(p - Q'*(1 - H*dxi) + P' * dr);
                dr = Sigma * P * beta;
                dxi = Q*beta;
            end
            xi = xi + dxi;
            %
            if isfndderiv
                N = diag(eHxi)*D*U*diag(beta);
                A = [ H - H*diag(dxi)*H,     -H*Q  , -H*N*Sigma;...
                     -(H*Q)'           , zeros(n)  ,   -P*Sigma;...
                     -(H*N*Sigma)'     ,-(P*Sigma)',      Sigma];
                A = 0.5*(A+A');
                ddr = A\[ H*dxidr, H*dxidp;...
                         zeros(n), -eye(n);...
                           eye(n), zeros(n)];
                dxidr = ddr(1:m,1:n);
                dxidp = ddr(1:m,n+1:end);
%                 drdr = Sigma*ddr((end-n+1):end,:);
%                 HH = inv(inv(H) - diag(xi)); HH = 0.5 * (HH + HH');
%                 dxidp = (H\HH)*(Q/(Q'*HH*Q + Sigma));
            end
        end
    end
    function value = l2(dr)
        value = 0.5*log(dr'*invDltSq*dr) - log(delta);
    end
    function value = linf(dr)
        value = log(max(abs(sqrtinvDltSq*dr))) - log(delta);
    end
%     function flag = R1()
%         flag = false;
%         error('incomplete rule');
%     end
%     function flag = Sense()
%         flag = false;
%         error('incomplete rule');
%     end
    %
    [ m,n,p,U,D,Q0,q0,H] = getInitData(h);
    nsubiter = h.method.nsubiter;
    nIter = h.method.niter;
    DltSq = h.method.DeltaSq;
    invDltSq = inv(DltSq);
    sqrtinvDltSq = sqrtm(invDltSq);
    delta = h.rule.delta;
    lambda = h.rule.lambda;
    isfndderiv = h.method.fndderiv;
    %
    switch h.rule.name
        case 'l2'
            events = @l2;
        case 'linf'
            events = @linf;
        case 'Lepskii'
            events = @Lepskii;
        case 'R1'
            events = @R1;
        case 'Sense'
            events = @Sense;
        otherwise
            error('undefined rule');
    end
    %
    tic;
    options = optimset('TolX',1e-9,'Display','notify');
    if isempty(lambda)
        lambda = exp(fzero(@(l)events(fun(exp(l))),log([1e-2 1e4]),options));
%         lambda = fzero(@(l)events(fun(l)),[1e-3 1e5],options);
    end
    [dr,xi,dxidr,dxidp] = fun(lambda);
    %
    h.rule.lambda = lambda;
    h.result.xi = xi;
    h.result.dxidr = dxidr;
    h.result.dxidp = dxidp;
    h.result.r = h.method.r0 + dr;
    h.result.time = toc;
    end