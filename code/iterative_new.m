function h = iterative_new(h)
    function [dr,xi,dxidr] = fun(lambda)
        Sigma = lambda*DltSq;
        xi = zeros(m,1);
        dxidr = zeros(m,n);
        for i = 1:nIter
            dxi = zeros(m,1);
            dr = zeros(n,1);
            for j = 1:nsubiter
                eHxi = exp(H*(xi+dxi));
                Q = diag(eHxi)*(Q0 + D*U*diag(dr));
                P = diag(U'*D*eHxi);
                A = (Q'*H*Q+P*Sigma*P);
                A = 0.5*(A+A');
                beta = A\(p - Q'*(1 - H*dxi) + P' * dr);
%                 fprintf('%2i(%i), rcond = %4.2e : dr = %4.2f, dxj = %4.2f\n',...
%                     j,ishermitian(A),rcond(Q'*H*Q),...
%                     1e4*norm(dr - Sigma*P*beta), norm(dxi - Q*beta));
                dr = Sigma * P * beta;
                dxi = Q*beta;
            end
            xi = xi + dxi;
            %
%             A = 0.5*(inv(H) - diag(dxi) + Q*inv(P*Sigma*P)*Q') + Q*inv(P)*diag(beta)*U'*D*eHxi;
%             R = chol(A+A');
            %
            if isfndderiv
                N = diag(eHxi)*D*U*diag(beta);
                A = [ H - H*diag(dxi)*H,     -H*Q  , -H*N*Sigma;...
                     -(H*Q)'           , zeros(n)  ,   -P*Sigma;...
                     -(H*N*Sigma)'     ,-(P*Sigma)',      Sigma];
                A = 0.5*(A+A');
                ddr = A\[H*dxidr; zeros(n); eye(n)];
                dxidr = ddr(1:m,:);
                drdr = Sigma*ddr((end-n+1):end,:);
            end
        end
    end
    function value = l2(dr)
        value = sqrt(dr'*invDltSq*dr) - delta;
        %flag = (dr'*invDltSq*dr_prev)>delta*sqrt(dr_prev'*invDltSq*dr_prev);
    end
    function value = linf(dr)
        value = max(abs(sqrtinvDltSq*dr)) - delta;
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
        lambda = fzero(@(l)events(fun(l)),[1e-1 1e5],options);
    end
    [dr,xi,dxidr] = fun(lambda);
    h.rule.lambda = lambda;
    h.result.xi = xi;
    h.result.dxi = dxidr;
    h.result.r = h.method.r0 + dr;
    h.result.time = toc;
    end