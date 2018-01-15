function h = implicit_new(h)
    function flag = l2()
        flag = sqrt(dr'*invDltSq*dr)>delta;
    end
    function flag = linf()
        flag = max(abs(sqrtinvDltSq*dr))>delta;
    end
    function flag = R1()
        flag = false;
        error('incomplete rule');
    end
    function flag = analR1()
        flag = false;
        error('incomplete rule');
    end
    function flag = Sense()
        flag = false;
        error('incomplete rule');
    end
    %
    [ m,n,p,U,D,Q0,q0,H] = getInitData( h);
    nsubiter = h.method.nsubiter;
    nMaxIter = h.method.niter;
    DltSq = h.method.DeltaSq;
    invDltSq = inv(DltSq);
    sqrtinvDltSq = sqrtm(invDltSq);
    delta = h.rule.delta;
    %
    lambda = h.rule.lambda;
    if isempty(lambda)
        lambda = 1e4;
    end
    Sigma = lambda*DltSq;
    isfndderiv = h.method.fndderiv;
    %
    switch h.rule.name
        case 'l2'
            events = @l2;
        case 'linf'
            events = @linf;
        case 'R1'
            events = @R1;
        case 'analR1'
            events = @analR1;
        case 'Sense'
            events = @Sense;
        otherwise
            error('undefined rule');
    end
    %
    tic;
    xi = zeros(m,1);
    dxidr = zeros(m,n);
    nIter = 0;
    flag = true;
    while flag
        nIter = nIter+1;
        dxi = zeros(m,1);
        dr = zeros(n,1);
        for j = 1:nsubiter
            eHxi = exp(H*(xi+dxi));
            Q = diag(eHxi)*(Q0 + D*U*diag(dr));
            P = diag(U'*D*eHxi);
            A = (Q'*H*Q+P*Sigma*P);
            A = 0.5*(A+A');
            beta = A\(p - Q'*(1 - H*dxi) + P' * dr);
            dr = Sigma * P * beta;
            dxi = Q*beta;
        end
        xi = xi + dxi;
        if isfndderiv
            N = diag(eHxi)*D*U*diag(beta);
            A = [H - H*diag(dxi)*H,     -H*Q  , -H*N*Sigma;...
                -(H*Q)'           , zeros(n)  ,   -P*Sigma;...
                -(H*N*Sigma)'     ,-(P*Sigma)',      Sigma];
            A = 0.5*(A+A');
            ddr = A\[H*dxidr; zeros(n); eye(n)];
            dxidr = ddr(1:m,:);
            drdr = Sigma*ddr((end-n+1):end,:);
        end
        flag = (nIter < nMaxIter) && events();
    end
    %
    h.method.niter = nIter;
    h.result.xi = xi;
    h.result.dxi = dxidr;
    h.result.r = h.method.r0 + dr;
    h.result.time = toc;
end
