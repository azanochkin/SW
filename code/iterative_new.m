function h = iterative_new(h)
    function flag = l2()
        flag = sqrt(dr'*invDltSq*dr)>delta;
        %flag = (dr'*invDltSq*dr_prev)>delta*sqrt(dr_prev'*invDltSq*dr_prev);
    end
    function flag = linf()
        flag = max(abs(sqrtm(invDltSq)*dr))>delta;
    end
    function flag = R1()
        flag = false;
        error('incomplete rule');
    end
    function flag = Sense()
        flag = false;
        error('incomplete rule');
    end
    %
    [ m,n,p,U,D,Q0,q0,H] = getInitData(h);
    nsubiter = h.method.nsubiter;
    nIter = h.method.niter;
    DltSq = h.method.DeltaSq;
    invDltSq = inv(DltSq);
    delta = h.rule.delta;
    snsfnc = sensefnc(h);
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
    lambda = 1e5;
    q = 0.99;
    flag = true;
    tic;
    while flag
        xi = zeros(m,1);
        dxj = zeros(m,1);
        dxidr = zeros(m,n);
        Q = Q0;
        %
        for i = 1:nIter
            for j = 1:nsubiter
                eHxj = exp(H*(xi+dxj));
                Q0_xj = diag(eHxj)*Q0;
                Q_xj = diag(eHxj)*Q;
                tld = diag(U'*D*eHxj);
                beta = (Q0_xj'*H*Q_xj+lambda*tld*DltSq*tld)\...
                    (p - Q0_xj'*(1-H*dxj));
                dr = lambda * DltSq * tld * beta;
                %fprintf('%i > %i: norm = %6.4f\n',i,j,norm( dxj - Q_xj*beta));
                dxj = Q_xj*beta;
                Q = Q0 + D*U*diag(dr);          
            end
            xi = xi + dxj;
            eHxi = exp(H*xi);
            Q_xi = diag(eHxi)*Q;
            M = diag(U'*D*eHxi);
            N = diag(eHxi)*D*U*diag(beta);
            A = [ H - H*diag(dxj)*H,  -H*Q_xi, -H*N;...
                 -(H*Q_xi)'       , zeros(n),   -M;...
                 -(H*N)'          ,      -M', invDltSq/lambda];
            A = 0.5*(A+A');
            ddr = A\[H*dxidr; zeros(n); invDltSq/lambda];
            dxidr = ddr(1:m,:);
            drdr = ddr((end-n+1):end,:);
        end
        flag = events();
        lambda = lambda*q;
    end
    h.result.xi = xi;
    h.result.dxi = dxidr;
    h.result.r = h.method.r0 + dr;
    h.result.time = toc;
    end