function h = implicit(h)
    function flag = l2()
        flag = sqrt(dr'*invDltSq*dr)>delta;
    end
    function flag = linf()
        flag = max(abs(sqrtm(invDltSq)*dr))>delta;
    end
    function flag = R1()
        tmp = tld\(invDltSq*dr);
        flag = nIter*(tmp'*Q'*H*Q*tmp)>delta^2*lambda/(2*exp(1));
    end
    function flag = analR1()
        flag = nIter*(dxi'*H*dxi)>delta^2/(2*lambda*exp(1));
    end
    function flag = Sense()
        grad = dxidr'*ann_vec;
        flag = snsfnc(grad,DltSq,ann_add + xi'*ann_vec)< delta;
    end
    %
    [ m,n,p,U,D,Q0,q0,H,ann_vec,ann_add] = getInitData( h);
    lambda = 1e5;
    nsubiter = h.method.nsubiter;
    DltSq = h.method.DeltaSq;
    invDltSq = inv(DltSq);
    delta = h.rule.delta;
    snsfnc = sensefnc(h.method.sensefncname); 
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
    xi = zeros(m,1);
    dxidr = zeros(m,n);
    Q = Q0;
    %tld = diag(U'*D*(1+H*xi));
    tld = diag(U'*D*(1+H*Q0*(Q0'*H*Q0\(p-q0))));
    %
    tic;
    flag = true;
    nIter = 1;
    %
    while flag
        for j = 1:nsubiter
            beta = (Q0'*H*Q+lambda*tld*DltSq*tld)\(p - q0 - Q0'*H*xi);
            dr = lambda * DltSq * tld * beta;
            dxi = Q*beta;
            Q = Q0 + D*U*diag(dr);
            tld = diag(U'*D*(1+H*(xi+dxi)));
        end
        xi = xi + dxi;       
        Du = D*(1+H*xi);
        M = diag(U'*Du)+Q'*H*D*U*diag(beta);
        N = -(Q/(Q'*H*Q))*diag(U'*Du);
        dxidr = N*(eye(n) - lambda*DltSq*M'*inv(Q'*H*Q+lambda*M*DltSq*M')*M)*(eye(n)+lambda*DltSq*N'*H*dxidr);
        %    
        flag = events();
        nIter = nIter+1;
    end
    h.method.niter = nIter-1;
    h.result.xi = xi;
    h.result.dxi = dxidr;
    h.result.r = h.method.r0 + dr;
    h.result.grad = dxidr'*ann_vec;
    h.result.annuity = ann_add + xi'*ann_vec;
    h.result.sense = snsfnc(h.result.grad,DltSq,h.result.annuity);
    h.method.name = 'implicit';
    h.result.time = toc;
end
