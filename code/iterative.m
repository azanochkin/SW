function h = iterative(h)
    function flag = l2()
        flag = sqrt(dr'*invDltSq*dr)>delta;
        %flag = (dr'*invDltSq*dr_prev)>delta*sqrt(dr_prev'*invDltSq*dr_prev);
    end
    function flag = linf()
        flag = max(abs(sqrtm(invDltSq)*dr))>delta;
    end
    function flag = R1()
        tmp = tld\(invDltSq*dr);
        tn = 2*nIter;
        flag = (tmp'*Q'*H*Q*tmp)/lambda>delta^2*(1/tn)*(1-1/tn)^(tn-1);
    end
%     function flag = Sense()
%         grad = dxidr'*ann_vec;
%         flag = snsfnc(grad,DltSq,ann_add + xi'*ann_vec)< delta;
%     end
%     function flag = Lepskii()
%         tau = 1+1/nIter;
%         gamma_p = sqrt(0.5/nIter)*(1-0.5/nIter)^(nIter-0.5);
%         c = nIter * gamma_p;
%         flag = dxi'*H*dxi < c^2*(((1-q)^2)/q)*delta^2/lambda;
%     end
    %
    [ m,n,p,U,D,Q0,q0,H] = getInitData(h);
    nsubiter = h.method.nsubiter;
    nIter = h.method.niter;
    DltSq = h.method.DeltaSq;
    invDltSq = inv(DltSq);
    delta = h.rule.delta;
    %snsfnc = sensefnc(h); 
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
    %dr_prev = zeros(n,1);
    lambda = 1e5;
    q = 0.99;
    flag = true;
    tic;
    while flag
        xi = zeros(m,1);
        dxidr = zeros(m,n);
        dxidp = zeros(m,n);
        Q = Q0;
        tld = diag(U'*D*(1+H*Q0*(Q0'*H*Q0\(p-q0))));
        %
        for i = 1:nIter
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
            drdr = (eye(n) - lambda*DltSq*M'*inv(Q'*H*Q+lambda*M*DltSq*M')*M)*(eye(n)+lambda*DltSq*N'*H*dxidr);
            dxidr = N*drdr;
        end
        flag = events();
        lambda = lambda*q;
        %dr_prev = dr;
    end
    %
    h.result.xi = xi;
    h.result.dxidr = dxidr;
    h.result.dxidp = dxidp;
    h.result.r = h.method.r0 + dr;
    h.result.time = toc;
end

