function h = explicit(h)
    function flag = l2()
        flag = sqrt(dr'*invDltSq*dr)>delta;
    end
    function flag = linf()
        flag = max(abs(sqrtm(invDltSq)*dr))>delta;
    end
    function flag = R1()
        tmp = tld\(invDltSq*dr);
        flag = nIter*(tmp'*Q'*H*Q*tmp)>delta^2/(2*exp(1)*mu);
    end
    function flag = analR1()
        flag = nIter*(dxi'*H*dxi)>delta^2*mu/(2*exp(1));
    end
%     function flag = Sense()
%         grad = dxidr'*ann_vec;
%         flag = snsfnc(grad,DltSq,ann_add + xi'*ann_vec)< delta;
%         error('not finished implementation');
%     end
    %
    [ m,n,p,U,D,Q0,q0,H] = getInitData( h);
    nsubiter = h.method.nsubiter;
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
    dxidp = zeros(m,n);
    Q = Q0;
    invtld = diag(1./(U'*D*(1+H*Q0*(Q0'*H*Q0\(p-q0)))));
    mu = 9.9e-1/norm(invtld*invDltSq*invtld*Q'*H*Q);
    %
    tic;
    flag = true;
    nIter = 1;
    while flag
        for j = 1:nsubiter
            beta = mu*(invtld*invDltSq*invtld*(p - q0 - Q0'*H*xi));
            dxi = Q*beta;
            dr = invtld*(p-q0 - Q0'*H*(xi+dxi));
            Q = Q0 + D*U*diag(dr);
            invtld = diag(1./(U'*D*(1+H*(xi+dxi))));
        end
        xi = xi + dxi;        
        %    
        flag = events();
        nIter = nIter+1;
    end
    %
    h.method.niter = nIter-1;
    h.result.xi = xi;
    h.result.dxidr = dxidr;
    h.result.dxidp = dxidp;
    h.result.r = h.method.r0 + dr;
    h.result.time = toc;
end

