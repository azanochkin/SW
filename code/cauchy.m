function h = cauchy(h)
    function dXi = odefun(~,Xi)
        xi = Xi(1:m);
        dxidr = reshape(Xi(m+1:end),[m,n])/denormdxidr;
        %
        tld = U'*D*(1+H*xi);
        dr = (p-q0 - Q0'*H*xi)./tld;
        Q = Q0 + D*U*diag(dr);
        %
        b = (invDltSq*((p - q0 - Q0'*H*xi)./tld))./tld;
        dxi = Q*b;
        %
        invtld = diag(1./tld);
        drdr0 = -invtld*(Q'*H*dxidr);
        tmp = -(eye(n)+(invtld*Q'+diag(DltSq*b)*U'*D)*H*dxidr);
        dbdr0 = invtld*invDltSq*tmp;
        ddxidr = D*U*diag(b)*drdr0+Q*dbdr0;
        dXi = [dxi ; ddxidr(:)*denormdxidr];
    end
    function [value,isterminal,direction] = l2(~,Xi)
        xi = Xi(1:m);
        tld = U'*D*(1+H*xi);
        dr = (p-q0 - Q0'*H*xi)./tld;
        value = sqrt(dr'*invDltSq*dr) - delta;
        isterminal = 1;
        direction = 0;
    end
    function [value,isterminal,direction] = linf(~,Xi)
        xi = Xi(1:m);
        tld = U'*D*(1+H*xi);
        dr = (p-q0 - Q0'*H*xi)./tld;
        value = max(abs(sqrtm(invDltSq)*dr)) - 1;%delta;????
        isterminal = 1;
        direction = 0;
    end
    function [value,isterminal,direction] = R1(s,Xi)
        dXi = odefun(s,Xi);
        dxi = dXi(1:m);
        value = s*(dxi'*H*dxi) - delta^2/(2*exp(1));
        isterminal = s>1e-4;
        direction = 0;
    end
    function [value,isterminal,direction] = Sense(~,Xi)
        xi = Xi(1:m);
        dxidr = reshape(Xi(m+1:end),[m,n])/denormdxidr;
        grad = dxidr'*ann_vec;
        sense = snsfnc(grad,DltSq,ann_add + xi'*ann_vec);
        value = sense - delta;
        isterminal = 1;
        direction = 0;
    end
    %
    [m,n,p,U,D,Q0,q0,H,ann_vec,ann_add] = getInitData(h);
    denormdxidr = 1e-4;
    Xi = zeros(m+m*n,1);
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
        case 'Sense'
            events = @Sense;
        otherwise
            error('undefined rule');
    end
    options = odeset('Events',events,'BDF','on','abstol',1e-4,'reltol',1e-2);
    tic;
    %
    [~,Xi] = ode15s(@odefun,[0 1e-1],Xi,options);
    xi = Xi(end,1:m)';
    dxidr = reshape(Xi(end,(m+1):end),[m,n])/denormdxidr;
    %
    h.result.xi = xi;
    h.result.dxi = dxidr;
    h.result.r = h.method.r0 + (p-q0 - Q0'*H*xi)./(U'*D*(1+H*xi));
    h.result.grad = dxidr'*ann_vec;
    h.result.annuity = ann_add + xi'*ann_vec;
    h.result.sense = snsfnc(h.result.grad,DltSq,h.result.annuity);
    h.method.name = 'cauchy';
    h.result.time = toc;
end

