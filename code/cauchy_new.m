function h = cauchy_new(h)
    function dXi = odefun(~,Xi)
        xi = Xi(1:m);
        %
        eHxi = exp(H*xi);
        tld = U'*D*eHxi;
        dr = (p - Q0'*eHxi)./tld;
        Q_xi = diag(eHxi)*(Q0 + D*U*diag(dr));
        %
        b = (invDltSq*dr)./tld;
        dXi = Q_xi*b;
        %
        if isfndderiv
            dxidr = reshape(Xi(m+1:end),[m,n])/denormdxidr;
            iTld = diag(1./tld);
            N = diag(eHxi)*D*U*diag(b);
            %
            A = diag(dXi) - N*iTld*Q_xi' - Q_xi*iTld*N' -...
                Q_xi * iTld * invDltSq * iTld * Q_xi';
            A = 0.5*(A+A');
            ddxidr = A*H*dxidr - Q_xi*iTld*invDltSq;
            dXi = [dXi ; ddxidr(:)*denormdxidr];
        end
    end
    function [value,isterminal,direction] = l2(~,Xi)
        xi = Xi(1:m);
        eHxi = exp(H*xi);
        tld = U'*D*eHxi;
        dr = (p - Q0'*eHxi)./tld;
        value = sqrt(dr'*invDltSq*dr) - delta;
        isterminal = 1;
        direction = 0;
    end
    function [value,isterminal,direction] = linf(~,Xi)
        xi = Xi(1:m);
        eHxi = exp(H*xi);
        tld = U'*D*eHxi;
        dr = (p - Q0'*eHxi)./tld;
        value = max(abs(sqrtinvDltSq*dr)) - delta;
        isterminal = 1;
        direction = 0;
    end
%     function [value,isterminal,direction] = R1(s,Xi)
%         dXi = odefun(s,Xi);
%         dxi = dXi(1:m);
%         value = s*(dxi'*H*dxi) - delta^2/(2*exp(1));
%         isterminal = s>1e-4;
%         direction = 0;
%     end
%     function [value,isterminal,direction] = Sense(~,Xi)
%         xi = Xi(1:m);
%         dxidr = reshape(Xi(m+1:end),[m,n])/denormdxidr;
%         grad = dxidr'*ann_vec;
%         sense = snsfnc(grad,DltSq,ann_add + xi'*ann_vec);
%         value = sense - delta;
%         isterminal = 1;
%         direction = 0;
%     end
    %
    [m,n,p,U,D,Q0,q0,H] = getInitData(h);
    denormdxidr = 1e-3;
    DltSq = h.method.DeltaSq;
    invDltSq = inv(DltSq);
    sqrtinvDltSq = sqrtm(invDltSq);
    delta = h.rule.delta;
    isfndderiv = h.method.fndderiv;
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
    tic;
    Xi = zeros(m,1);
    if isfndderiv
        Xi = [Xi; zeros(m*n,1)];
    end
    options = odeset('BDF','on','abstol',1e-6,'reltol',1e-4);
    if isempty(h.rule.lambda)
        options = odeset(options,'Events',events);
        [h.rule.lambda,Xi,~,~,ie] = ode15s(@odefun,[0 2e1],Xi,options);
        if isempty(ie)
            warning(['Cauchy_new : condition is not fulfilled. ',...
                     'Change reg. parametr interval'])
        end
    else
        [~,Xi] = ode15s(@odefun,[0 h.rule.lambda],Xi,options);
    end
    xi = Xi(end,1:m)';
    if isfndderiv
        dxidr = reshape(Xi(end,(m+1):end),[m,n])/denormdxidr;
    else
        dxidr = zeros(m,n);
    end
    %
    h.result.xi = xi;
    h.result.dxidr = dxidr;
    h.result.dxidp = zeros(m,n);
    %
    eHxi = exp(H*xi);
    tld = U'*D*eHxi;
    dr = (p - Q0'*eHxi)./tld;
    h.result.r = h.method.r0 + dr;
    h.result.time = toc;
end

