function h = original(h)
    tic;
    [m,n,p,U,D,Q0,q0,H] = getInitData(h);
    beta = (Q0'*H*Q0)\(p-q0);
    xi = Q0*beta;
    if h.method.fndderiv
        P = diag(U'*D*(1+H*xi));
        N = D*U*diag(beta);
%         dxidr = -(Q0/(Q0'*H*Q0))*P + ...
%             (eye(m) - Q0*((Q0'*H*Q0)\(Q0'*H)))*N ;
%         dxidp = Q0/(Q0'*H*Q0);
        A = [       H,    -H*Q0 ;...
               -Q0'*H,  zeros(n)];
        A = 0.5*(A+A');
        ddr = A\[H*N,  zeros(m,n); P, -eye(n)];
        dxidr = ddr(1:m,1:n);
        dxidp = ddr(1:m,n+1:end);
    else
        dxidr = zeros(m,n);
        dxidp = zeros(m,n);
    end
    %
    h.result.r = h.method.r0;
    h.result.xi = xi;
    h.result.dxidr = dxidr;
    h.result.dxidp = dxidp;
    h.result.time = toc;
end