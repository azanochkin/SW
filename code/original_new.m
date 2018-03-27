function h = original_new(h)
    tic;
    [m,n,p,U,D,Q0,q0,H] = getInitData(h);
    nsubiter = h.method.nsubiter;
    %%
    xi = zeros(m,1);
    for j = 1:nsubiter
        eHxi = exp(H*xi);
        Q = diag(eHxi)*Q0;
        QHQ = (Q'*H*Q);
        QHQ = 0.5*(QHQ+QHQ');
        beta = QHQ\(p - Q'*(1 - H*xi));
        xi = Q*beta;
    end
    if h.method.fndderiv
        P = diag(U'*D*eHxi);
        N = diag(eHxi)*D*U*diag(beta);
        A = [H - H*diag(xi)*H , -H*Q ;...
                    -Q'*H ,  zeros(n)];
        A = 0.5*(A+A');
        ddr = A\[H*N,  zeros(m,n); P, -eye(n)];
        dxidr = ddr(1:m,1:n);
        dxidp = ddr(1:m,n+1:end);
    else
        dxidr = zeros(m,n);
        dxidp = zeros(m,n);
    end
    %%
    h.result.r = h.method.r0;
    h.result.xi = xi;
    h.result.dxidr = dxidr;
    h.result.dxidp = dxidp;
    h.result.time = toc;
end

