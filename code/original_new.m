function h = original_new(h)
    tic;
    [m,n,p,U,D,Q0,q0,H] = getInitData(h);
    nsubiter = h.method.nsubiter;
    %%
    xi = zeros(m,1);
    for j = 1:nsubiter
        Q0_xi = diag(exp(H*xi))*Q0;
        beta = (Q0_xi'*H*Q0_xi)\(p - Q0_xi'*(1 - H*xi));
        xi = Q0_xi*beta;
    end
    eHxi = exp(H*xi);
    Q0_xi = diag(eHxi)*Q0;
    M = diag(U'*D*eHxi);
    N = diag(eHxi)*D*U*diag(beta);
    A = [H - H*diag(xi)*H , -H*Q0_xi ;...
                -Q0_xi'*H ,  zeros(n)];
    A = 0.5*(A+A');
    ddr = A\[H*N; M];
    dxidr = ddr(1:m,:);
    %
    h.result.r = h.method.r0;
    h.result.xi = xi;
    h.result.dxi = dxidr;
    h.result.time = toc;
end

