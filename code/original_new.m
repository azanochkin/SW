function h = original_new(h)
    tic;
    [m,n,p,U,D,Q0,q0,H] = getInitData(h);
    nsubiter = h.method.nsubiter;
    %%
    xi = zeros(m,1);
    for j = 1:nsubiter
        Q0_xi = diag(exp(H*xi))*Q0;
        beta = (Q0_xi'*H*Q0_xi)\(p - Q0_xi'*(1 - H*xi));
        %fprintf('%i: norm = %6.4f\n',j,norm( xi - Q0_xi*beta));
        xi = Q0_xi*beta;
    end
    A = eye(m) + ((diag(1./xi) - H)\H);
    M = diag(U'*D*exp(H*xi));
    N = diag(exp(H*xi))*D*U*diag(beta);
    if n==m
        dxidr = A*Q0_xi*(-(Q0_xi'*H*A*Q0_xi)\M);
    else
        dxidr = A*(N - Q0_xi*((Q0_xi'*H*A*Q0_xi)\(M+Q0_xi'*H*A*N)) );
    end
    %
    ddr = [eye(m) - diag(xi)*H , -Q0_xi; -Q0_xi'*H, zeros(n)]\[N; M];
    %
    h.result.r = h.method.r0;
    h.result.xi = xi;
    h.result.dxi = dxidr;
    h.result.time = toc;
end

