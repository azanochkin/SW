function h = original(h)
    tic;
    %h.method.r0 = h.data.r_mid;
    [m,n,p,U,D,Q0,q0,H,ann_vec,ann_add] = getInitData(h);
    beta = (Q0'*H*Q0)\(p-q0);
    xi = Q0*beta;
    dxidr = -(Q0/(Q0'*H*Q0))*diag(U'*D*(1+H*xi));
    %
    h.result.r = h.method.r0;
    h.result.xi = xi;
    h.result.dxi = dxidr;
    h.result.grad = dxidr'*ann_vec;
    h.result.annuity = ann_add + xi'*ann_vec;
    h.result.time = toc;
end

