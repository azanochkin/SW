function h = original(h)
    tic;
    h.method.r0 = h.data.r_mid;
    [m,n,p,U,D,Q0,q0,H,ann_vec,ann_add] = getInitData(h);
    beta = (Q0'*H*Q0)\(p-q0);
    xi = Q0*beta;
    dxidr = -Q0*inv(Q0'*H*Q0)*diag(U'*D*(1+H*xi));
    %
    h.result.r = h.data.r_mid;
    h.result.xi = xi;
    h.result.dxi = dxidr;
    h.result.grad = dxidr'*ann_vec;
    h.result.annuity = ann_add + xi'*ann_vec;
    snsfnc = sensefnc(h);
    h.result.sense = snsfnc(h.result.grad,h.method.DeltaSq,...
        h.result.annuity);
    h.method.name = 'original';
    h.result.time = toc;
end

