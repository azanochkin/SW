function h = original_new(h)
    tic;
    [m,n,p,U,D,Q0,q0,H,ann_vec,ann_add] = getInitData(h);
    nsubiter = h.method.nsubiter;
    %%
    xi = zeros(m,1);
    for j = 1:nsubiter
        Q0_xi = diag(exp(H*xi))*Q0;
        beta = (Q0_xi'*H*Q0_xi)\(p - Q0_xi'*(1 - H*xi));
        fprintf('%i: norm = %6.4f\n',j,norm( xi - Q0_xi*beta));
        xi = Q0_xi*beta;
    end
%    dxidr = -Q0*inv(Q0'*H*Q0)*diag(U'*D*(1+H*xi));
    %
    h.result.r = h.method.r0;
    h.result.xi = xi;
%    h.result.dxi = dxidr;
%     h.result.grad = dxidr'*ann_vec;
%     h.result.annuity = ann_add + xi'*ann_vec;
%     snsfnc = sensefnc(h);
%     h.result.sense = snsfnc(h.result.grad,h.method.DeltaSq,...
%        h.result.annuity);
%     h.method.name = 'original';
    h.result.time = toc;
end

