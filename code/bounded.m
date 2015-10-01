function h = bounded(h)
    [ m,n,p,U,D,Q0,q0,H,ann_vec,ann_add ] = getInitData(h);
    nsubiter = h.method.nsubiter;
    DltSq = h.method.DeltaSq;
    delta = h.rule.delta;
    Q = Q0;
    tld = diag(U'*D*(1+H*Q0*(Q0'*H*Q0\(p-q0))));
    Dlt = sqrtm(DltSq);
    options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
    %
    tic;
    for j = 1:nsubiter
        QHQ = Q'*H*Q;
        QHQ = (QHQ'+QHQ)/2;
        dltVec = delta*tld*Dlt*ones(n,1);
        b = quadprog(QHQ,[],[Q0'*H*Q;-Q0'*H*Q],[dltVec+(p-q0),dltVec-(p-q0)],[],[],[],[],[],options);
        xi = Q*b;
        dr = tld\(p-q0 - Q0'*H*xi);
        Q = Q0 + D*U*diag(dr);
        tld = diag(U'*D*(1+H*xi));
    end
    h.result.xi = xi;
    h.result.r = h.method.r0 + dr;
    h.result.annuity = ann_add + xi'*ann_vec;
    h.method.name = 'bounded';
    h.result.time = toc;
end

