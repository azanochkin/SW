function [ m,n,p,U,D,Q0,q0,H,ann_vec,ann_add] = getInitData(h)
    a = h.method.alpha;
    w = h.method.ufr;
    kernel = h.method.kernel;
    u = h.data.u;
    r0 = h.method.r0;
    m = length(u);
    n = length(r0);

    p = ones(n,1);
    U = reshape((repmat(triu(ones(n),0),1,m/n))',n,m)'; 
    E = (n/m)*reshape(([repmat(zeros(n),1,m/n-1),eye(n)])',n,m)';
    C0 = E + U * diag(r0);

    d = exp(-w*u);
    D = diag(d);
    Q0 = D*C0;
    q0 = C0'*d;

    H = Wilson_Heart(a,u,u,kernel);
    
    v = h.data.u_ann;
    pr = h.data.profile;
    ann_vec = Wilson_Heart(a,u,v,kernel)*(exp(-w*v).*pr);
    ann_add = exp(-w*v)'*pr;
end

