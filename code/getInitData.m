function [ m,n,p,U,D,Q0,q0,H,ann_vec,ann_add] = getInitData(h)
    a = h.method.alpha;
    w = h.method.ufr;
    kernel = h.method.kernel;
    u = h.data.u(:);
    tenor = h.data.tenor(:);
    r0 = h.method.r0;
    m = length(u);
    n = length(r0);

    p = ones(n,1);
    U = zeros(m,n);
    U(u<=tenor') = 1;
    E = zeros(m,n);
    E(u==tenor') = 1;
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

