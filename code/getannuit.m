function [ h ] = getannuit( h )
%GETANNUIT_NEW Compute annuity price and delta's
v = h.data.u_ann;
pr = h.data.profile;
u = h.data.u;
w = h.method.ufr;
a = h.method.alpha;
kernel = h.method.kernel;
H = Wilson_Heart(a,v,u,kernel);
%
xi = h.result.xi;
dxi = h.result.dxidr;
%
switch h.method.functional
    case 'original'
        d = exp(-w*v);
        h.result.annuity = (1 + H*xi)'*(d.*pr);
    case 'new'
        d = exp(-w*v + H*xi);
        h.result.annuity = d'*pr;
end
h.result.grad = (H*dxi)'*(d.*pr);
end

