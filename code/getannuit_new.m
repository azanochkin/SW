function [ h ] = getannuit_new( h )
%GETANNUIT_NEW Compute annuity price and delta's
v = h.data.u_ann;
pr = h.data.profile;
u = h.data.u;
w = h.method.ufr;
a = h.method.alpha;
H = Wilson_Heart(a,v,u);
%
xi = h.result.xi;
dxi = h.result.dxi;
%
d = exp(-w*v + H*xi);
h.result.annuity = d'*pr;
h.result.grad = (H*dxi)'*(d.*pr);
%
snsfnc = sensefnc(h);
h.result.sense = snsfnc(h.result.grad,h.method.DeltaSq,...
    h.result.annuity);
end

