function G = Wilson_G( a,v,u )
    [U,V] = meshgrid(u,v);
    G = U;
    ind = U>V;
    T = a-a*exp(-a*U).*cosh(a*V);
    G(ind) = T(ind);
    T = a*exp(-a*V).*sinh(a*U);
    G(~ind) = T(~ind);
    
end

