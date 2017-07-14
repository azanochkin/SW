function G = Wilson_G( a,v,u,kernel )
    [U,V] = meshgrid(u,v);
    G = U;
    switch kernel
        case 'SW'
            T1 = a-a*exp(-a*U).*cosh(a*V);
            T2 = a*exp(-a*V).*sinh(a*U);
        case 'Z'
            T1 = a - 0.5 * a * (exp(-a*V)+exp(-a*(U-V)));
            T2 = - 0.5 * a * (exp(-a*V)-exp(-a*(V-U)));
        otherwise
            error('unknow kernel type');
    end
    ind = U>V;
    G(ind) = T1(ind); 
    G(~ind) = T2(~ind); 
end

