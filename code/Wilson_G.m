function G = Wilson_G( a,v,u,kernel,ufrvar )
    [U,V] = meshgrid(u,v);
    if nargin < 4
        kernel =  'SW';
    end
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
    if nargin < 5
        ufrvar = 0;
    end
    G = ufrvar*U;
    ind = U>V;
    G(ind) = G(ind) + T1(ind); 
    G(~ind) = G(~ind) + T2(~ind);
end

