function H = Wilson_Heart( a, v, u, kernel )
    [U,V] = meshgrid(u,v);
    Min = min(U,V);
    Max = max(U,V);
    switch kernel
        case 'SW'
            %Wilson
            H = a * Min - exp(-a*Max).* sinh(a*Min);
        case 'Z'
            %Zanochkin
            H = a * Min - 0.5*(1 - exp(-a*U) - exp(-a*V) + exp(-a*(Max-Min)));
        otherwise
            error('unknow kernel type');
    end
end

