function H = Wilson_Heart( a, v, u, kernel, ufrvar )
    [U,V] = meshgrid(u,v);
    Min = min(U,V);
    Max = max(U,V);
    if nargin < 4
        kernel =  'SW';
    end
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
    if nargin < 5
        ufrvar = 0;
    end
    H = H + ufrvar * Min .* Max;
end

