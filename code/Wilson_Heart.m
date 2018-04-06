function H = Wilson_Heart( a, v, u, kernel )
    [U,V] = meshgrid(u,v);
    Min = min(U,V);
    Max = max(U,V);
    if nargin < 4
        kernel =  'SW';
    end
    if isnumeric(kernel)
        H = kernel * Min .* Max;
        kernel =  'SW';
    else
        H = zeros(size(U));
    end
    switch kernel
        case 'SW'
            %Wilson
            H = H + a * Min - exp(-a*Max).* sinh(a*Min);
        case 'Z'
            %Zanochkin
            H = H + a * Min - 0.5*(1 - exp(-a*U) - exp(-a*V) + exp(-a*(Max-Min)));
        otherwise
            error('unknow kernel type');
    end
end

