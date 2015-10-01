function [ v, y, f] = getrates( arrh, T)
    N_grid = 300;
    switch numel(T)
        case 0
            v = arrh(1).data.u(1);
        case 1
            v = linspace(0,T,N_grid)';
        case 2
            v = linspace(T(1),T(2),N_grid)';
        otherwise
            v = T(:);
    end
    v(v==0) = 1e-5;
    
    N = numel(arrh);
    y = zeros(length(v),N);
    f = zeros(length(v),N);
    for i=1:N
        u = arrh(i).data.u;
        xi = arrh(i).result.xi;
        w = arrh(i).method.ufr;
        a = arrh(i).method.alpha;

        %%%%%%%%%%%%%%%%%%%%% spot rate
        %H = Wilson_Heart(a,u,u);
        %spot = w - log(1+H*xi)./u;

        %%%%%%%%%%%%%%%%%%%%% yield
        H = Wilson_Heart(a,v,u);
        y(:,i) = w - log(1+H*xi)./v;

        %%%%%%%%%%%%%%%%%%%%% forvard intensity
        G = Wilson_G(a,v,u);
        f(:,i) = w - (G*xi)./(1+H*xi);
    end
end

