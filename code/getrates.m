function [ v, y, f] = getrates( arrh, T)
    N_grid = 300;
    switch numel(T)
        case 0
            v = arrh(1).data.u;
        case 1
            v = linspace(0,T,N_grid)';
        case 2
            v = linspace(T(1),T(2),N_grid)';
        otherwise
            v = T(:);
    end
    N = numel(arrh);
    y = zeros(length(v),N);
    f = zeros(length(v),N);
    for i=1:N
        u = arrh(i).data.u;
        xi = arrh(i).result.xi;
        w = arrh(i).method.ufr;
        a = arrh(i).method.alpha;
        H = Wilson_Heart(a,v,u);
        G = Wilson_G(a,v,u);

        switch arrh(i).method.functional
            case 'original'
                %%%%%%%%%%%%%%%%%%%%% forvard intensity
                f(:,i) = w - (G*xi)./(1+H*xi);
                %%%%%%%%%%%%%%%%%%%%% yield
                y(:,i) = w - log(1+H*xi)./v;        
                %%%%%%%%%%%%%%%%%%%%% spot rate
                %H = Wilson_Heart(a,u,u);
                %spot = w - log(1+H*xi)./u;
            case 'new'
                f(:,i) = w - G*xi;
                y(:,i) = w - (H*xi)./v;
        end
        indNaN = v==0;
        y(indNaN,i) = f(indNaN,i);
    end
end

