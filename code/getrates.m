function [ y, f] = getrates( arrh, v)
    N = numel(arrh);
    y = zeros(length(v),N);
    f = zeros(length(v),N);
    for i=1:N
        u = arrh(i).data.u;
        xi = arrh(i).result.xi;
        w = arrh(i).method.ufr;
        a = arrh(i).method.alpha;
        kernel = arrh(i).method.kernel;
        H = Wilson_Heart(a,v,u,kernel);
        G = Wilson_G(a,v,u,kernel);

        switch arrh(i).method.functional
            case 'original'
                f(:,i) = w - (G*xi)./(1+H*xi);
                y(:,i) = w - log(1+H*xi)./v;
            case 'new'
                f(:,i) = w - G*xi;
                y(:,i) = w - (H*xi)./v;
        end
        indNaN = v==0;
        y(indNaN,i) = f(indNaN,i);
    end
end

