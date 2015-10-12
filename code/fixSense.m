function [sen,arrh,res] = fixSense( data,dayGrid,u,mres)
    function meanres = getmeanres(sen)
        res = zeros(1,N);
        disp(sen);
        for j = 1:N
            day = dayGrid(j);
            h1 = SW(data,day,u,'method','implicit','rule','Sense','delta',sen);
            arrh(j) = h1;
            dr = h1.result.r - h1.method.r0;
            %res(j) = (sqrt(dr'*(h1.method.DeltaSq\dr)))/sqrt(numel(u));
            %res(j) = mean(abs(sqrtm(h1.method.DeltaSq)\dr));
            res(j) = max(abs(sqrtm(h1.method.DeltaSq)\dr));
        end
        meanres = mean(res);
    end
    N = numel(dayGrid);
    arrh(N) = getstruct();
    sen1 = 0.002;
    sen2 = 0.003767;
    mres1 = getmeanres(sen1);
    mres2 = getmeanres(sen2);
    while abs(mres1-mres)>1e-2
        sen = ((log(mres2) - log(mres))*sen1 - (log(mres1) - log(mres))*sen2)/(log(mres2)-log(mres1));
        sen2 = sen1;
        sen1 = sen;
        mres2 = mres1;
        mres1 = getmeanres(sen);
    end 
end

