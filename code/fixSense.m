function [sen,arrh,res] = fixSense( data,dayGrid,u,mres)
% funtiob find sensitivity level for mean residual == mres in interval
% dayGrid
    function meanres = getmeanres(sen)
        % direct problem
        res = zeros(1,N);
        disp(sen);
        fprintf('\b--> %6.2f%%\n',0 )
        for j = 1:N
            fprintf('\b\b\b\b\b\b\b\b%6.2f%%\n',j/N*100)
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
    sen2 = 0.003;
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

