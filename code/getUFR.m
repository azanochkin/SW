function ufr = getUFR(h)
%GETUFR Give a posteriori estimate of UFR
%   give SW structuge - get UFR estimate 
ufr = h.method.ufr - h.method.ufrvar*h.data.u*h.result.xi;
end

