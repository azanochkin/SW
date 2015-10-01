function S = Wilson_Smooth( a, v, u )
    [U,V] = meshgrid(u,v);
    y1 = a * min(U,V);
    y2 = a * max(U,V);
    S = exp(-y2).*(-cosh(y1).*y1 + sinh(y1).*(y2+1))/2;
end

