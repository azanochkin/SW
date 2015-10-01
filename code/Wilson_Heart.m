function H = Wilson_Heart( a, v, u )
    [U,V] = meshgrid(u,v);
    Min = min(U,V);
    Max = max(U,V);
    H = a * Min - exp(-a*Max).* sinh(a*Min);
end

