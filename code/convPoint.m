function h = convPoint( h, method )
    function f = F(h)
        [~,f] = getrates( h, convpnt);
        df = abs(ufr - f);
        f = log(tau/df);
    end
    tau = h.method.tau;
    tautol = h.method.tautol;
    convpnt = h.method.convpnt;
    ufr = h.method.ufr;
    llp = max(h.data.u);
    %
    a = h.method.alpha;
    Fa = F(h);
    h.method.alpha = a - Fa/(convpnt - llp);
    isnconv = true;
    itcnt = 1;
    while isnconv
        h = method(h);
        b = h.method.alpha;
        Fb = F(h);
        isnconv = abs(Fb) > tautol; %abs(exp(Fb) - 1) > tautol
        fprintf('conv iter %2i: alpha = %7.5f\n',itcnt,h.method.alpha);
        if isnconv
            h.method.alpha = b - (b - a) * Fb/(Fb - Fa);
            a = b;
            Fa = Fb;
            itcnt = itcnt + 1;
        end
    end
end

