function h = convPoint( h, method )
    tau = h.method.tau;
    tautol = h.method.tautol;
    convpnt = h.method.convpnt;
    ufr = h.method.ufr;
    llp = max(h.data.u);
    itcnt = 0;  
    isnConv = true;
    while isnConv
        [~,f] = getrates( h, convpnt);
        df = abs(ufr - f);
        isnConv = abs(tau/df - 1) > tautol;
        if isnConv
            itcnt = itcnt + 1;
            h.method.alpha = h.method.alpha - log(tau/df)/(convpnt - llp);
            h = method(h);
            fprintf('conv iter %2i: alpha = %7.5f\n',itcnt,h.method.alpha);
        end
    end
end

