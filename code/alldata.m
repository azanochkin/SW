N = 10;
arrh(1,N) = getstruct();
for day = 1:N
    if (day/10) == round(day/10)
        disp(day);
    end
    main;
    h.rule.name = 'Sense';
    h.rule.delta = 0.003;
    arrh(1,day) = cauchy(h);
end
%%
time(1,N) = 0;
for day = 1:N
    time(day) = array_h_impl_30_linf(day).result.time;
end