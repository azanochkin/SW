function [u] = plotrates(h)
    u = h.data.u;
    r_min = h.data.r_bid;
    r_max = h.data.r_ask;
    r_mid = h.data.r_mid;
    r0 = h.method.r0;
    r = h.result.r(:);
    
    subplot(2,2,2)
    hold on
    vec = sqrt(diag(h.method.DeltaSq));
    plot(u,r_min-r0,'b',u,r_max-r0,'r',u,r_mid-r0,'g+',u,r-r0,'k*',u,vec,'c--',u,-vec,'c--')
    %title(['clear resudual = ',num2str(1e4*norm(r-r0)),' bps; relev residual = ',num2str(delta)]);
    subplot(2,2,4)
    hold on
    plot(u,r_min,'b*',u,r_max,'r*',u,r_mid,'g*',u,r,'k')
    %legend('bid','ask','regularized','location','southeast')
    xlabel('Term, years')
    ylabel('Fixed leg rate, %')
    title (['Quotes of 1 Y-fixed vs 6 M-float EUR IRS on ',h.data.date])
end