function plotrates(h)
    u = h.data.u;
    r_min = h.data.r_bid;
    r_max = h.data.r_ask;
    r_mid = h.data.r_mid;
    r0 = h.method.r0;
    r = h.result.r(:);
    subplot(1,2,1)
    hold on
    vec = sqrt(diag(h.method.DeltaSq));
    plot(u,1e4*(r_min-r0),'b',u,1e4*(r_max-r0),'r',u,1e4*(r_mid-r0),'g+',u,1e4*(r-r0),'k*');%,u,-1e4*vec,'c--')
    hold on;
    plot(u,1e4*[vec,-vec],'color',[0.3,0.8,0.3],'Linestyle','--','linewidth',2);
    hold off;
    xlabel 'Term, years'
    ylabel 'Centered quotes, bp'
    legend('location','best')
    legend('bid','ask','mid','regularized','volatility bounds')
    %
    dr = r - r0;
    %title(['max resudual = ',num2str(1e4*max(dr)),' bps; relev max residual = ',num2str(max(dr./vec))]);
    title(['max residual = ',num2str(1e4*max(dr)),' bps']);
    subplot(1,2,2)
    title(['Quotes on ',h.data.date])
    hold on
    plot(u,1e2*r_min,'b*',u,1e2*r_max,'r*',u,1e2*r_mid,'g*',u,1e2*r,'k')
    %legend('bid','ask','regularized','location','southeast')
    xlabel('Term, years')
    ylabel('Fixed leg rate, %')
    legend('bid','ask','mid','regulirized')
	legend('location','southeast')
end