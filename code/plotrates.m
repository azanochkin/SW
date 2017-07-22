function plotrates(h)
    tenor = h.data.tenor;
    r_min = h.data.r_bid;
    r_max = h.data.r_ask;
    r_mid = h.data.r_mid;
    r_last = h.data.r_mid;
    r0 = h.method.r0;
    r = h.result.r(:);
    subplot(2,2,2)
    hold on
    vec = sqrt(diag(h.method.DeltaSq));
    plot(tenor,1e4*(r_min-r0),'b^-',tenor,1e4*(r_max-r0),'rv-',...
         tenor,1e4*(r_mid-r0), 'g+',tenor,1e4*(r_last-r0),'co',...
         tenor,1e4*(r-r0),'k*');%,u,-1e4*vec,'c--')
    hold on;
    %plot(u,1e4*[vec,-vec],'color',[0.3,0.8,0.3],'Linestyle','--','linewidth',2);
    hold off;
    grid minor
    xlabel 'Term, years'
    ylabel 'Centered quotes, bp'
    legend('location','best')
    legend('bid','ask','mid','last','regularized');%,'volatility bounds')
    %
    dr = r - r0;
    %title(['max resudual = ',num2str(1e4*max(dr)),' bps; relev max residual = ',num2str(max(dr./vec))]);
    title(['max residual = ',num2str(1e4*max(abs(dr))),' bps']);
    subplot(2,2,4)
    title(['Quotes on ',h.data.date])
    plot(tenor,1e2*r_min,'b*',tenor,1e2*r_max,'r*',...
         tenor,1e2*r_mid,'g+',tenor,1e2*r_last,'co',...
         tenor,1e2*r,'k')
    grid on
    legend('bid','ask','mid','last','regularized','location','southeast')
    xlabel('Term, years')
    ylabel('Fixed leg rate, %')
    legend('bid','ask','mid','regulirized')
	legend('location','southeast')
end