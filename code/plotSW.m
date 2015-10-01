function plotSW(arrh,varargin)
    if nargin>1
        T = varargin{1};
    else
        T = 60;
    end
    
    [ v, y, f] = getrates( arrh, T );
    [ u, spot] = getrates( arrh, []);
    name = [arrh(1).method.name,'; ',arrh(1).rule.name];
    %
    h = gca;
    nextplot = get(h,'NextPlot');
    %
    p1 = plot(v,100*y,'linestyle','-.','linewidth',1);
    hold on;
    p2 = plot(u,100*spot,'Marker','+','linestyle','none');
    p3 = plot(v,100*f,'linestyle','-','linewidth',2);
    if numel(arrh) == 1
        set(p1,'DisplayName',['spot curve ',name]);
        set(p2,'DisplayName',['yield ', name]);
        set(p3,'DisplayName',['forward rate ', name]);
    end
    if nargin>2
        color = varargin{2};
        set(p1,'color',color);
        set(p2,'color',color);
        set(p3,'color',color);
    end
    %
    set(h,'NextPlot',nextplot);
    %legend('location','southeast')
    xlabel 'Term, years'
    ylabel 'Rate, %'
    grid on
end

