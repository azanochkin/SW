function plotSW(arrh,varargin)
    T = 60;
    colorarg = {};
    i = 1;
    while nargin>i
        switch varargin{i}
            case 'time'
                T = varargin{i+1};
            case 'color'
                colorarg = varargin(i:(i+1));
            otherwise
                error('undefined parametr name')
        end
        i = i+2;
    end
    [ v, y, f] = getrates( arrh, T );
    [ u, spot] = getrates( arrh, []);
    %
    h = gca;
    nextplot = get(h,'NextPlot');
    %
    p1 = plot(v,100*y,'linestyle','-.','linewidth',1.5,colorarg{:});
    hold on;
    p2 = plot(u,100*spot,'Marker','+','linestyle','none',colorarg{:});
    p3 = plot(v,100*f,'linestyle','-','linewidth',2,colorarg{:});
    if numel(arrh) == 1
        %name = arrh.data.date;
        if strcmp(arrh(1).method.name,'original')
            name = 'original';
        else
            name = 'regularized';
        end
        %name = ['regul,; ', arrh(1).method.name,'; ',arrh(1).rule.name];
        set(p1,'DisplayName',['spot curve ',name]);
        set(p2,'DisplayName',['traded maturities ', name]);
        set(p3,'DisplayName',['forward rate ', name]);
    end
    %
    set(h,'NextPlot',nextplot);
    %legend('location','southeast')
    xlabel 'Term, years'
    ylabel 'Rate, %'
    grid on
end

