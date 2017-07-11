function h = plotSW(arrh,varargin)
    T = 60;
    colorarg = {};
    marker = {};
    MarkerIndices = 1;
    i = 1;
    while nargin>i
        switch varargin{i}
            case 'time'
                T = varargin{i+1};
            case 'color'
                colorarg = varargin(i:(i+1));
            case 'marker'
                marker = varargin(i:(i+1));
            case 'MarkerIndices'
                MarkerIndices = varargin{i+1};
            otherwise
                error('undefined parametr name')
        end
        i = i+2;
    end
    [ v, y, f] = getrates( arrh, T );
    %[ u, spot, fy] = getrates( arrh, []);
    [ u, spot, fy] = getrates( arrh, 1:2:T);
    %
    nextplot = get(gca,'NextPlot');
    %
    h(1) = plot(v,100*y,'linestyle','-.','linewidth',0.5,colorarg{:},...
        marker{:},'MarkerIndices',1:MarkerIndices:length(v));
    hold on;
    h(2) = plot(v,100*f,'linestyle','-','linewidth',1,colorarg{:},...
        marker{:},'MarkerIndices',1:MarkerIndices:length(v));
    %h(3) = plot(u,100*spot,'linestyle','none',marker{:},colorarg{:});
    %h(4) = plot(u,100*fy,'linestyle','none',marker{:},colorarg{:});
    if numel(arrh) == 1
        %name = arrh.data.date;
        if strcmp(arrh(1).method.name,'original')
            name = 'original';
        else
            name = 'regularized';
        end
        %name = ['regul,; ', arrh(1).method.name,'; ',arrh(1).rule.name];
        set(h(1),'DisplayName',['spot curve ',name]);
        %set(h(2),'DisplayName',['traded maturities ', name]);
        set(h(2),'DisplayName',['forward rate ', name]);
    end
    %
    set(gca,'NextPlot',nextplot);
    %legend(,'location','southeast')
    xlabel 'Term, years'
    ylabel 'Rate, %'
    grid on
end

