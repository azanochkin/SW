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
    N_grid = 300;
    switch numel(T)
        case 1
            v = linspace(0,T,N_grid)';
        case 2
            v = linspace(T(1),T(2),N_grid)';
        otherwise
            v = T(:);
    end
    [ y, f] = getrates( arrh,v );
    %
    nextplot = get(gca,'NextPlot');
    h(:,1) = plot(v,100*y,'linestyle','-.','linewidth',0.5,colorarg{:},...
        marker{:},'MarkerIndices',1:MarkerIndices:length(v));
    hold on;
    h(:,2) = plot(v,100*f,'linestyle','-','linewidth',1,colorarg{:},...
        marker{:},'MarkerIndices',1:MarkerIndices:length(v));
    for i = 1:numel(arrh)
        h(i,2).Color = h(i,1).Color;
        func = arrh(i).method.functional(1);
        name = arrh(i).method.name(1:4);
        tnr = arrh(i).data.tenor(end);
        mask = arrh(i).data.mask;
        if any(~mask(1:end-1) & mask(2:end))
            ful = '?';
        else
            ful = ':';
        end
        fulname = sprintf('spot curve %s_%c 1%c%i',name,func,ful,tnr);
        set(h(i,1),'DisplayName',fulname);
        fulname = sprintf('forward rate %s_%c 1%c%i',name,func,ful,tnr);
        set(h(i,2),'DisplayName',fulname);
    end
    %
%     yl = ylim();
%     if yl(2)>30
%         yl(2) = 100*max(max([y(v<30,:) f(v<30,:)]));
%     end
%     ylim(yl);
    set(gca,'NextPlot',nextplot);
    %legend(,'location','southeast')
    xlabel 'Term, years'
    ylabel 'Rate, %'
    grid on
end

