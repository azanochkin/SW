function h = SW(data,day,u,varargin)
    h = parser(varargin{:});
    h.data.day_number = day;
    h.data.u = u(:);
    h.data.date = data.Date{day};
    h.data.r_bid = 1e-2*(data.PX_BID(day,u) - data.CRA(day))';
    h.data.r_ask = 1e-2*(data.PX_ASK(day,u) - data.CRA(day))';
    h.data.r_mid = 1e-2*(data.PX_MID(day,u) - data.CRA(day))';
    %
    %spread = 1e-2*median(data.PX_ASK((day-30:day),u) - data.PX_BID((day-30)-1:day-1,u));
    %spread = 1e-2*mean(data.PX_ASK((day-30:day),u) - data.PX_BID((day-30)-1:day-1,u));
    %spread = max(h.data.r_ask-h.data.r_bid,1e-4); spread(20:end) = 3e-4;
    %spread = 1e-2*std(data.PX_MID((day-30:day),u) - data.PX_MID((day-30)-1:day-1,u));
    %spread = ones(size(u))*3e-4;
    h.method.DeltaSq = diag(spread/2)^2;
    %
    h.method.r0 = (h.data.r_bid+h.data.r_ask)/2;
    %h.method.r0 = h.data.r_mid;
    switch h.method.name
        case 'original'
            h = original(h);
        case 'bounded'
            h = bounded(h);
        case 'Tikhonov'
            h.method.niter = 1;
            h = iterative(h);
        case 'iterative'
            h = iterative(h);
        case 'explicit'
            h = explicit(h);
        case 'implicit'
            h = implicit(h);
        case 'cauchy'
            h = cauchy(h);
        otherwise
            error('not correct method name')
    end
end

