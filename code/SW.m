function h = SW(data,day,u,varargin)
    h = parser(varargin{:});
    h.data.day_number = day;
    h.data.u = u(:);
    h.data.date = data.Date{day};
    h.data.r_bid = 1e-2*(data.PX_BID(day,u) - data.CRA(day))';
    h.data.r_ask = 1e-2*(data.PX_ASK(day,u) - data.CRA(day))';
    h.data.r_mid = 1e-2*(data.PX_MID(day,u) - data.CRA(day))';
    % determing normalization matrix
    wind = 30;
    switch h.method.normname
        case 'spread'
            normval = max(h.data.r_ask-h.data.r_bid,1e-4);
        case 'spread+'
            normval = max(h.data.r_ask-h.data.r_bid,1e-4); 
            normval(20:end) = 3e-4;
        case 'volatility'
            normval = 1e-2*mean(abs(diff(data.PX_MID(min(1,day-wind):day,u))));
        case 'volatility+'
            normval = 1e-2*std(data.PX_MID((day-wind:day),u));
        case 'avrgspread'
            normval = 1e-2*mean(data.PX_ASK((day-wind:day),u) - ...
                data.PX_BID((day-wind)-1:day-1,u));
        case 'avrgspread+'
            normval = 1e-2*median(data.PX_ASK((day-wind:day),u) - ...
                data.PX_BID((day-wind)-1:day-1,u));
        case 'simple'
            normval = ones(size(u))*3e-4;
        otherwise
            error('undefined')
    end
    h.method.DeltaSq = diag(normval/2)^2;
    %
    switch h.method.r0_type
        case 'mid'
            h.method.r0 = h.data.r_mid;
        case 'mean'
            h.method.r0 = (h.data.r_bid+h.data.r_ask)/2;
        otherwise
            error('undefined r0 data type')
    end
    switch h.method.functional
        case 'original'
            switch h.method.name
                case 'original'
                    h = original(h);
                case 'bounded'
                    h = bounded(h);
                case 'Tikhonov'
                    h.method.niter = 1;
                    h = iterative(h);
                    h.method.name = 'Tikhonov';
                case 'iterative'
                    h = iterative(h);
                case 'explicit'
                    h = explicit(h);
                case 'implicit'
                    h = implicit(h);
                case 'cauchy'
                    h = cauchy(h);
                otherwise
                    error('not a method name')
            end
        case 'new'
            if h.method.nsubiter < 5
                warning('Small number of subiterations');
            end
            switch h.method.name
                case 'original'
                    h = original_new(h);
                case 'Tikhonov'
                    h.method.niter = 1;
                    h = iterative_new(h);
                    h.method.name = 'Tikhonov';
                case 'iterative'
                    h = iterative_new(h);
                otherwise
                    error('not a method name')
            end
            h = getannuit_new(h);
        otherwise
            error('Unknown functional type')
    end
    
end

