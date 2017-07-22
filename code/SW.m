function h = SW(data,date,varargin)
    h = parser(varargin{:});
    day=find(ismember(data.Date,date));
    mask = data.liquid_mask;
    h.data.day_number = day;
    h.data.tenor = data.tenor(mask);
    h.data.u = 1:max(h.data.tenor);
    h.data.date = data.Date{day};
    h.data.r_bid = 1e-2*(data.PX_BID(day,mask) - data.CRA(day))';
    h.data.r_ask = 1e-2*(data.PX_ASK(day,mask) - data.CRA(day))';
    h.data.r_mid = 1e-2*(data.PX_MID(day,mask) - data.CRA(day))';
    h.data.r_last = 1e-2*(data.PX_LAST(day,mask) - data.CRA(day))';
    % determing normalization matrix
    wind = 30;
    
    switch h.method.normname
        case 'spread'
            normval = max(h.data.r_ask-h.data.r_bid,1e-4);
        case 'spread+'
            normval = max(h.data.r_ask-h.data.r_bid,1e-4); 
            normval(20:end) = 3e-4;
        case 'volatility'
            normval = 1e-2*mean(abs(diff(data.PX_LAST(min(1,day-wind):day,mask))));
        case 'volatility+'
            normval = 1e-2*std(data.PX_LAST((day-wind:day),mask));
        case 'avrgspread'
            normval = 1e-2*mean(data.PX_ASK((day-wind:day),mask) - ...
                data.PX_BID((day-wind)-1:day-1,mask));
        case 'avrgspread+'
            normval = 1e-2*median(data.PX_ASK((day-wind:day),mask) - ...
                data.PX_BID((day-wind)-1:day-1,mask));
        case 'simple'
            normval = ones(size(h.data.r_ask))*3e-4;
        otherwise
            error('undefined')
    end
    h.method.DeltaSq = diag(normval/2)^2;
    %
    switch h.method.r0_type
        case 'last'
            h.method.r0 = h.data.r_last;
        case 'mid'
            h.method.r0 = h.data.r_mid;
        case 'mean'
            h.method.r0 = (h.data.r_bid+h.data.r_ask)/2;
        otherwise
            error('undefined r0 data type')
    end
    switch h.method.functional
        case 'original'
            getannuit = @(h)h;
            switch h.method.name
                case 'original'
                    method = @original;
                case 'bounded'
                    method = @bounded;
                case 'Tikhonov'
                    h.method.niter = 1;
                    method = @iterative;
                case 'iterative'
                    method = @iterative;
                case 'explicit'
                    method = @explicit;
                case 'implicit'
                    method = @implicit;
                case 'cauchy'
                    method = @cauchy;
                otherwise
                    error('not a method name')
            end
        case 'new'
            if h.method.nsubiter < 3
                warning('Small number of subiterations');
            end
            getannuit = @getannuit_new;
            switch h.method.name
                case 'original'
                    method = @original_new;
                case 'Tikhonov'
                    h.method.niter = 1;
                    method = @iterative_new;
                case 'iterative'
                    method = @iterative_new;
                case 'implicit'
                    method = @implicit_new;
                case 'cauchy'
                    method = @cauchy_new;
                otherwise
                    error('not a method name')
            end
        otherwise
            error('Unknown functional type')
    end
    h = method(h);
    if ~isempty(h.method.tau) && ~isempty(h.method.convpnt)
        h = convPoint( h, method );
    end
    h = getannuit(h);
end

