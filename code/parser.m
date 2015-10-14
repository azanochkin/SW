function h = parser( varargin )
    h = getstruct();
    i = 1;
    while nargin>i
        switch varargin{i}
            case 'method'
                h.method.name = varargin{i+1};
            case 'UFR'
                h.method.ufr = varargin{i+1};
            case 'alpha'
                h.method.alpha = varargin{i+1};
            case 'delta'
                h.rule.delta = varargin{i+1};
            case 'nsubiter'
                h.method.nsubiter = varargin{i+1};
            case 'sensefnc'
                h.method.sensefncname = varargin{i+1};
            case 'denormsense'
                h.method.denormsense = varargin{i+1};
            case 'niter'
                h.method.niter = varargin{i+1};
            case 'rule'
                h.rule.name = varargin{i+1};
            otherwise
                error(['unknown field : "',varargin{i},'"'])
        end
        i = i+2;
    end
end

