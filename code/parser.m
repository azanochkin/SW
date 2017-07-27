function h = parser( varargin )
    h = getstruct();
    i = 1;
    while nargin>i
        switch varargin{i}
            case 'functional'
                h.method.functional = varargin{i+1};
            case 'method'
                h.method.name = varargin{i+1};
            case 'kernel'
                h.method.kernel = varargin{i+1};
            case 'UFR'
                h.method.ufr = varargin{i+1};
            case 'tau'
                h.method.tau = varargin{i+1};
            case 'tautol'
                h.method.tautol = varargin{i+1};
            case 'convpnt'
                h.method.convpnt = varargin{i+1};
            case 'fixalpha'
                h.method.fixalpha = varargin{i+1};                
            case 'alpha'
                h.method.alpha = varargin{i+1};
            case 'lambda'
                h.rule.lambda = varargin{i+1};
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
            case 'norm'
                h.method.normname = varargin{i+1};
            case 'fndderiv'
                h.method.fndderiv = varargin{i+1};
            case 'profile'
                profMat = varargin{i+1};
                if size(profMat,2)~=2
                    error('profMat will be k \times 2 matrix.')
                end
                h.data.u_ann = profMat(:,1);
                h.data.profile = profMat(:,2);
            otherwise
                error(['unknown field : "',varargin{i},'"'])
        end
        i = i+2;
    end
end

