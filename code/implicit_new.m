function h = implicit_new(h)
    function flag = l2()
        flag = sqrt(dr'*invDltSq*dr)>delta;
    end
    function flag = linf()
        flag = max(abs(sqrtm(invDltSq)*dr))>delta;
    end
    function flag = R1()
        flag = false;
        error('incomplete rule');
    end
    function flag = analR1()
        flag = false;
        error('incomplete rule');
    end
    function flag = Sense()
        flag = false;
        error('incomplete rule');
    end
    %
    [ m,n,p,U,D,Q0,q0,H] = getInitData( h);
    nsubiter = h.method.nsubiter;
    DltSq = h.method.DeltaSq;
    invDltSq = inv(DltSq);
    delta = h.rule.delta;
    %
    switch h.rule.name
        case 'l2'
            events = @l2;
        case 'linf'
            events = @linf;
        case 'R1'
            events = @R1;
        case 'analR1'
            events = @analR1;
        case 'Sense'
            events = @Sense;
        otherwise
            error('undefined rule');
    end
    tic;
    flag = true;
    lambda = 1e5;
    while flag
        xi = zeros(m,1);
        dxi = zeros(m,1);
        dxidr = zeros(m,n);
        Q = Q0;
        %
        nIter = 0;
        while flag
            for j = 1:nsubiter
                eHxj = exp(H*(xi+dxi));
                Q0_xj = diag(eHxj)*Q0;
                Q_xj = diag(eHxj)*Q;
                tld = diag(U'*D*eHxj);
                beta = (Q0_xj'*H*Q_xj+lambda*tld*DltSq*tld)\...
                    (p - Q0_xj'*(1-H*dxi));
                dr = lambda * DltSq * tld * beta;
%                 fprintf('%4i > %3i: norm = %6.10f\n',nIter,j,...
%                     norm( dxj - Q_xj*beta)/norm(dxj));
                dxi = Q_xj*beta;
                Q = Q0 + D*U*diag(dr);  
            end
            xi = xi + dxi;
            eHxi = exp(H*xi);
            Q_xi = diag(eHxi)*Q;
            M = diag(U'*D*eHxi);
            N = diag(eHxi)*D*U*diag(beta);
            A = [H - H*diag(dxi)*H,  -H*Q_xi, -H*N;...
                 -(H*Q_xi)'       , zeros(n),   -M;...
                 -(H*N)'          ,      -M', invDltSq/lambda];
            A = 0.5*(A+A');
            ddr = A\[H*dxidr; zeros(n); invDltSq/lambda];
            dxidr = ddr(1:m,:);
            drdr = ddr((end-n+1):end,:);
            %    
            flag = events();
            nIter = nIter+1;
            if toc>3
                lambda = lambda/10;
                warning('lambda')
                tic;
                break
            end
        end
    end
    h.method.niter = nIter;
    h.result.xi = xi;
    h.result.dxi = dxidr;
    h.result.r = h.method.r0 + dr;
    h.result.time = toc;
end
