function sense = sensefnc(h)
    switch h.method.sensefncname
        case 'S1'
            f = @(grad,DltSq)sum(1e-4*abs(grad));
        case 'S2'
            f = @(grad,DltSq)sum(abs(sqrtm(DltSq)*grad));
        case 'S3'
            f = @(grad,DltSq)sqrt(grad'*DltSq*grad); 
        otherwise
            error('not correct name')
    end
    switch h.method.denormsense
        case 'y'
            f = @(grad,DltSq,ann)f(grad,DltSq)/ann;
        case 'n'
            f = @(grad,DltSq,ann)f(grad,DltSq);
        otherwise
            error(['undefined value "',h.method.denormsense,'"']);
    end
%     sense = f(h.result.grad, h.method.DeltaSq, h.result.annuity);
    sense = f(h.result.grad, (4e-8)*eye(size(h.method.DeltaSq)), h.result.annuity);
end

