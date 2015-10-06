function fnc = sensefnc(h)
    switch h.method.sensefncname
        case 'S1'
            fnc = @(grad,DltSq)sum(1e-4*abs(grad));
        case 'S2'
            fnc = @(grad,DltSq)sum(abs(sqrtm(DltSq)*grad));
        case 'S3'
            fnc = @(grad,DltSq)sqrt(grad'*DltSq*grad); 
        otherwise
            error('not correct name')
    end
    switch h.method.denormsense
        case 'y'
            fnc = @(grad,DltSq,ann)fnc(grad,DltSq)/ann;
        case 'n'
            fnc = @(grad,DltSq,ann)fnc(grad,DltSq);
        otherwise
            error(['undefined value "',h.method.denormsense,'"']);
    end
end

