function fnc = sensefnc( name )
    switch name 
        case 'S1'
            fnc = @(grad,DltSq,ann)sum(1e-4*abs(grad));
        case 'S2'
            fnc = @(grad,DltSq,ann)sum(abs(sqrtm(DltSq)*grad));
        case 'S3'
            fnc = @(grad,DltSq,ann)sqrt(grad'*DltSq*grad); 
        case 'S4'
            fnc = @(grad,DltSq,ann)sqrt(grad'*DltSq*grad)/ann;
        otherwise
            error('not correct name')
    end
end

