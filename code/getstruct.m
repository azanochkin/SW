function h = getstruct()
    data = struct(...
        'day_number',[],...
        'date','',...
        'r_bid',[],...
        'r_ask',[],...
        'r_mid',[],...
        'r_last',[],...
        'tenor',[],...
        'mask',[],...
        'u',[],...
        'u_ann',(1:100)',...
        'profile',ones(100,1)...
    );
    method = struct(...
        'functional','original',...
        'name','original',...
        'DeltaSq',[],...
        'normname','simple',...
        'r0_type','last',...
        'kernel','SW',...
        'tau',1e-4,...
        'tautol',1e-5,...
        'convpnt',60,...
        'fixalpha',true,...
        'alpha',0.1,...
        'ufr',log(1.0405),...
        'ufrvar',0,...
        'r0',[],...
        'sensefncname','S3',...
        'denormsense','y',...
        'nsubiter',5,...
        'niter',10,...
        'fndderiv',true...
    );
    rule = struct(...
        'lambda',[],...
        'name','l2',...
        'delta',1 ...
    );
    result = struct(...
        'r',[],...
        'xi',[],...
        'dxidr',[],...
        'dxidp',[],...
        'grad',[],...
        'sense',[],...
        'time',[],...
        'annuity',[] ...
    );
    h = struct(...
        'data',data,...
        'method',method,...
        'rule',rule,...
        'result',result...
    );
end

