function h = getstruct()
    data = struct('day_number',[],'date','','r_bid',[],'r_ask',[],'r_mid',[],...
        'u',[],'u_ann',(1:100)','profile',ones(100,1));
    method = struct('name','original','DeltaSq',[],'normname','spread','alpha',0.1,'ufr',log(1.042),'r0',...
        [],'sensefncname','S3','denormsense','y','nsubiter',2,'niter',10);
    rule = struct('name','l2','delta',1);
    result = struct('r',[],'xi',[],'dxi',[],'grad',[],'sense',[],'time',[],...
        'annuity',[]);
    h = struct('data',data,'method',method,'rule',rule,'result',result);
end

