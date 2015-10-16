% get statistic for annuity,smoothness, residual, volatility. Require arrh1
% be a result structure(impl_*). 
N = length(arrh1);
ann = zeros(1,N);
smt = zeros(1,N);
clres = zeros(1,N);
res = zeros(1,N);
vol = zeros(1,N);
for j = 1:N
    %snsbeg(j) = arrh0(j).result.sense;
    dr = arrh1(j).result.r - arrh1(j).method.r0; %residual
    %
    fnc = @(x)max(abs(x));%30*mean(abs(x)); % normalization function
    vec = sqrtm(arrh1(j).method.DeltaSq)\dr; % normalized residual
    %res(j) = max(abs(vec));
    %res(j) = mean(abs(vec));
    %res(j) = norm(vec);
    res(j) = fnc(vec);
    %
    %vec = 1e-4*(arrh1(j).method.DeltaSq)\dr;
    vec = 1e4*dr;
    %clres(j) = max(abs(vec));
    %clres(j) = mean(abs(vec));
    %clres(j) = norm(vec);
    clres(j) = fnc(vec);
    %
    vec = sqrt(diag(arrh1(j).method.DeltaSq)); %volatility vector
    %val(j) = max(abs(vec));
    %val(j) = mean(abs(vec));
    %val(j) = norm(vec);
    vol(j) = 1e4*fnc(vec);
    %
    ann(j) = arrh1(j).result.annuity;
    xi = arrh1(j).result.xi;
    %smt(j) = xi'*Wilson_Heart(0.1,1:30,1:30)*xi;
    smt(j) = xi'*Wilson_Smooth(0.1,1:30,1:30)*xi;
end
[~,~,fvc] = getrates(arrh1,[30,30,30]);
fvc = fvc(1,:);
%% plot statistics 
wndSize = 10; % filter window
plot(medianfilter(wndSize,res),'displayname','normalized residual')
hold on
plot(medianfilter(wndSize,clres),'g','displayname','non-normalized residual')
%plot(medianfilter(wndSize,(ann - sum((1.042).^(-(1:100))))/4),'m')
plot(medianfilter(wndSize,(ann)/10),'m','displayname','annuity')
plot(medianfilter(wndSize,vol),'r','displayname','volatility')
%plot(medianfilter(wndSize,1.2e2*fvc),'k')
%plot(medianfilter(wndSize,1e0*smt),'c','displayname','smootness')
%plot(medianfilter(wndSize,3e2*snsbeg),'c','displayname','original sense')
legend show
