%%
[sen,arrh,res] = fixSense(data,1600:1650,1:30,1);
%%
dayGrid = 2032:2450;
N = numel(dayGrid);
res = zeros(1,N);
val = zeros(1,N);

senGrid = 0.003;%linspace(0.0001,0.005,10);
meanres = zeros(size(senGrid));
clear arrh0 arrh1;
arrh0(N) = getstruct();
arrh1(N) = getstruct();
for i = 1:numel(senGrid)
    sen = senGrid(i);
    for j = 1:N
        disp(j)
        day = dayGrid(j);
        arrh0(j) = SW(data,day,1:30);
        %h1 = SW(data,day,1:30,'method','implicit','rule','R1','delta',sen);
        h1 = SW(data,day,1:30,'method','implicit','rule','Sense','delta',sen,'denormsense','y');
        %h1 = SW(data,day,1:30,'method','implicit','rule','Sense','delta',0.1*h0.result.sense);
        arrh1(j) = h1;
        dr = h1.result.r - h1.method.r0;
        %res(j) = (sqrt(dr'*(h1.method.DeltaSq\dr)))/sqrt(30);
        %res(j) = mean(abs(sqrtm(h1.method.DeltaSq)\dr));
        res(j) = max(abs(sqrtm(h1.method.DeltaSq)\dr));
        val(j) = max(abs(sqrt(diag(h1.method.DeltaSq))));
        
    end
    meanres(i) = mean(res);
end

%%
wndSize = 10;
plot(medianfilter(wndSize,res),'displayname','normalized residual')
hold on
plot(medianfilter(wndSize,clres),'g','displayname','non-normalized residual')
%plot(medianfilter(wndSize,(ann - sum((1.042).^(-(1:100))))/4),'m')
%plot(medianfilter(wndSize,(ann)),'m')
plot(medianfilter(wndSize,2e3*val),'r','displayname','volatility')
%plot(medianfilter(wndSize,1.2e2*fvc),'k')
plot(medianfilter(wndSize,1e1*smt),'c','displayname','smootness')
legend show
%%
N = length(arrh1);
ann = zeros(1,N);
smt = zeros(1,N);
clres = zeros(1,N);
for j = 1:N
    dr = arrh1(j).result.r - arrh1(j).method.r0;
    %
    fnc = @(x)30*mean(abs(x)); 
    vec = sqrtm(arrh1(j).method.DeltaSq)\dr;
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
    vec = sqrt(diag(arrh1(j).method.DeltaSq));
    %val(j) = max(abs(vec));
    %val(j) = mean(abs(vec));
    %val(j) = norm(vec);
    val(j) = fnc(vec);
    %
    ann(j) = arrh1(j).result.annuity;
    xi = arrh1(j).result.xi;
    %smt(j) = xi'*Wilson_Heart(0.1,1:30,1:30)*xi;
    smt(j) = xi'*Wilson_Smooth(0.1,1:30,1:30)*xi;
end
[~,~,fvc] = getrates(arrh1,[30,30,30]);
fvc = fvc(1,:);