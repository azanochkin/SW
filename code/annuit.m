% script for producing a result array of structures
daygrid = 2432:2450;
N = length(daygrid);
clear arrh0 arrh1;
arrh0(1,N) = getstruct();
arrh1(1,N) = getstruct();
fprintf('%6.2f%%\n',0 )
for i = 1:N
    fprintf('\b\b\b\b\b\b\b\b%6.2f%%\n',i/N*100)
    day = daygrid(i);
    arrh0(i) = SW(data,day,1:20); % 1:20 for EIOPA strategy 
    arrh1(i) = SW(data,day,1:20,'method','implicit','rule','l2','delta',sqrt(20),'norm','volatility');
%     arrh1(i) = SW(data,day,1:20,'method','implicit','rule','Sense',...
%          'norm','volatility','delta',0.1);%0.04*arrh0(i).result.sense); % also one can use the sense value 0.003.
end
%% compare annuity for original and new 
ann = zeros(2,N);
sen = zeros(2,N);
for i=1:N
    ann(1,i) = arrh0(i).result.annuity;
    ann(2,i) = arrh1(i).result.annuity;
    sen(1,i) = arrh0(i).result.sense;
    sen(2,i) = arrh1(i).result.sense;
    if sen(1,i)<0
        sen(:,i) = sen(:,i-1);
        ann(:,i) = ann(:,i-1);
    end
    if sen(2,i)>0.007
        sen(:,i) = sen(:,i-1);
        ann(:,i) = ann(:,i-1);
    end
end
%plot(daygrid,ann(1,:),'r',daygrid,ann(2,:),'b')
grd = 1:N;
dt = datenum(data.Date(daygrid(grd)),'dd.mm.yyyy');
subplot(2,1,1)
plot(dt',ann(1,grd),'r',dt',ann(2,grd),'b')
xlabel 'Date'
ylabel 'Annuity, m euro'
title 'Annuity'
legend( 'original', 'regularized')
legend('location','best')
datetick
subplot(2,1,2)
sen = 1e2*sen;
plot(dt',sen(1,grd),'r',dt',sen(2,grd),'b')
xlabel 'Date'
ylabel 'Sensitivity, %'
title 'Sensitivity'
legend( 'original', 'regularized')
datetick
grid on
%plot(daygrid(1:end-1),diff(ann(1,:),1),'.',daygrid(1:end-1) , diff(ann(2,:)),'+r')
%%
std(diff(ann(1,1:end)))
std(diff(ann(2,1:end)))
% period = 1:1400;
% std(diff(ann(1,period)))
% std(diff(ann(2,period)))