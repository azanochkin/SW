daygrid = 1:2450;
N = length(daygrid);
arrh0(1,N) = getstruct();
arrh1(1,N) = getstruct();
for i = 1:N
    if (i/10) == round(i/10)
        disp(i);
    end
    day = daygrid(i);
    arrh0(i) = SW(data,day,1:30);
    arrh1(i) = SW(data,day,1:30,'method','implicit','rule','Sense',...
        'delta',0.003);%0.1*h0.result.sense);
end
%%
ann = zeros(2,N);
for i=1:N
    ann(1,i) = arrh0(i).result.annuity;
    ann(2,i) = arrh1(i).result.annuity;
end
plot(daygrid,ann)
grid on
%plot(daygrid(1:end-1),ann(:,2:end) - ann(:,1:end-1))