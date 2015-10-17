% script for producing a result array of structures
daygrid = 32:2450;
N = length(daygrid);
clear arrh0 arrh1;
arrh0(1,N) = getstruct();
arrh1(1,N) = getstruct();
fprintf('%6.2f%%\n',0 )
for i = 1:N
    fprintf('\b\b\b\b\b\b\b\b%6.2f%%\n',i/N*100)
    day = daygrid(i);
    arrh0(i) = SW(data,day,1:30); % 1:20 for EIOPA strategy 
    arrh1(i) = SW(data,day,1:30,'method','implicit','rule','Sense',...
        'delta',0.1*arrh0(i).result.sense); % also one can use the sense value 0.003.
end
%% compare annuity for original and new 
ann = zeros(2,N);
for i=1:N
    ann(1,i) = arrh0(i).result.annuity;
    ann(2,i) = arrh1(i).result.annuity;
end
%plot(daygrid,ann)
grid on
plot(daygrid(1:end-1),diff(ann(1,:),1),'.',daygrid(1:end-1) , diff(ann(2,:)),'+r')