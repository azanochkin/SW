daygrid = 1:20;
N = length(daygrid);
arrh0(1,N) = getstruct();
arrh1(1,N) = getstruct();
u = 1:30;
for i = 1:length(daygrid)
    if (i/10) == round(i/10)
        disp(i);
    end
    day = daygrid(i);
    arrh0(i) = SW(data,day,u);
    arrh1(i) = SW(data,day,u,'method','implicit','rule','Sense',...
        'delta',0.003);%0.1*h0.result.sense);
end