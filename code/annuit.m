daygrid = 1600:1650;
N = length(daygrid);
arrh0(1,N) = h;
arrh1(1,N) = h;
for i = 1:length(daygrid)
    disp(i)
    day = daygrid(i);
    go
    %ann(:,i) = [h0.result.annuity;h1.result.annuity];
    arrh0(i) = h0;
    arrh1(i) = h1;
end
%%
ann = zeros(2,N);
for i = 1:N
    ann(:,i) = [arrh0(i).result.annuity;arrh1(i).result.annuity];
end
plot(daygrid,ann);
%%
hold on;
for i = 1:N
    plotSW(arrh1(i),[i/N 1-i/N 0]);
end
