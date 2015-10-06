med = 1:30;
for ind=1:30
    spread = data.PX_ASK(2338:end,ind) - data.PX_BID(2338:end,ind);
    med(ind) = median(spread(spread>0));
end
%%
diap = 2430:2450;
dr = data.PX_MID(diap,:) - data.PX_MID(diap-1,:);
c = zeros(1,30);
for i=1:30
    co = cov(dr(1:end-1,i),dr(2:end,i));
    c(i) = co(2,1);
end
2*sqrt(max(0,-c));