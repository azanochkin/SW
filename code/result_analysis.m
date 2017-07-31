%% Value dynamics
ann_val=nan(size(results));
dmr = nan(size(results));
dr = nan(length(results),length(results{1,2}.result.r));
for i=1:size(results,1)
    dr(i,results{i,2}.data.tenor) = results{i,2}.result.r - results{i,2}.data.r_last;
    for j=1:size(results,2)
        ann_val(i,j)=results{i,j}.result.annuity;
        dmr(i,j) = max(abs(results{i,j}.result.r - results{i,j}.data.r_last));
    end
end
colors=[0 0 0;...
        0 0 1;...
        1 .5 0;...
        0 .5 0;...
        0 1 1];
lgd={'SW\_original','SW+Tikhonov\_linf','SW\_modified'};

figure

subplot(1,2,1)
for i=1:size(ann_val,2)
plot(datenum(data.Date,'dd.mm.yyyy'),ann_val(:,i),'color',colors(i,:))
hold on
end
legend(lgd);
xlim([min(datenum(data.Date,'dd.mm.yyyy')) max(datenum(data.Date,'dd.mm.yyyy'))])
datetick('keeplimits')
legend('location','southeast')

subplot(1,2,2)
diff_val=100*log(ann_val(2:end,:)./ann_val(1:end-1,:));

for i=1:size(ann_val,2)
    plot(diff_val(:,i),'color',colors(i,:),'DisplayName',[lgd{i},': std=',num2str(std(diff_val(:,i)),'%3.3f'),'%'])
    hold on
end
legend('location','southeast')
xlim([1 length(diff_val(:,1))])
ylabel('%')
%% Sensitivity 
mat=nan(size(sens)); 
mat2=nan(size(sens)); 
% mat3=nan(size(sens)); 
dmr = [dr1;dr2;dr3;dr4];
for day=1:558
for i=1:12
    k = fix((i-1)/4)+1;
    j = mod(i-1,4)+1;
    mat(day,i)=100*(sens{day,i}.result.annuity/results{day,k}.result.annuity-1);
    ism = ismember(data.tenor,results{day,k}.data.tenor);
    mat2(day,i) = 100*(dmr(j,ism)*results{day,k}.result.grad)/results{day,k}.result.annuity;
%     mat3(day,i) = 100*(sens{day,i}.method.alpha / results{day,k}.method.alpha - 1);
end
end
%% Impact of data changes
dmr=[dr1;dr2;dr3;dr4];
for i=1:4
   subplot(2,4,i)
   plot(mat(:,[i,i+4,i+8]))
   legend({'Orig\_Original','New\_Tikhonov','New\_Orig'})
   ylabel('%')  
   yl = ylim();
   
   subplot(2,4,i+4)
   plot(mat2(:,[i,i+4,i+8]))
   legend({'Orig\_Original','New\_Tikhonov','New\_Orig'})
   ylabel('%')  
   ylim(yl);
end