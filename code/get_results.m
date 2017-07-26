%% Load data
load ./data/data_update_may2015.mat
%% Different methods results
alpha = 0.12;
T = 30;
tau = 1e-4;
convpnt = 60;
UFR = 0.0365;
kernel = 'SW';
prof = data.profile;
opt1 = {'functional','original','nsubiter',5,...
    'UFR',UFR,'kernel',kernel,...
    'profile',prof};
opt2 = {'method','cauchy','rule','l2','delta',sqrt(0),...
    'functional','new','nsubiter',5,'norm','simple',...
    'UFR',UFR,'kernel',kernel,...
    'profile',prof};
opt3 = {'functional','new','nsubiter',5,...
    'UFR',UFR,'kernel',kernel,...
    'profile',prof};
cvp = {'tau',tau,'convpnt',convpnt};
results=cell(length(data.Date),3);
for i=1:length(data.Date)
    dt=data.Date{i};
    results{i,1} = SW(data,dt,opt1{:},cvp{:},'alpha',alpha);
    results{i,2} = SW(data,dt,opt2{:},cvp{:},'alpha',results{i,1}.method.alpha);
    results{i,3} = SW(data,dt,opt3{:},cvp{:},'alpha',results{i,2}.method.alpha);
    fprintf('Step %d of %d.\n',i,length(data.Date));
%     clf
%     subplot(1,2,1);
%     plotSW(results{i,3},'time',T,'color',[0 0 0]);
%     ylim('manual');
%     hold on;
%     plotSW(results{i,2},'time',T,'color',[1 0.4 0]);
%     legend('location','southeast')
%     title(results{i,3}.data.date)
%     plotrates(results{i,2});
%     pause(0.1)
end
%% Sensitivity results 
close all
x=[];
x(~logical(mod(1:1:sum(data.liquid_mask),2)))=-0.0001;
x(logical(mod(1:1:sum(data.liquid_mask),2)))=0.0001;
dr=zeros(size(data.tenor));
dr(data.liquid_mask)=x;
dr1=dr';
dr2=-dr';
dr3=0.0001*double(dr~=0)';
dr4=-0.0001*double(dr~=0)';
coef = 1e-2;
dr1 = coef*dr1; dr2 = coef*dr2;
dr3 = coef*dr3; dr4 = coef*dr4;
data1=data;
data2=data;
data3=data;
data4=data;

sens=cell(length(data.Date),12);
for i=1:length(data.Date)
    dt=data.Date{i};
    %
    data1.PX_LAST(i,:)=data1.PX_LAST(i,:)+100*dr1;
    data1.PX_ASK(i,:)=data1.PX_ASK(i,:)+100*dr1;
    data1.PX_BID(i,:)=data1.PX_BID(i,:)+100*dr1;
    data1.PX_MID(i,:)=data1.PX_MID(i,:)+100*dr1;

    data2.PX_LAST(i,:)=data2.PX_LAST(i,:)+100*dr2;
    data2.PX_ASK(i,:)=data2.PX_ASK(i,:)+100*dr2;
    data2.PX_BID(i,:)=data2.PX_BID(i,:)+100*dr2;
    data2.PX_MID(i,:)=data2.PX_MID(i,:)+100*dr2;

    data3.PX_LAST(i,:)=data3.PX_LAST(i,:)+100*dr3;
    data3.PX_ASK(i,:)=data3.PX_ASK(i,:)+100*dr3;
    data3.PX_BID(i,:)=data3.PX_BID(i,:)+100*dr3;
    data3.PX_MID(i,:)=data3.PX_MID(i,:)+100*dr3;

    data4.PX_LAST(i,:)=data4.PX_LAST(i,:)+100*dr4;
    data4.PX_ASK(i,:)=data4.PX_ASK(i,:)+100*dr4;
    data4.PX_BID(i,:)=data4.PX_BID(i,:)+100*dr4;
    data4.PX_MID(i,:)=data4.PX_MID(i,:)+100*dr4;
    %
    sens{i,1} = SW(data1,dt,opt1{:},cvp{:},'alpha',results{i,1}.method.alpha);
    sens{i,2} = SW(data2,dt,opt1{:},cvp{:},'alpha',results{i,1}.method.alpha);
    sens{i,3} = SW(data3,dt,opt1{:},cvp{:},'alpha',results{i,1}.method.alpha);
    sens{i,4} = SW(data4,dt,opt1{:},cvp{:},'alpha',results{i,1}.method.alpha);

    sens{i,5} = SW(data1,dt,opt2{:},cvp{:},'alpha',results{i,2}.method.alpha);
    sens{i,6} = SW(data2,dt,opt2{:},cvp{:},'alpha',results{i,2}.method.alpha);
    sens{i,7} = SW(data3,dt,opt2{:},cvp{:},'alpha',results{i,2}.method.alpha);
    sens{i,8} = SW(data4,dt,opt2{:},cvp{:},'alpha',results{i,2}.method.alpha);

    sens{i,9} = SW(data1,dt,opt3{:},cvp{:},'alpha',results{i,3}.method.alpha);
    sens{i,10} = SW(data2,dt,opt3{:},cvp{:},'alpha',results{i,3}.method.alpha);
    sens{i,11} = SW(data3,dt,opt3{:},cvp{:},'alpha',results{i,3}.method.alpha);
    sens{i,12} = SW(data4,dt,opt3{:},cvp{:},'alpha',results{i,3}.method.alpha);
    % 
    fprintf('Step %d of %d.\n',i,length(data.Date));
end
