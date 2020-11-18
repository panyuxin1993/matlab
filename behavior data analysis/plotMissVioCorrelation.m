%比较miss和violation session by session的correlation
clear;
%  dir='D:\xulab\behavior\CNO experiment\infusion\bilateral\';
dir='D:\xulab\behavior\';
% group={'saline', 'CNO'};
group={'pyx017'};
color={'k','r'}; 
A=figure;
for n=1:length(group)
    fPlotMissVioCorrelation(A,dir,group{n},color{n});
end

function [curve]=fPlotMissVioCorrelation(fig,dir,group,color)%one session one point
    [animal_name,dataChoiceResult]=fFindChoiceMat(strcat(dir,group));
    nsample=size(dataChoiceResult,1)-1;
    miss_vio=zeros(nsample,2);%第一列是miss,第二列是vio rate
    for i=1:nsample
        miss=(dataChoiceResult{i}(:,2)==3);
        violation=(dataChoiceResult{i}(:,2)==4);
        miss_vio(i,1)=nanmean(miss);
        miss_vio(i,2)=nanmean(violation);
        if miss_vio(i,2)==0 %vio=0 means allow violation
            miss_vio(i,1)=nan;
            miss_vio(i,2)=nan;
        end
    end
%     ffun=fittype('a+b.*x','independent','x');
%     f =fit(miss_vio(:,1),miss_vio(:,2),ffun,'Startpoint',[1,1]);
    x=min(miss_vio(:,1)):0.01:max(miss_vio(:,1));
    X=[ones(size(miss_vio,1),1) miss_vio(:,1)];
    [b,bint,r,rint,stats] = regress(miss_vio(:,2),X);
    y=b(1)+b(2)*x;
    r2=stats(1);
    p=stats(3);
    figure(fig);
    curve=scatter(miss_vio(:,1),miss_vio(:,2),color);
    hold on;
    plot(x,y,color,'LineWidth',1);
    hold on;
    xlabel('Miss');
    ylabel('Violation');
    set(gca,'FontName','Arial','FontSize',14);
    box off;
    text(0.3,0.3,strcat('r2=',num2str(r2),' p=',num2str(p)),'Color',color);
end