function [ A ] = fPlotCorrectRate( animal_name,dataChoice,sessionNum )
%FPLOTCORRECTRATE Summary of this function goes here
%   Detailed explanation goes here
%根据时间（月份日期）画出各种正确率和bias
performance = fCorrectRate( dataChoice );
% answerPeriod=dataChoice(:,1001,2);
if sessionNum<size(performance,1)
    performanceL_R=performance(size(performance,1)-sessionNum+1:end,:);
    date=performance(size(performance,1)-sessionNum+1:end,11);
    clickRate_L=dataChoice(size(performance,1)-sessionNum+1:end,1002,1);%列向量
    performanceCircle=zeros(sessionNum,size(performance,2));
    performanceDot=zeros(sessionNum,size(performance,2));
    temp=zeros(size(performance,2));
else
    date=performance(:,11);
    clickRate_L=dataChoice(:,1002,1);%列向量
    performanceL_R=performance;
    performanceCircle=zeros(size(performance));
    performanceDot=zeros(size(performance));
    temp=zeros(size(performance,2));
end
for i=1:size(performanceL_R,1)
    if clickRate_L(i)==125
        temp=performance(i,:);
        performanceL_R(i,5:8)=temp(9:12);
        performanceL_R(i,9:12)=temp(5:8);
        performanceDot(i,:)=performanceL_R(i,:);
        performanceCircle(i,:)=nan;
    else
        performanceCircle(i,:)=performanceL_R(i,:);
        performanceDot(i,:)=nan;
    end
end
date_d=mod(date,100);
date_m=(date-date_d)/100;
x=strcat(num2str(date_m),'-',num2str(date_d));
biasNoMiss=performance(:,13);

%画两张图，一张正确率，一张miss_rate
A=figure(1);
hold on;
%suptitle(animal_name);
%ha = tight_subplot(1,2,[1 1],[1 1],[1 1]);
set(gcf, 'position', [0 0 1200 400]);%控制fig尺寸
%set(gca,'position',[0.1,0.1,0.9,0.9]);
%颜色统一：左边绿色，右边红色，总蓝色
%规则标注：左边20-右边125（空心点），左边125-右边20（实心点）
subplot(1,3,1);
%axes(ha(1)); 
plot(performanceL_R(:,2),'k','LineWidth',2);	%总正确率，不含miss
%set(gca,'XTickLabel',x);     %不能用char型的值作为横坐标，而应该先画出来，把横坐标改成字符型
hold on;
plot(performanceL_R(:,6),'b','LineWidth',2);    %左边正确率，不含miss
hold on;
plot(performanceL_R(:,10),'r','LineWidth',2);    %右边正确率，不含miss
hold on;
scatter(1:size(performanceDot,1),performanceDot(:,2),50,'k','filled');%实心点
hold on;
scatter(1:size(performanceDot,1),performanceDot(:,6),50,'b','filled');
hold on;
scatter(1:size(performanceDot,1),performanceDot(:,10),50,'r','filled');
hold on;
scatter(1:size(performanceCircle,1),performanceCircle(:,2),50,'k');%空心点
hold on;
scatter(1:size(performanceCircle,1),performanceCircle(:,6),50,'b');%空心点
hold on;
scatter(1:size(performanceCircle,1),performanceCircle(:,10),50,'r');%空心点
hold on;
%legend('% correct','% correct left','% correct right','Location','EastOutside' );
%plot(performanceL_R(:,5)./performanceL_R(:,9),'k','LineWidth',2);%左边/右边trial数量
%plot(biasNoMiss,'g','LineWidth',2);           %不含miss的bias
%legend('correct rate with miss','bias with miss','correct rate without miss','bias without miss','Location','BestOutside');
xlabel('Training day','FontName','Arial','FontSize',14);
ylabel('% Correct','FontName','Arial','FontSize',14);
title(['Performance of mouse ', animal_name],'FontName','Arial','FontWeight','Bold','FontSize',16)
set(gca,'FontName','Arial','FontSize',14)%设置坐标轴刻度字体名称，大小
plot([1,size(performanceL_R,1)],[0.8,0.8],'--k','LineWidth',1);
plot([1,size(performanceL_R,1)],[0.5,0.5],'--k','LineWidth',1);
plot([1,size(performanceL_R,1)],[1,1],'-k','LineWidth',1);
% plot([0,size(performance,1)],[0,0],'-k','LineWidth',1);
% plot([answerPeriod(1),answerPeriod(1)],[0,1],'--k','LineWidth',1);解决训练时的bug，标注出让answerPeriod变成300的情况
% plot([answerPeriod(end),answerPeriod(end)],[0,1],'--k','LineWidth',1);

subplot(1,3,2);
%axes(ha(2)); 
plot(performanceL_R(:,3),'k','LineWidth',2);	%总miss_rate
hold on;
plot(performanceL_R(:,7),'b','LineWidth',2);	%左miss_rate
hold on;
plot(performanceL_R(:,11),'r','LineWidth',2);	%右miss_rate
hold on;
legend('Total','Left','Right');
xlabel('Training day','FontName','Arial','FontSize',14);
ylabel('% Miss','FontName','Arial','FontSize',14);
title(['Miss rate of mouse ' ,animal_name],'FontName','Arial','FontWeight','Bold','FontSize',16)
set(gca,'FontName','Arial','FontSize',14)%设置坐标轴刻度字体名称，大小

subplot(1,3,3);
%axes(ha(2)); 
plot(performanceL_R(:,4),'k','LineWidth',2);	%总violation_rate
hold on;
plot(performanceL_R(:,8),'b','LineWidth',2);	%左violation_rate
hold on;
plot(performanceL_R(:,12),'r','LineWidth',2);	%右violation_rate
hold on;
%legend('Total','Left','Right');
xlabel('Training day','FontName','Arial','FontSize',14);
ylabel('% Violation','FontName','Arial','FontSize',14);
title(['Violation rate of mouse ' ,animal_name],'FontName','Arial','FontWeight','Bold','FontSize',16)
set(gca,'FontName','Arial','FontSize',14)%设置坐标轴刻度字体名称，大小

end

function [ performance ] = fCorrectRate( dataChoice )
%FCORRECTRATE Summary of this function goes here
%   Detailed explanation goes here
% 输入选择矩阵数据，计算某日的正确率（包括两侧），相对偏好bias.
n_day=size(dataChoice,1);     %获得训练天数
%每一列分别表示总（7000Hz侧/28000Hz侧）trial数量，正确率（不包含miss），miss_rate，violation
%rate(不包括miss）
%以及bias,日期
performance=zeros(n_day,14);    
for i=1:n_day
    n_trial=find(dataChoice(i,:,1)==0); 
    if isempty(n_trial)%做完1000trial时n_tial=null;
        performance(i,1)=1000;
    else
        performance(i,1)=n_trial(1)-1;  %总trial数
    end
%     m_7k=double(dataChoice(i,:,1)==7000);  %同时要把逻辑矩阵转化为数值矩阵
%     m_28k=double(dataChoice(i,:,1)==28000);
    m_20=double((dataChoice(i,:,1)<50).*(dataChoice(i,:,1)>0));%行向量用50做标准的话，就能使得probe也进行计算。注意防止出现0为没做的trial
    m_125=double(dataChoice(i,:,1)>50);
    m_correct=double(dataChoice(i,:,2)==1);
    m_error=double(dataChoice(i,:,2)==2);
    m_miss=double(dataChoice(i,:,2)==3);
    m_violation=double(dataChoice(i,:,2)==4);
    m_nomiss=m_correct+m_error+m_violation;
    m_total=m_nomiss+m_miss;
    performance(i,2)=sum(m_correct)/sum(m_correct+m_error);   %正确率
    performance(i,3)=sum(m_miss)/sum(m_total);     %miss_rate
    performance(i,4)=sum(m_violation)/sum(m_nomiss);   %不包括miss的violation
    performance(i,5)=sum(m_20);     %总20数    
    performance(i,6)=sum(m_correct*m_20')/((m_correct+m_error)*m_20');      %20正确率, .'才是转置
    performance(i,7)=sum(m_miss*m_20')/performance(i,5);        %20miss_rate
    performance(i,8)=sum(m_violation*m_20')/sum(m_nomiss*m_20');   %不包括miss的violation
    performance(i,9)=sum(m_125);    %总125数
    performance(i,10)=sum(m_correct*m_125')/((m_correct+m_error)*m_125');     %125正确率
    performance(i,11)=sum(m_miss*m_125')/performance(i,9);       %125miss_rate
    performance(i,12)=sum(m_violation*m_125')/sum(m_nomiss*m_125');   %不包括miss的violation
%     performance(i,4)=sum(m_7k);     %总7k数
%     performance(i,7)=sum(m_28k);    %总28k数
%     performance(i,5)=sum(m_correct*m_7k')/performance(i,4);      %7k正确率, .'才是转置
%     performance(i,6)=sum(m_error*m_7k')/performance(i,4);        %7k错误率
%     performance(i,8)=sum(m_correct*m_28k')/performance(i,7);     %28k正确率
%     performance(i,9)=sum(m_error*m_28k')/performance(i,7);       %28k错误率
    %归一化的bias 指数
    %performance(i,10)=(performance(i,5)-performance(i,8))/(performance(i,5)+performance(i,8));  %bias;
    %简单相减得到的bias
    performance(i,13)=(performance(i,5)-performance(i,8));
    performance(i,14)=dataChoice(i,1001,1);             %日期
end

end

