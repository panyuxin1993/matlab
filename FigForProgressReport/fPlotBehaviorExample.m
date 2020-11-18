function [ fig ] = fPlotBehaviorExample( fig, position )
%FPLOTBEHAVIOREXAMPLE plot an example session
%   include cor rate, vio rate, RT, curve, licking raster, etc.
ndelaygroup=4;
path='D:\FP\pyx172_20190508';
filename='D:\FP\pyx172_20190508\PYX172_20190508_FP.mat'; 
animal='pyx172';
load(filename);    
dataChoice=zeros(length(SessionResults),4);%第1-3列分别表示当前刺激，Action_choice，Delay_duration，reaction time  
n=1; %n个trial 
clickright=SessionSettings{1}.clickRate_;  
for i=1:length(SessionResults)
    if  isfield(SessionResults{i},'Stim_clickRate') && ~isempty(SessionResults{i}.Stim_clickRate)
        if SessionResults{i}.Stim_clickRate<126%存在125170500这种类型，和后面volume合并的
            dataChoice(i,1)= SessionResults{i}.Stim_clickRate;
        else 
            dataChoice(i,1)= floor(SessionResults{i}.Stim_clickRate/10^6);%存在125170这种类型
        end
    elseif isfield(SessionResults{i},'Stim_clickRate') && isempty(SessionResults{i}.Stim_clickRate)
        continue;
    elseif ~isfield(SessionResults{i}  ,'Stim_clickRate') %'clicks20'这种类型
        dataChoice(i,1)=  str2double(SessionResults{i}.Stim_Type(7:end));
    end
    dataChoice(i,2)=SessionResults{i}.Action_choice;%0左，1右，2miss，3violation
    dataChoice(i,3)=SessionResults{i}.Delay_duration;  
    dataChoice(i,4)=SessionResults{i}.Time_answer-SessionResults{i}.Time_delayOffset;
end
if dataChoice(1,1)>125 %有些文件记录出Bug,把之后的volume的数据直接连到click rate后面，导致该数字多了三位。
    dataChoice(:,1)=floor(dataChoice(:,1)/1000);  
end
maxdelay=max(dataChoice(:,3));
mindelay=min(dataChoice(:,3));
part=(maxdelay-mindelay)/ndelaygroup;
delay=dataChoice(:,3)<=(mindelay+part);
for i=1:ndelaygroup
    delaycol=(dataChoice(:,3)>(mindelay+i*part)).*(dataChoice(:,3)<=(mindelay+i*part+part));
    delay=[delay, delaycol];%汇总成二维数组，每一列代表不同长度的delay
end
boundary=50;
%boundary=200000;
if clickright==20
    left=(dataChoice(:,1)>boundary);
    right=(dataChoice(:,1)<boundary);
else
    left=(dataChoice(:,1)<boundary);
    right=(dataChoice(:,1)>boundary);
end
all=left+right;
side=[all,left,right];%汇总成二维数组，每一列代表不同边
correct=left.*(dataChoice(:,2)==0)+right.*(dataChoice(:,2)==1);
error=right.*(dataChoice(:,2)==0)+left.*(dataChoice(:,2)==1);
miss=(dataChoice(:,2)==2);
violation=(dataChoice(:,2)==3);
do=correct+error;
nomiss=do+violation;
indRT=correct.*(dataChoice(:,4)>0);%只看正确的trial

performDiffDelay=zeros(ndelaygroup,4,3);%第一维控制delay长度，第二维控制正确率/miss/violation/RT，第三维控制总，左，右
for i=1:ndelaygroup
    for j=1:3
        performDiffDelay(i,1,j)=100*sum(double(correct.*side(:,j).*delay(:,i)))/sum(double(do.*side(:,j).*delay(:,i)));
        performDiffDelay(i,2,j)=100*sum(double(miss.*side(:,j).*delay(:,i)))/sum(double(side(:,j).*delay(:,i)));
        performDiffDelay(i,3,j)=100*sum(double(violation.*side(:,j).*delay(:,i)))/sum(nomiss.*double(side(:,j).*delay(:,i)));
        temp=double(dataChoice(:,4).*indRT.*side(:,j).*delay(:,i));
        temp2=temp;
        temp2(temp==0)=nan;
        performDiffDelay(i,4,j)=nanmedian(temp2);
    end
end
 
figure(fig);
hold on;
%suptitle(animal_name);
%ha = tight_subplot(1,2,[1 1],[1 1],[1 1]);
set(fig, 'PaperPosition',position);%控制fig尺寸
%颜色统一：左边绿色，右边红色，总蓝色
%规则标注：左边20-右边125（空心点），左边125-右边20（实心点）
% suptitle(strcat(animal,'-',date));
y_label={'% Correct';'% Miss';'%Violation';'reaction time(ms)'};%{}是cell类型，不会出现矩阵长度不齐的问题。
% x_axis={strcat(num2str(mindelay),' - ',num2str(mindelay+part));...
%     strcat(num2str(mindelay+part),' - ',num2str(mindelay+2*part));...
%     strcat(num2str(mindelay+2*part),' - ',num2str(maxdelay));...
%     };
x_plot=mindelay+part/2:part:maxdelay-part/2;
subplot_ind=[2,0,3,5];
for i=[1,3,4]
    subplot(2,3,subplot_ind(i));
    H1=plot(x_plot,performDiffDelay(:,i,1),'k','LineWidth',2);	%总正确率，不含miss
    %set(gca,'XTickLabel',x);     %不能用char型的值作为横坐标，而应该先画出来，把横坐标改成字符型
    hold on;
%     H2=plot(x_plot,performDiffDelay(:,i,2),'b','LineWidth',2);    %左边正确率，不含miss
%     hold on;
%     H3=plot(x_plot,performDiffDelay(:,i,3),'r','LineWidth',2);    %右边正确率，不含miss
%     hold on;
%     H=[H1,H2,H3];
    if SessionSettings{1}.clickRate_==20 %画实心点
        scatter(x_plot,performDiffDelay(:,i,1),50,'k','filled');%实心点
        hold on;
%         scatter(x_plot,performDiffDelay(:,i,2),50,'b','filled');
%         hold on;
%         scatter(x_plot,performDiffDelay(:,i,3),50,'r','filled');
%         hold on;
    else
        scatter(x_plot,performDiffDelay(:,i,1),50,'k');%空心点
        hold on;
%         scatter(x_plot,performDiffDelay(:,i,2),50,'b'); 
%         hold on;
%         scatter(x_plot,performDiffDelay(:,i,3),50,'r');
%         hold on;
    end
    if i<4
        plot([mindelay,maxdelay],[50,50],'--k','LineWidth',1);
        hold on;
        set(gca, 'YLim',[0 100]);
        ylimmax=100;
    else
        ylimrange=get(gca,'Ylim');
        ylimmax=ylimrange(end);
    end
    set(gca,'ylim',[0,ylimmax]);
    xlabel('Delay(ms)');
    ylabel(y_label(i));
    box off; %取消图标外边框
    set(gca, 'XLim',[mindelay maxdelay]);  
    set(gca,'XTick',mindelay:part:maxdelay);%x轴只标注几个边界点
    set(gca,'XTickLabel',num2str(round(get(gca,'xTick')'/100)*100,'% .f'));%保留到百位数的精度
    %set(gca,'XTickLabel',{mindelay,mindelay+part,maxdelay-part,maxdelay});
%     hl=legend(H,'Total','Left','Right','Location','Best');%legend需要在Plot后面，否则之后的plot都会自动加上legend
%     set(hl,'Box','off');%关掉legend的边框
%     set(gca,'FontName','Arial','FontSize',14);
end

%plot licking raster
resample=[];%'resample';%trial数量太多的话可以可考虑resample
rasterSize=2;%每个lick点的横向长度，因为分辨率为1ms小鼠不可能舔水这么密集，此值甚至可以到上百
[lickTime{n},maxdelay]=fGetLickTimes(SessionResults,'go');%align={'stim','go'};
%排序
subplot(2,3,[1,4]);%这种设定能画指定数量张图，无violation的1000ms和首次加上300ms的violation
fPlotLickRasterOneSession(lickTime{n},[0 1],'withVio','sort',filename,animal,maxdelay,resample,rasterSize);

%plot curve
[animal_name,dataChoice]=fFindChoice(path,1);
b=dataChoice;
%correctProbe=fCorrectRateProbePerSession(choice);
%获取实验组和对照组的两组数据，并画图 
 p_ctrl=find(dataChoice(1,:,1)==0);%这里有bug，当坐满1000trial时是空值
 if isempty(p_ctrl) 
     p_ctrl=1000;
 end
 choice_ctrl=squeeze(dataChoice(1,1:p_ctrl-1,:));  
[correctProbe_ctrl,toneOct_ctrl,choicep_ctrl,f_ctrl,x_ctrl,yfit_ctrl]=fProbePerformance(choice_ctrl);
subplot(2,3,6);
ctrl_curve=fPlotProbePerformance(fig,'k',correctProbe_ctrl,toneOct_ctrl,choicep_ctrl,x_ctrl,yfit_ctrl,animal_name); 
end

function [strtitle]=fPlotLickRasterOneSession(lickTime,trialType,resultType,sortType,filename,animal,maxdelay,nresample,rasterSize)
%trialType选择trial为左0右1,均选择则用[0 1];ismember(lickTime.trialType,trialType)
%resultType选择result是‘withVio'所有trial，’noVio'不画violation
%sortType选择是否sort，'sort'表示sort，'raw'表示不sort
if strcmp(nresample,'resample')
    nsample=unidrnd(length(lickTime),1,200);%trial 数量太多，抽取其中一部分
    nsample=sort(nsample);
    lickTime=lickTime(nsample);
else
    nsample=1:length(lickTime);
end
if strcmp(sortType,'sort')
    [~,index]=sort([lickTime.delay]);
elseif strcmp(sortType,'raw')
    index=1:1:length(lickTime);%when sort don't need
end
yi=1;%控制lick raster在第几行，如果只画左侧或者右侧的trial， 可以把中间空余的行去掉，向下堆叠
for i=1:length(nsample) %i个trial
    if ismember(lickTime(index(i)).trialType,trialType)
%         line([lickTime(index(i)).stimOn,lickTime(index(i)).stimOff],[i-0.5,i+0.5],'color',[0.9,0.9,0.9],'linewidth',2);%stim, indicated with grey
        line([lickTime(index(i)).stimOn,lickTime(index(i)).stimOn],[yi-0.5,yi+0.5],'color','k','linewidth',2);%stimOn
        hold on;
        line([lickTime(index(i)).stimOff,lickTime(index(i)).stimOff],[yi-0.5,yi+0.5],'color','k','linewidth',2);%stimOff
        hold on;
        line([lickTime(index(i)).go,lickTime(index(i)).go],[yi-0.5,yi+0.5],'color','k','linewidth',2);%go cue
        hold on;
        if strcmp(resultType,'noVio') && (lickTime(index(i)).choice==3)%不画vio的图
           continue;
        else     
            for jl=1:length(lickTime(index(i)).lickleft)
                line([lickTime(index(i)).lickleft(jl),lickTime(index(i)).lickleft(jl)],[yi-0.5,yi+0.5],'color','b','linewidth',rasterSize);%left lick
                hold on;
            end
            for jr=1:length(lickTime(index(i)).lickright)
                line([lickTime(index(i)).lickright(jr),lickTime(index(i)).lickright(jr)],[yi-0.5,yi+0.5],'color','r','linewidth',rasterSize);%right lick
                hold on;
            end
            plot(lickTime(index(i)).answer,yi,'k.');%answer lick
            hold on;
        end
        yi=yi+1;
    end
end
%分析文件名，获得date,animal信息
if  contains(filename,'pyx019') || contains(filename,'pyx046')|| contains(filename,'pyx051')
    pdate=strfind(filename,'_');
    pend=strfind(filename,'.');
    date=filename(pdate(1)+1:pend(end)-1);
    strtitle=strcat(animal,'-',date);
else
    pdate=strfind(filename,'_');
    date=filename(1:pdate(end)-1);
    date(date=='_')='';%将下划线取消掉
    strtitle=strcat(animal,'-',date);
end
set(gca, 'YLim',[0 yi]);
set(gca, 'XLim',[-maxdelay maxdelay-500]);
xlabel('time(ms) from go');
% xlabel('time(ms) from go cue','FontName','Arial','FontSize',14);
ylabel('trials');
% title(strtitle,'FontName','Arial','FontSize',14);
plot([0,0],[0,1000],'k','linewidth',1);
hold on;
% set(gca,'FontName','Arial','FontSize',14);
end

function [curve] = fPlotProbePerformance(fig,color,correctProbe,toneOct,choicep,x,yfit,animal_name)
%x1 = linspace(min(toneOct)-0.1, max(toneOct)+0.1, 100);
    figure(fig);
    set(gcf, 'position', [100 100 1200 300]);%控制fig尺寸
    
    box off;
    points_legend=scatter(2.^choicep(:,1)*20,choicep(:,2)*100,correctProbe(1:size(correctProbe,1)-3,6),color,'filled');
    hold on;
    pointsize=sort(correctProbe(1:size(correctProbe,1)-3,6));
    %pointsizelevel=ceil(pointsize(end-2)/10)*10;%find probe trial numbers,choose the largest and decide legend markersize.
%     pointsizelevel=ceil(pointsize(1)/10)*10;%find probe trial numbers,choose the smallest and decide legend markersize.
%     [lgd,icons,~,~]=legend(points_legend,strcat(num2str(pointsizelevel),' trials'),'Location','Best');
%     set(lgd,'Box','off');
%     icons(2).Children.MarkerSize = sqrt(pointsizelevel);%legend as example how big the marker is
    %errorbar(choicep(:,1),choicep(:,2),correctProbe(:,4),'LineWidth',0);
    %hold on;
    curve=plot(x ,yfit*100 ,color,'LineWidth',2);%标记该曲线，为legend提供名称
    hold on;
%     ylabel('% Choice of High Click Rate Side','FontName','Arial','FontSize',14);
ylabel('% Choice Right');
    xlabel('clicks/s');
%     title(['Psychometric curve'],'FontName','Arial','FontWeight','Bold','FontSize',16);
%     set(gca,'FontName','Arial','FontSize',14);
    set(gca, 'YLim',[0 100]);
    set(gca, 'XLim',[x(1) x(end)]); 
    set(gca,'XScale','log');
    set(gca,'XTick',[20,125],'LineWidth',1);
    set(gca,'XTickLabel',[20,125]);
    set(gca,'YTick',[0,50,100]);
    set(gca,'YTickLabel',[0,50,100]);
    plot([x(1),x(end)],[50,50],'--k','LineWidth',2);
    hold on;   
%     plot([50,50],[0,100],'--k','LineWidth',1);
    hold on;
end
function [correctProbe,toneOct,choicep,f,x,yfit] = fProbePerformance(choice)
    %click=[20 28 40 62 90 125];
    click=unique(choice(:,1));
    if length(click)==2
        warning('Not enough stimuli');
    end
    freq_n=length(click);
    correctProbe=zeros(freq_n+3,6);%每行表示一种频率，每列分别表示频率，正确率，miss率，正确率的std，violation率，correct+error的trial数量之和；最后另外加三行，分别表示end trial,和probe trial, all trials

    pclick=zeros(size(choice,1),freq_n);%用于记录每个click rate的trial位置
    for i=1:freq_n
        pclick(:,i)=double(choice(:,1)==click(i));
    end
    pclick(:,freq_n+1)=pclick(:,1)+pclick(:,freq_n);
    pclick(:,freq_n+2)=sum(pclick(:,2:freq_n-1),2);%对行求和
    pclick(:,freq_n+3)=pclick(:,freq_n+1)+pclick(:,freq_n+2);
    pperformance=zeros(size(choice,1),4);%用于记录各种类型的trial位置，包括correct,error,miss三列
    performance=[1 2 3 4];%分别表示CORRECT,ERROR,MISS,VIOLATION
    for i=1:4
        pperformance(:,i)=double(choice(:,2)==performance(i));
    end
    correctProbe(1:freq_n,1)=click;%最后两行由于是end 和probez，暂时不赋值了
    for i=1:freq_n+3
        do=(pperformance(:,1)+pperformance(:,2)).*pclick(:,i);
        notMiss=(do+pperformance(:,4)).*pclick(:,i);
        correctProbe(i,2)=sum(pperformance(:,1).*pclick(:,i))/sum(do);
        correctProbe(i,3)=sum(pperformance(:,3).*pclick(:,i))/sum(pclick(:,i));
        %计算各个频率的正确率的std
        correctNoMiss=pperformance(:,1).*pclick(:,i);
        correctNoMiss(~notMiss)=nan;
        correctProbe(i,4)=nanstd(correctNoMiss);
        correctProbe(i,5)=sum(pperformance(:,4).*pclick(:,i))/sum(notMiss);
        correctProbe(i,6)=sum(do);
    end
    
    %拟合
    toneOct  = log2(correctProbe(1:freq_n,1)/20);
    choicep =correctProbe(1:freq_n,1:2);
    choicep(:,1)=toneOct ;
    choicep(1:freq_n/2,2)=1-correctProbe(1:freq_n/2,2);
    ffun=fittype('a+b./(1+exp(-k.*(x-c)))','independent','x');
    f =fit(choicep(:,1),choicep(:,2),ffun,'Startpoint',[0,1,1,1]);
    slope=f.b*f.k/4;
    x =toneOct(1)-0.5:0.01:toneOct(end)+0.5;
    yfit =f(x);
    x=2.^x*20;
end
