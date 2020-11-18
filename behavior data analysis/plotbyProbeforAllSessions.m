%   读取.mat文件，绘制population 结果
%   A包含不同delay的curve，cor/miss/vio rate，
%   B比较saline和CNO 条件下的cor/miss/vio rate等
%   C比较不同delay/difficulty条件下saline与CNO之间的差异
%   delay-difficulty interaction可以在各个参数上进行包括cor, RT,miss, vio等，
%   常见bug：运行时确认fFindChoiceMat中click为确定的6个值，否则会出现nan无法拟合；已经解决，但是不排除还会发生bug
clear;
% close all;%close all figure
%  dir='D:\xulab\behavior\CNO experiment\CNO ctrl\infusion\bilateral\';
% dir='D:\xulab\behavior\CNO experiment\CNO ctrl\infusion\unilateral\';
% dir='D:\xulab\behavior\CNO experiment\CNO ctrl\ip\';
dir='D:\xulab\behavior\CNO experiment\infusion\bilateral\';
%  dir='D:\xulab\behavior\CNO experiment\infusion\unilateral\';
% dir='D:\xulab\behavior\CNO experiment\ip\';
group={'saline', 'CNO'};
color={'k','r'}; 
type=  {'single case','population'};
curve=cell(length(group),1);%用于存储图像句柄  
transform_method={'top-down','left-right'};
% transform_method={'top-down'};
% transform_method={'left-right'};
% transform_method={};
if contains(dir,'uni')
    transform_method={'uni'};
end  
%%
%绘制每个session以及所有trial合并后的总的curve
slope=cell(length(transform_method),length(group));%第一维因为slope和旋转方式有关，故存多套
for tn=1:length(transform_method)%两种curve变换方式都进行
    A(tn)=figure;   
    for n=1:length(group)
        [animal_name,dataChoiceResult]=fFindChoiceMat(strcat(dir,group{n}));
        nsample=size(dataChoiceResult,1)-1;
        curve{n}=cell(size(dataChoiceResult,1),1);
        slope{tn,n}=cell(nsample+1,1);%最后一列表示的是总curve的slope
        for i=1:nsample
            [correctProbe,toneOct,choicep,f]=fProbePerformance(dataChoiceResult{i},transform_method{tn});
            curve{n}{i}=fPlotProbePerformance(A(tn),color{n},correctProbe,toneOct,choicep,f,animal_name{end},type{1},transform_method{tn});
            slope{tn,n}{i}=cellfun(@(x) x.b.*x.k/4, f,'UniformOutput', false);
        end
        %population mean
        %%%%%%%%%%%%%%%%%%此处datachoice的数据将不同session汇总以后会丢失rule信息，左右转的情况下无法旋转；考虑重新设计datachoice
        [correctProbe,toneOct,choicep,f]=fProbePerformance(dataChoiceResult{nsample+1},transform_method{tn});
        curve{n}{nsample+1}=fPlotProbePerformance(A(tn),color{n},correctProbe,toneOct,choicep,f,animal_name{end},type{2},transform_method{tn});
        slope{tn,n}{end}=cellfun(@(x) x.b.*x.k/4, f,'UniformOutput', false);
    end
%     for i=1:3%比较saline和CNO curve的slope变化
%         subplot(3,5,5*i); 
%         for k=1:length(slope{tn,1})-1
%             scatter(slope{tn,1}{k}{i},slope{tn,2}{k}{i},50,'k');  
%             hold on;
%         end
%         scatter(slope{tn,1}{end}{i},slope{tn,2}{end}{i},50,'k','filled');
%         hold on;
%         if i==1
%             title('Slop of curve','FontName','Arial','FontSize',14);
%         end
%         if i==3
%             xlabel('Saline','FontName','Arial','FontSize',14);
%         end
%         ylabel('CNO','FontName','Arial','FontSize',14);
%         plot([0,10],[0,10],'k','LineWidth',2);
%         hold on;
%         set(gca,'FontName','Arial','FontSize',14);
%     end
end
  
%add legend to psychometric curve
%figure(A);
%subplot(1,4,4);
%hlegend=legend([curve{1}{end}, curve{2}{end}],'ctrl','ip CNO','Location','Best');%curve{1}是ctrl，curve{2}是exp
%hlegend=legend([curve{1}{end}, curve{2}{end}],group,'Location','Best');
%set(hlegend,'Box','off');
%%
%%二维点图，比较correct rate, violation rate, miss rate, curve slope，reaction time from go cue(for correct trials),区分长短delay(3*5 subplot)
%paraB=cell(5,3,length(group));%第一维代表五种指标cor,vio,miss,RT,cor_diff(endtrial-probe)，第二维代表总、短delay、长delay，第三维代表两种实验条件（saline,CNO)
% if strcmp(transform_method,'uni')
%     %group={'ipsi saline', 'ipsi CNO', 'contra saline', 'contra CNO'};
%     group={'saline', 'CNO', 'saline', 'CNO'};  
% end

%    定义delay组数或者自定义delay分组标准获得不同delay长度的参数paraB
ndelaygroup=2;
delaysetting=[300 900 900 1500];
ndifficultygroup=2;
[paraB,pparaB,paraBdiff,pparaBdiff,delaybylength]=fGetpara(dir, transform_method, ndelaygroup ,group, ndifficultygroup);
B=figure;%二维点图，比较correct rate, violation rate, miss rate, curve slope，reaction time from go cue(for correct trials),区分长短delay，或者难度(3*5 subplot)
C=figure;%画CNO-saline的差值，比较不同情况下的区别
curveB=fPlotOverallPerformance(B,C,paraB,pparaB,paraBdiff,pparaBdiff,delaybylength,ndelaygroup);
%%
%比较不同参数的correlation，3*6图，第1-3行总trial，长短delay，难易；1-6列：cor-vio,cor-miss,cor-RT,vio-miss,vio-RT,miss-RT
colname={'Correct','Violation','Miss','RT'};
titlestr={'saline','CNO','CNO-saline'};
data=cell(1,3);%分别放saline,CNO,saline-CNO数据
for ifig=1:2%分别画saline,cno
    data{ifig}=cell(1,5);
    for i=1:5
        data{ifig}{i}=zeros(length(paraB{1,1,1,1}),4);
    end
    for icol=1:4
        data{ifig}{1}(:,icol)=paraB{icol,1,ifig,1};%all trial
        data{ifig}{2}(:,icol)=paraB{icol,2,ifig,1};%shortest delay
        data{ifig}{3}(:,icol)=paraB{icol,end,ifig,1};%longest delay
        data{ifig}{4}(:,icol)=paraB{icol,1,ifig,2};%easy trials
        data{ifig}{5}(:,icol)=paraB{icol,1,ifig,3};%hard trials
    end
end
data{3}=cell(1,5);
for i=1:5
    data{3}{i}=zeros(length(paraBdiff{1,1,1,1}),4);
end
for icol=1:4
    data{3}{1}(:,icol)=paraBdiff{icol,1,1,1};
    data{3}{2}(:,icol)=paraBdiff{icol,2,1,1};
    data{3}{3}(:,icol)=paraBdiff{icol,end,1,1};
    data{3}{4}(:,icol)=paraBdiff{icol,1,1,2};
    data{3}{5}(:,icol)=paraBdiff{icol,1,1,3};
end
for ifig=1:3
     fPlotCorr(data{ifig},titlestr{ifig},colname);
end
%%
% plot difficulty-delay interaction for each parameters, e.g. correct rate,
% miss rate, violation rate, reaction time, etc.
% can assign whether scatter plot or dot-point plot;which parameters to
% show
paraname='cor rate'; %this variable need to be assigned according to different problems
rowname_interaction={'all difficulty','easy','hard'};
colname_interaction={'all delay','short delay','long delay'};
parafield={'cor rate','violation rate','miss rate', 'RT'};
paraindex=cellfun(@(x) strcmp(x,paraname),parafield);
data_interaction=cell(2,ndifficultygroup+1,ndelaygroup+1);
for i=1:2
    for j=1:ndifficultygroup+1
        for k=1:ndelaygroup+1
            data_interaction{i,j,k}=paraB{paraindex,k,i,j}';
        end
    end
end
fig=fPlot2FactorInteraction(data_interaction,'3D',rowname_interaction,colname_interaction);
%%
%调用的函数
function [curve] = fPlotProbePerformance(fig,color_population,correctProbe,toneOct,choicep,f,animal_name,type,transform_method)
%x1 = linspace(min(toneOct)-0.1, max(toneOct)+0.1, 100);
    colormode={[0.7,0.7,0.7],[1,0.7,0.7]};
    colormodeabbr={[0,0,0],[1,0,0]};
    colorp=strfind('kr',color_population);
    if strcmp(type,'single case')
        color=colormode{colorp};
    else
        color=colormodeabbr{colorp};
    end
    figure(fig);
    set(gcf, 'position', [0 0 1500 900]);%控制fig尺寸
    for i=1:3
        subplot(3,5,5*i-4);
        if strcmp(type,'population')
    %         scatter(2.^toneOct*20,correctProbe(1:6,2)*100,50,color_population);
            scatter([2.^(toneOct(end)+0.4)*20,2.^(toneOct(end)+0.7)*20,2.^(toneOct(end)+1)*20],correctProbe(i,7:9,2)*100,50,color_population,'filled');
    %         scatter(1,correctProbe(9,2)*100,50,color_population,'filled');
            hold on;    
            plot(2.^toneOct*20,correctProbe(i,1:6,2)*100,'Color',color,'LineWidth',2);
            hold on;
        else
            plot(2.^toneOct*20,correctProbe(i,1:6,2)*100,'Color',color,'LineWidth',0.5);
            hold on;
        end
        ylabelstr={'All','Short delay','Long delay'};
        ylabel(strcat(ylabelstr{i},' trials'),'FontName','Arial','FontSize',14);
        if i==1
            title('% Correct','FontName','Arial','FontSize',14);
        end
        set(gca,'XTick',1,'LineWidth',1);
        set(gca,'XTickLabel',' ');
    %     title(['Performance of ', animal_name],'FontName','Arial','FontWeight','Bold','FontSize',16);
    %     plot([2.^(toneOct(1)-0.5)*20, 2.^(toneOct(end)+0.5)*20],[80,80],'Color',color,'LineWidth',1);
    %     hold on;
        set(gca, 'YLim',[0 100]);
        set(gca,'FontName','Arial','FontSize',14);
        set(gca,'XScale','log');
        set(gca, 'XLim',[2.^(toneOct(1)-0.5)*20, 2.^(toneOct(end)+1.5)*20]);
        if i==3
            if strcmp(transform_method,'top-down')
                xlabel('clicks/s','FontName','Arial','FontSize',14);
                set(gca,'XTick',[20,125],'LineWidth',1);
                set(gca,'XTickLabel',[20,125]);
            elseif strcmp(transform_method,'left-right')
                xlabel('stimulus side','FontName','Arial','FontSize',14);
                set(gca,'XTick',[20,125],'LineWidth',1);
                set(gca,'XTickLabel',{'left','right'});
            elseif strcmp(transform_method,'uni')
                xlabel('stimulus side','FontName','Arial','FontSize',14);
                set(gca,'XTick',[20,125],'LineWidth',1);
                set(gca,'XTickLabel',{'ipsi','contra'});
            end
        end

        subplot(3,5,5*i-3);

        if strcmp(type,'population')
            plot(2.^toneOct*20 ,correctProbe(i,1:6,3)*100,'Color',color,'LineWidth',2);
            hold on;
        else
            plot(2.^toneOct*20 ,correctProbe(i,1:6,3)*100,'Color',color,'LineWidth',0.5);
            hold on;
        end
%         ylabel('% Miss','FontName','Arial','FontSize',14);
        if i==1
            title('% Miss','FontName','Arial','FontSize',14);
        end
    %     title(['Miss rate of ' ,animal_name],'FontName','Arial','FontWeight','Bold','FontSize',16);
        hold on;
        set(gca,'XScale','log');
        set(gca, 'XLim',[2.^(toneOct(1)-0.5)*20, 2.^(toneOct(end)+0.5)*20]);
        set(gca,'XTick',[20,125],'LineWidth',1);
        set(gca,'XTickLabel',' ');
        if i==3
            if strcmp(transform_method,'top-down')
                xlabel('clicks/s','FontName','Arial','FontSize',14);
                set(gca,'XTick',[20,125],'LineWidth',1);
                set(gca,'XTickLabel',[20,125]);
            elseif strcmp(transform_method,'left-right')
                xlabel('stimulus side','FontName','Arial','FontSize',14);
                set(gca,'XTick',[20,125],'LineWidth',1);
                set(gca,'XTickLabel',{'left','right'});
            elseif strcmp(transform_method,'uni')
                xlabel('stimulus side','FontName','Arial','FontSize',14);
                set(gca,'XTick',[20,125],'LineWidth',1);
                set(gca,'XTickLabel',{'ipsi','contra'});
            end
        end
        set(gca,'FontName','Arial','FontSize',14);

        subplot(3,5,5*i-2);
        if strcmp(type,'population')
            plot(2.^toneOct*20 ,correctProbe(i,1:6,5)*100,'Color',color,'LineWidth',2);
            hold on;
        else
            plot(2.^toneOct*20 ,correctProbe(i,1:6,5)*100,'Color',color,'LineWidth',0.5);
            hold on;
        end
%         ylabel('% Violation','FontName','Arial','FontSize',14);
        if i==1
            title('% Violation','FontName','Arial','FontSize',14);
        end
        if i==3
            xlabel('clicks/s','FontName','Arial','FontSize',14);
        end
    %     title(['Violation rate of ' ,animal_name],'FontName','Arial','FontWeight','Bold','FontSize',16);
        set(gca,'XScale','log');
        set(gca, 'XLim',[2.^(toneOct(1)-0.5)*20, 2.^(toneOct(end)+0.5)*20]);
        set(gca,'XTick',[20,125],'LineWidth',1);
        set(gca,'XTickLabel',' ');
        if i==3
            if strcmp(transform_method,'top-down')
                xlabel('clicks/s','FontName','Arial','FontSize',14);
                set(gca,'XTick',[20,125],'LineWidth',1);
                set(gca,'XTickLabel',[20,125]);
            elseif strcmp(transform_method,'left-right')
                xlabel('stimulus side','FontName','Arial','FontSize',14);
                set(gca,'XTick',[20,125],'LineWidth',1);
                set(gca,'XTickLabel',{'left','right'});
            elseif strcmp(transform_method,'uni')
                xlabel('stimulus side','FontName','Arial','FontSize',14);
                set(gca,'XTick',[20,125],'LineWidth',1);
                set(gca,'XTickLabel',{'ipsi','contra'});
            end
        end
        set(gca,'FontName','Arial','FontSize',14);

        subplot(3,5,5*i-1);
        if strcmp(type,'population')
            plot(2.^toneOct*20 ,correctProbe(i,1:6,6),'Color',color,'LineWidth',2);
            hold on;
        else
            plot(2.^toneOct*20 ,correctProbe(i,1:6,6),'Color',color,'LineWidth',0.5);
            hold on;
        end
%         ylabel('% Violation','FontName','Arial','FontSize',14);
        if i==1
            title('Reaction Time(ms)','FontName','Arial','FontSize',14);
        end
        if i==3
            xlabel('clicks/s','FontName','Arial','FontSize',14);
        end
    %     title(['Violation rate of ' ,animal_name],'FontName','Arial','FontWeight','Bold','FontSize',16);
        set(gca,'XScale','log');
        set(gca, 'XLim',[2.^(toneOct(1)-0.5)*20, 2.^(toneOct(end)+0.5)*20]);
        set(gca,'XTick',[20,125],'LineWidth',1);
        set(gca,'XTickLabel',' ');
        if i==3
            if strcmp(transform_method,'top-down')
                xlabel('clicks/s','FontName','Arial','FontSize',14);
                set(gca,'XTick',[20,125],'LineWidth',1);
                set(gca,'XTickLabel',[20,125]);
            elseif strcmp(transform_method,'left-right')
                xlabel('stimulus side','FontName','Arial','FontSize',14);
                set(gca,'XTick',[20,125],'LineWidth',1);
                set(gca,'XTickLabel',{'left','right'});
            elseif strcmp(transform_method,'uni')
                xlabel('stimulus side','FontName','Arial','FontSize',14);
                set(gca,'XTick',[20,125],'LineWidth',1);
                set(gca,'XTickLabel',{'ipsi','contra'});
            end
        end
        set(gca,'FontName','Arial','FontSize',14);
        
        subplot(3,5,5*i);
        x =toneOct(1)-0.5:0.01:toneOct(end)+0.5;
        yfit =f{i}(x);
        x=2.^x*20;
        box off;
        if  strcmp(type,'population')
            scatter(2.^choicep{i}(:,1)*20,choicep{i}(:,2)*100,50,color_population,'filled');
            hold on;
            curve=plot(x ,yfit*100 ,'Color',color,'LineWidth',2);%标记该曲线，为legend提供名称
            hold on;
        else
            curve=plot(x ,yfit*100 ,'Color',color,'LineWidth',0.5);%标记该曲线，为legend提供名称
            hold on;
        end
        %errorbar(choicep(:,1),choicep(:,2),correctProbe(:,4),'LineWidth',0);
        %hold on;
        set(gca, 'XLim',[2.^(toneOct(1)-0.5)*20, 2.^(toneOct(end)+0.5)*20]); 
        set(gca,'XScale','log');
    %     title(['Psychometric curve'],'FontName','Arial','FontWeight','Bold','FontSize',16);
        set(gca, 'YLim',[0 100]);
        set(gca,'YTick',[0,50,100]);
        set(gca,'YTickLabel',[0,50,100]);
        plot([x(1),x(end)],[50,50],'--k','LineWidth',2);
        hold on;
    %    plot([50,50],[0,100],'--k','LineWidth',1);
    %    hold on;
        if strcmp(transform_method,'top-down')
            ylabel('% Choice of High','FontName','Arial','FontSize',14);
%             if i==1
%                 title('% Choice of High Click Rate Side','FontName','Arial','FontSize',14);
%             end
%             if i==3
                xlabel('clicks/s','FontName','Arial','FontSize',14);
%             end
            set(gca,'XTick',[20,125],'LineWidth',1);
            set(gca,'XTickLabel',[20,125]);
        elseif strcmp(transform_method,'left-right')
            ylabel('% Choice of Right','FontName','Arial','FontSize',14);
%             if i==1
%                 title('% Choice of Right','FontName','Arial','FontSize',14);
%             end
%             if i==3
%                 xlabel('stimulus side','FontName','Arial','FontSize',14);
%             end
            set(gca,'XTick',[20,125],'LineWidth',1);
            set(gca,'XTickLabel',{'left','right'});
        elseif strcmp(transform_method,'uni')
            if i==1
                title('% Choice of Contra Side','FontName','Arial','FontSize',14);
            end
            if i==3
                xlabel('stimulus side','FontName','Arial','FontSize',14);
            end
            set(gca,'XTick',[20,125],'LineWidth',1);
            set(gca,'XTickLabel',{'ipsi','contra'});
        end
        set(gca,'FontName','Arial','FontSize',14);
    end
end
function [correctProbe,toneOct,choicep,f] = fProbePerformance(choice,transform_method)
     %delay先确定好
    maxdelay=max(choice(:,4));
    mindelay=min(choice(:,4));
    part=(maxdelay-mindelay)/3;%将delay分成哪几种长度
    shortdelay=(choice(:,4)<=(mindelay+part));
    longdelay=(choice(:,4)>(maxdelay-part));
    delay=[ones(size(choice,1),1),shortdelay,longdelay];%汇总成二维数组，每一列代表不同长度的delay,% logical 相加就变成double
    correctProbe=zeros(size(delay,2),9,7);%第一维区分delay长度，之后第二维表示一种频率，第三维分别表示频率(如果是左右调换的话就是声音对应的左右)，正确率，miss率，正确率的std，violation;RT,LCI率；第二维最后另外加三行，分别表示end trial,和probe trial, all trials
    click=[20 28 40 62 90 125];
    if strcmp(transform_method,'left-right')
        choice(:,1)=choice(:,3);%左右调换的话就用第六列的声音对应的左右替换第一列频率
    elseif strcmp(transform_method,'top-down')

    elseif strcmp(transform_method,'uni')
        choice(:,1)=choice(:,3);%左右调换的话就用第六列的声音对应的左右替换第一列频率
        protate=logical((choice(:,6)==2)+(choice(:,6)==3));%这些trials需要旋转，也就是上下，左右各颠倒一次
        choice(protate,1)=7-choice(protate,1);      
    end
    pclick=zeros(size(choice,1),9);%用于记录每个click rate的trial位置
    for i=1:6
        pclick(:,i)=double(choice(:,1)==click(i));
    end
    pclick(:,7)=pclick(:,1)+pclick(:,6);
    pclick(:,8)=pclick(:,2)+pclick(:,3)+pclick(:,4)+pclick(:,5);
    pclick(:,9)=pclick(:,7)+pclick(:,8);
    pperformance=zeros(size(choice,1),4);%用于记录各种类型的trial位置，包括correct,error,miss，violation 4列
    performance=[1 2 3 4];%分别表示CORRECT,ERROR,MISS,VIOLATION ，right(注意action choice 中右为1)
    for i=1:4
        pperformance(:,i)=double(choice(:,2)==performance(i));
    end
    choicep=cell(size(delay,2),1);
    fitdata=cell(size(delay,2),1);
    f=cell(size(delay,2),1);
    for j=1:size(delay,2)%区分不同delay
        correctProbe(j,1:6,1)=click;%最后两行由于是end 和probez，暂时不赋值了
        for i=1:9
            do=pperformance(:,1)+pperformance(:,2); 
            notMiss=do+pperformance(:,4);
            correctProbe(j,i,2)=sum(pperformance(:,1).*pclick(:,i).*delay(:,j))/sum(do.*pclick(:,i).*delay(:,j));
            correctProbe(j,i,3)=sum(pperformance(:,3).*pclick(:,i).*delay(:,j))/sum(pclick(:,i).*delay(:,j));
            %计算各个频率的正确率的std
            correctNoMiss=pperformance(:,1).*pclick(:,i).*delay(:,j);
            errorNoMiss=pperformance(:,2).*pclick(:,i).*delay(:,j);
            pcor=sum(correctNoMiss)/sum(correctNoMiss+ errorNoMiss);
            correctProbe(j,i,4)=sqrt(sum(correctNoMiss+ errorNoMiss)*pcor*(1-pcor)/(sum(correctNoMiss+ errorNoMiss)-1));
%             correctNoMiss(~notMiss)=nan;
%             correctProbe(j,i,4)=nanstd(correctNoMiss);
            correctProbe(j,i,5)=sum(pperformance(:,4).*pclick(:,i).*delay(:,j))/sum(notMiss.*pclick(:,i).*delay(:,j));%%%此处有误，population比任意个都大，或者偏小     
            correctProbe(j,i,6)=sum(choice(:,5).*pperformance(:,1).*pclick(:,i).*delay(:,j))/sum(pperformance(:,1).*pclick(:,i).*delay(:,j));%RT only include correct trials
            correctProbe(j,i,7)=sum(choice(:,5).*pperformance(:,1).*pclick(:,i).*delay(:,j))/sum(pperformance(:,1).*pclick(:,i).*delay(:,j));%LCI also only include correct trials
%             correctProbe(j,i,7)=sum(choice(:,5).*pclick(:,i).*delay(:,j))/sum(pclick(:,i).*delay(:,j));%LCI include all trials(cor,err,miss,vio)
        end

        %计算choicep，表示选择某个方向的比例，画在curve处，为散点
        %计算toneOct，表示绘图所用的横坐标
        %为了方便两种情况下标度的一致性，先预设为刺激相同尺度
        if strcmp(transform_method,'left-right')
            correctProbe(j,1:6,1)=[20;28;40;62;90;125];
        end
        toneOct  = log2(correctProbe(j,1:6,1)/correctProbe(j,1,1));
        choicep{j} =squeeze(correctProbe(j,1:6,1:2));%6*2矩阵
        choicep{j}(:,1)=toneOct;
        choicep{j}(1:3,2)=1-correctProbe(j,1:3,2);
        
        %拟合
        ffun=fittype('a+b./(1+exp(-k.*(x-c)))','independent','x');
        %方法一，采用average过后总的正确率计算，即每个stimulus一个点
%         f{j} =fit(choicep{j}(:,1),choicep{j}(:,2),ffun,'Startpoint',[0,1,1,1]);
        %方法二，采用单个trial来拟合，一个trial一个点，自变量为stimulus，因变量0-1；
        %这种方法会对end point更加倚重 
        fitdata{j}=choice(logical((choice(:,2)<3).*delay(:,j)),1:2);%第一列1-6代表六种刺激，第二列分别表示CORRECT,ERROR,MISS,VIOLATION，只保留correct,error，即1,2
        if strcmp(transform_method,'left-right')
            fitdata{j}((fitdata{j}(:,1)<50),2)=fitdata{j}((fitdata{j}(:,1)<50),2)-1;%前三个stimulus，应当选左，转化为0（对，左）-1（错，右）
            fitdata{j}((fitdata{j}(:,1)>50),2)=2-fitdata{j}((fitdata{j}(:,1)>50),2);%后三个stimulus，应当选右，转化为0（错，左）-1（对，右）    
        elseif strcmp(transform_method,'top-down')
            fitdata{j}((fitdata{j}(:,1)<50),2)=fitdata{j}((fitdata{j}(:,1)<50),2)-1;%前三个stimulus，应当选低频侧，转化为0（对，低）-1（错，高）
            fitdata{j}((fitdata{j}(:,1)>50),2)=2-fitdata{j}((fitdata{j}(:,1)>50),2);%后三个stimulus，应当选高频侧，转化为0（错，低）-1（对，高）
        end
        fitdata{j}(:,1)=log2(fitdata{j}(:,1)/min(fitdata{j}(:,1)));%stimulus先转化为对数
        f{j} =fit(fitdata{j}(:,1),fitdata{j}(:,2),ffun,'Startpoint',[0,1,1,1]);
    end
end
function [curve] = fPlotOverallPerformance(fig1,fig2,paraB,pparaB,paraBdiff,pparaBdiff,delaybylength,ndelaygroup)
%    分别在fig1和fig2画二维点图和分组比较差值的图
    figure(fig1);
    set(gcf, 'position', [0 0 1500 810]);%控制fig尺寸
%     ylabelstr={' all trials',' short delay trials',' long delay trials'};
    ylabelstr={'CNO','CNO','CNO'};
    %xlabelstr={};
    titlestr={'Correct rate','Violation rate','Miss rate','Reaction time','Licking Consistency'};
    for i=1:3
        for j=1:5
            subplot(3,5,5*i+j-5);%五列分别画正确率，violation ,miss,RT,LCI

            if j==4%画对角线,RT
                plot([0,1500],[0,1500],'k','LineWidth',2);
                hold on;
            elseif j==5 %LCI
                plot([0,1],[0,1],'k','LineWidth',2);
                hold on;
            else
                plot([0,100],[0,100],'k','LineWidth',2);
                hold on;
            end
            ytext=get(gca,'Ylim');%获取坐标范围，在特定处显示text
            if size(paraB, 3)==2%bilateral
                if i==1%第一行画所有trial的总体情况
                    curve=scatter(paraB{j,i,1,1},paraB{j,i,2,1},50,'k','filled');
                    hold on;                    
                    textstr=plabel(pparaB(j,i,1,1));%自动根据p值选择符号
                elseif i==2%第二行画区分delay长度
%                     if j~=5%第五列改动为区分长短delay和hard一共四种情况
                    curve_short=scatter(paraB{j,2,1,end},paraB{j,2,2,end},50,'k','filled');%最短delay
                    hold on;
                    curve_long=scatter(paraB{j,end,1,end},paraB{j,end,2,end},50,'b','filled');%最长delay
                    hold on;
%                     textstr={strcat('short:',plabelsymbol(pparaB(j,2,1,1)));strcat('long: ',plabelsymbol(pparaB(j,3,1,1)));strcat('short vs long: ',plabelsymbol(pparaBdiff(i-1,j)))};
                   textstr={strcat('short:',plabelsymbol(pparaB(j,2,1,end)));strcat('long: ',plabelsymbol(pparaB(j,end,1,end)))};
%                     else
%                         curve_short_difficult=scatter(paraB{1,2,1,3},paraB{1,2,2,3},30,'k');
%                         hold on;
%                         curve_long_difficult=scatter(paraB{1,3,1,3},paraB{1,3,2,3},30,'b');
%                         hold on;
%                         curve_short_easy=scatter(paraB{1,2,1,2},paraB{1,2,2,2},30,'k','filled');
%                         hold on;
%                         curve_long_easy=scatter(paraB{1,3,1,2},paraB{1,3,2,2},30,'b','filled');
%                         hold on;
% 
%                         textstr={};%strcat('short p=',num2str(pparaB(j,2,1,1)));strcat('long p=',num2str(pparaB(j,3,1,1)));strcat('long p=',num2str(pparaB(j,3,1,1)));strcat('long p=',num2str(pparaB(j,3,1,1)))};
%                     end
                    if j==5%仅需在之后一列画legend
                        legend([curve_short curve_long],strcat('shortest delay(',num2str(delaybylength(1)),'-',num2str(delaybylength(2)),'ms)'),strcat('longest delay(',num2str(delaybylength(end-1)),'-',num2str(delaybylength(end)),'ms)'),'Location','Best');
%                           legend([curve_short_easy curve_long_easy curve_short_difficult curve_long_difficult],strcat('short delay(300-900ms) easy p=',num2str(pparaB(1,2,1,2))),strcat('long delay(900-1500ms) easy p=',num2str(pparaB(1,3,1,2))),strcat('short delay(300-900ms) hard p=',num2str(pparaB(1,2,1,3))),strcat('long delay(900-1500ms) hard p=',num2str(pparaB(1,3,1,3))),'Location','Best');
                    end
                else%第三行画区分难度end/probe
                    curve_end=scatter(paraB{j,1,1,2},paraB{j,1,2,2},50,'k','filled');%end
                    hold on;
                    curve_probe=scatter(paraB{j,1,1,3},paraB{j,1,2,3},50,'k');%probe
                    hold on;
%                     textstr={strcat('easy:',plabelsymbol(pparaB(j,1,1,2)));strcat('hard:',plabelsymbol(pparaB(j,1,1,3)));strcat('easy vs hard:',plabelsymbol(pparaBdiff(i-1,j)))};
                    textstr={strcat('easy:',plabelsymbol(pparaB(j,1,1,2)));strcat('hard:',plabelsymbol(pparaB(j,1,1,end)))};
                    if j==5%仅需在之后一列画legend
                        legend([curve_end curve_probe],'easy','hard','Location','Best');
                    end
                end
                text(0,ytext(end),textstr,'FontSize',12);
            else %size(paraB, 3)==4；unilateral data
                if i==1%第一行画所有trial的总体情况
                    curve_ipsi=scatter(paraB{j,i,1,1},paraB{j,i,2,1},50,'k','filled');
                    hold on;
                    curve_contra=scatter(paraB{j,i,3,1},paraB{j,i,4,1},50,'k','d');
                    hold on;
                    textstr={strcat('ipsi',plabelsymbol(pparaB(j,i,1,1)));strcat('contra',plabelsymbol(pparaB(j,i,2,1)))};
                    if j==5%仅需在之后一列画legend
                        legend([curve_ipsi curve_contra],'ipsi','contra','Location','Best');
                    end
                elseif i==2%第二行画区分delay长度
                    curve_short_ipsi=scatter(paraB{j,2,1,1},paraB{j,2,2,1},50,'k','filled');%短delay
                    hold on;
                    curve_long_ipsi=scatter(paraB{j,3,1,1},paraB{j,3,2,1},50,'b','filled');%长delay
                    hold on;                 
                    curve_short_contra=scatter(paraB{j,2,1,1},paraB{j,2,2,1},50,'k','filled','d');%短delay,contra
                    hold on;
                    curve_long_contra=scatter(paraB{j,3,1,1},paraB{j,3,2,1},50,'b','filled','d');%长delay,contra
                    hold on;
                    textstr={strcat('ipsi short p=',plabelsymbol(pparaB(j,2,1,1)));strcat('ipsi long p=',plabelsymbol(pparaB(j,3,1,1)));strcat('contra short p=',plabelsymbol(pparaB(j,2,2,1)));strcat('contra long p=',plabelsymbol(pparaB(j,3,2,1)))};
                    if j==5%仅需在之后一列画legend
                        legend([curve_short_ipsi curve_long_ipsi curve_short_contra curve_long_contra],'ipsi short delay','ipsi long delay','contra short delay','contra long delay','Location','Best');
                    end
                else%第三行画区分难度end/probe
                    curve_end_ipsi=scatter(paraB{j,1,1,2},paraB{j,1,2,2},50,'k','filled');%end
                    hold on;
                    curve_probe_ipsi=scatter(paraB{j,1,1,3},paraB{j,1,2,3},50,'k');%probe
                    hold on;
                    curve_end_contra=scatter(paraB{j,1,1,2},paraB{j,1,2,2},50,'k','filled','d');%end,contra
                    hold on;
                    curve_probe_contra=scatter(paraB{j,1,1,3},paraB{j,1,2,3},50,'k','d');%probe,contra
                    hold on;
                    textstr={strcat('ipsi easy p=',plabelsymbol(pparaB(j,1,1,2)));strcat('ipsi difficult p=',plabelsymbol(pparaB(j,1,1,3)));strcat('contra easy p=',plabelsymbol(pparaB(j,1,2,2)));strcat('contra difficult p=',plabelsymbol(pparaB(j,1,2,3)))};
                    if j==5%仅需在之后一列画legend
                        legend([curve_end_ipsi curve_probe_ipsi curve_end_contra curve_probe_contra],'ipsi easy','ipsi difficult','contra easy','contra difficult','Location','Best');
                    end
                end
                text(0,ytext(end),textstr,'FontSize',12);
            end

            % set(gca, 'YLim',[0 2000]);
            % set(gca,'XLim',[0,2000]);
%             if j==1
                ylabel(ylabelstr{i},'FontName','Arial','FontSize',14);
%             end
%             if i==3
                %xlabel(strcat('Saline',xlabelstr{i}),'FontName','Arial','FontSize',14);
                if j==4
                    xlabel('Saline(ms)','FontName','Arial','FontSize',14);
                else
                    xlabel('Saline(%)','FontName','Arial','FontSize',14);
                end
%             end
            if i==1
                title(titlestr{j},'FontName','Arial','FontSize',14);
            end
            hold on;
            set(gca,'FontName','Arial','FontSize',14);
            box off;%去掉边框
        end
    end
    figure(fig2);%2*4张图，分别存correct rate, violation rate, miss rate, reaction time；第一行总结果，  第二行区分delay，不只是长短，而是分成200ms一组的6组,第三行区分difficulty
    set(gcf, 'position', [0 0 800 800]);%控制fig尺寸
    ylabelstr={'Correct rate change(%)';'Violation rate change(%)';'Miss rate change(%)';'Reaction time change(ms)'};
    titlestr={'Correct';'Violation';'Miss';'Reaction time'};
    a=suptitle('CNO-saline');
    set(a,'FontSize',14);
    %首行合并各种指标散点图于一张大图
    subplot(3,4,1:3);%首行合并
    pararaw=zeros(length(paraB{i,1,1,1}),4);%第一维表示样本数，值表示CNO-saline情况
    for i=1:4
        pararaw(:,i)=paraB{i,1,2,1}-paraB{i,1,1,1};
    end
    for k=1:3
        yyaxis left;
        scatter(ones(size(pararaw,1),1)*k,pararaw(:,k),30,'filled');
        hold on;
        [~,p]=ttest(pararaw(:,k));
        ytext=get(gca,'Ylim');
        text(k,ytext(end),plabelsymbol(p),'FontName','Arial','FontSize',14);
    end
    ylabel('\Delta (%)');
    yyaxis right;
    scatter(ones(size(pararaw,1),1)*4,pararaw(:,4),30,'filled');
    hold on;
%     for i=1:size(pararaw,1)%用线连接相同session的参数
%         plot(pararaw(i,:));
%         hold on;
%     end
    [~,p]=ttest(pararaw(:,4));
    ytext=get(gca,'Ylim');
    text(4,ytext(end),plabelsymbol(p),'FontName','Arial','FontSize',14);
    plot([0,size(pararaw,2)+1],[0,0],'--k','LineWidth',0.5);
    hold on;
%     set(gca,'Ylim',[-600,600]);
    ylabel('\Delta (ms)');
    set(gca,'XTick',[1 2 3 4]);
    set(gca,'XTickLabel',titlestr);
    set(gca,'FontName','Arial','FontSize',14);
    box off;
       
    for i=1:4 %控制四列       
       subplot(3,4,4+i);%第二行,delay
       parabydelay=zeros(length(paraBdiff{i,1,1,1}),ndelaygroup);
       for j=2:size(paraBdiff,2)
          parabydelay(:,j-1)=paraBdiff{i,j,1,end}';
       end
       for j=1:size(parabydelay,1)
           plot(parabydelay(j,:),'k','LineWidth',0.5);
           hold on;
       end
       for k=1:size(parabydelay,2)
           scatter(ones(size(parabydelay,1),1)*k,parabydelay(:,k),10,'k','filled');
           hold on;
       end
       plot([0,size(parabydelay,2)+1],[0,0],'--k','LineWidth',0.5);
       hold on;
       set(gca, 'XLim',[0 size(parabydelay,2)+1]);
%        set(gca, 'YLim',[0 100]);
       set(gca,'XTick',[1,2]);
%        set(gca,'XTickLabel',['1';'2']);
       set(gca,'XTickLabel',{'short','long'});
       xlabel('delay');
       ylabel(ylabelstr{i});
%        title(titlestr{i});
       [~,p]=ttest(parabydelay(:,1),parabydelay(:,end));
       xtext=get(gca,'Xlim');
       ytext=get(gca,'Ylim');
       text(xtext(end)/2,ytext(end),plabelsymbol(p),'FontName','Arial','FontSize',14);
       set(gca,'FontName','Arial','FontSize',14);
       box off;
       
       subplot(3,4,8+i);%第三行,difficulty
       easy=paraBdiff{i,1,1,2};
       hard=paraBdiff{i,1,1,3};
       parabydifficulty=[easy' hard'];
       for j=1:size(parabydifficulty,1)
           plot(parabydifficulty(j,:),'k','LineWidth',0.5);
           hold on;
       end
       for k=1:size(parabydifficulty,2)
           scatter(ones(size(parabydifficulty,1),1)*k,parabydifficulty(:,k),10,'k','filled');
       end
       plot([0,size(parabydifficulty,2)+1],[0,0],'--k','LineWidth',0.5);
       hold on;
       set(gca, 'XLim',[0 size(parabydifficulty,2)+1]);
%        set(gca, 'YLim',[0 100]);
       set(gca,'XTick',[1,2]);
%        set(gca,'XTickLabel',['1';'2']);
       set(gca,'XTickLabel',{'easy','hard'});
       xlabel('difficulty');
       ylabel(ylabelstr{i});
%        title(titlestr{i});
       [~,p]=ttest(parabydifficulty(:,1),parabydifficulty(:,2));
       xtext=get(gca,'Xlim');
       ytext=get(gca,'Ylim');
       text(xtext(end)/2,ytext(end),plabelsymbol(p),'FontName','Arial','FontSize',14);
       set(gca,'FontName','Arial','FontSize',14);
       box off;
    end
end
function [fig]=fPlotCorr(data,title,colname)%dara为1*5cell，分别包含总，长短delay和难易trial五类数据；data{i}为n*m矩阵，m种指标，n个重复画出两两指标的相关关系,colname表示各个字段的名称，cell,长度为m
npara=size(data{1},2);
if length(colname)~=npara
    disp('data size and colname length not match');
    return;
end
nfigcol=nchoosek(npara,2);
fig=figure;
set(gcf, 'position', [0 0 1800 900]);%控制fig尺寸
suptitle(title);
i=1;
for j=1:npara-1
    for k=j+1:npara
        subplot(3,nfigcol,i);
        scatter(data{1}(:,j),data{1}(:,k),20,'k','filled');
        hold on;
        lsline;%Add least-squares line to scatter plot
        xlabel(colname{j},'FontName','Arial','FontSize',14);
        ylabel(colname{k},'FontName','Arial','FontSize',14);
        [r,p]=corr(data{1}(:,j),data{1}(:,k));
        r=round(r,2);%保留两位即可
        text(max(data{1}(:,j)),max(data{1}(:,k)),{strcat('RHO=',num2str(r));plabelsymbol(p)},'FontSize',10);
        set(gca,'FontName','Arial','FontSize',14);
        
        subplot(3,nfigcol,i+nfigcol);
        curveshort=scatter(data{2}(:,j),data{2}(:,k),20,'k','filled');%short delay
        hold on;
        curvelong=scatter(data{3}(:,j),data{3}(:,k),20,'b','filled');%long delay
        hold on;
        lsline;%Add least-squares line to scatter plot
        if i==nfigcol
            legend([curveshort curvelong],'short delay','longdelay');
        end
        xlabel(colname{j},'FontName','Arial','FontSize',14);
        ylabel(colname{k},'FontName','Arial','FontSize',14);
        [r2,p2]=corr(data{2}(:,j),data{2}(:,k));
        r2=round(r2,2);%保留两位即可
        [r3,p3]=corr(data{3}(:,j),data{3}(:,k));
        r3=round(r3,2);%保留两位即可
        text(max(data{2}(:,j)),max(data{2}(:,k)),{strcat('short RHO=',num2str(r2));plabelsymbol(p2);strcat('long RHO=',num2str(r3));plabelsymbol(p3)},'FontSize',10);
        set(gca,'FontName','Arial','FontSize',14);
        
        subplot(3,nfigcol,i+2*nfigcol);
        curveeasy=scatter(data{4}(:,j),data{4}(:,k),20,'k','filled');%easy
        hold on;
        curvehard=scatter(data{5}(:,j),data{5}(:,k),20,'k');%hard
        hold on;        
        lsline;%Add least-squares line to scatter plot
        if i==nfigcol
            legend([curveeasy curvehard],'easy','hard');
        end
        xlabel(colname{j},'FontName','Arial','FontSize',14);
        ylabel(colname{k},'FontName','Arial','FontSize',14);
        [r4,p4]=corr(data{4}(:,j),data{4}(:,k));
        r4=round(r4,2);%保留两位即可
        [r5,p5]=corr(data{5}(:,j),data{5}(:,k));
        r5=round(r5,2);%保留两位即可
        text(max(data{4}(:,j)),max(data{4}(:,k)),{strcat('easy RHO=',num2str(r4));plabelsymbol(p4);strcat('hard RHO=',num2str(r5));plabelsymbol(p5)},'FontSize',10);
        set(gca,'FontName','Arial','FontSize',14);
        i=i+1;
        
    end
end

end
function [paraB,pparaB,paraBdiff, pparaBdiff,delaybylength]=fGetpara(dir, transform_method, vardelay,treatmentgroup, ndifficultygroup)
% paraB=cell(5,ndelaygroup+1,length(group),3);
%第一维代表五种指标（cor rate, violation rate, miss rate, RT, licking consistency index LCI)，
%第二维代表总、短、长delay，根据ndelaygroup确定组数;若vardelay为自定义划分标准，则默认分成两组
%第三维代表两种实验条件（saline,CNO)，
%第四维是难度（end/probe)
% pparaB=ones(size(paraB,1),size(paraB,2),size(paraB,3)/2,size(paraB,4));%保存各个paraB中的参数CNO-saline的p值
% paraBdiff=cell(5,ndelaygroup+1,length(group)/2,3);%区分CNO-saline的效果
% pparaBdiff=ones(2,size(paraB,1));%保存paraB中的参数CNO-saline的值在区分delay/difficulty时的p值,只为作图服务，索引对应fig位置,第一行不分组，不需比较，故此处第一维代表分delay，第二维代表分difficulty
% delayPart 根据ndelaygroup确定的每个部分的delay长度
if length(vardelay)==1
    ndelaygroup=vardelay;
elseif length(vardelay)==4
    ndelaygroup=2;
    delaysetting=vardelay;
end
if strcmp(transform_method,'uni')
    paraB=cell(5,ndelaygroup+1,2*length(treatmentgroup),ndifficultygroup+1);%第三维区分ipsi, contra来看
    paraBdiff=cell(5,ndelaygroup+1,length(treatmentgroup),ndifficultygroup+1);%区分CNO-saline的效果
else
    paraB=cell(5,ndelaygroup+1,length(treatmentgroup),ndifficultygroup+1);%第一维代表五种指标，第二维代表总、短、长delay，第三维代表两种实验条件（saline,CNO)，第四维是难度（end/probe)
    paraBdiff=cell(5,ndelaygroup+1,length(treatmentgroup)/2,ndifficultygroup+1);
end
for n=1:length(treatmentgroup)
    [animal_name,dataChoiceResult]=fFindChoiceMat(strcat(dir,treatmentgroup{n}));
    nsample=size(dataChoiceResult,1)-1;
    for i=1:nsample
        %delay先确定好
        if length(vardelay)==1
            [delaybylength,delay]=fDelayGroup(dataChoiceResult{i}(:,4),'ngroup',ndelaygroup);%汇总成二维数组，每一列代表不同长度的delay,% logical 相加就变成double
        elseif length(vardelay)==4
            [delaybylength,delay]=fDelayGroup(dataChoiceResult{i}(:,4),'user-defined',delaysetting);
        end
        %获取trial results
        correct=dataChoiceResult{i}(:,2)==1;%取值1-4分别代表cor,err,miss.err
        error=dataChoiceResult{i}(:,2)==2;
        miss=dataChoiceResult{i}(:,2)==3;
        violation=dataChoiceResult{i}(:,2)==4;
        novio=correct+error;
        do=novio+violation;
        %获取endtrial/probetrial
        click=unique(dataChoiceResult{i}(:,1));
        sort(click);
        pend=logical((dataChoiceResult{i}(:,1)==click(1))+(dataChoiceResult{i}(:,1)==click(end)));
        pprobe=logical(ones(size(dataChoiceResult{i},1),1)-pend);
        difficulty=[pend+pprobe,pend,pprobe];%汇总成二维数组，每一列代表不同难度
        %获取reaction time
        indRT=double((dataChoiceResult{i}(:,5)>0).*(dataChoiceResult{i}(:,2)==1));%只看正确的trial
        indLCI=double(dataChoiceResult{i}(:,2)==1);%只看正确的trial
%         indLCI=double(dataChoiceResult{i}(:,2)==2);%只看错误的trial
%         indLCI=double(dataChoiceResult{i}(:,2)~=3);%看所有trial，除去不舔水的miss
        for k=1:size(difficulty,2)%区分难度
            %获取ipsi/contra
            if strcmp(transform_method,'uni')
                dataChoiceResult{i}(:,1)=dataChoiceResult{i}(:,3);%左右调换的话就用第六列的声音对应的左右替换第一列频率
                protate=logical((dataChoiceResult{i}(:,6)==2)+(dataChoiceResult{i}(:,6)==3));%这些trials需要旋转，也就是上下，左右各颠倒一次
                dataChoiceResult{i}(protate,1)=7-dataChoiceResult{i}(protate,1);
                ipsi=(dataChoiceResult{i}(:,1)<=3);
                contra=(dataChoiceResult{i}(:,1)>=3);
                for j=1:size(delay,2)%各种长短的delay
                    paraB{1,j,n,k}(i)=100*sum(correct.*delay(:,j).*ipsi.*difficulty(:,k))/sum(novio.*delay(:,j).*ipsi.*difficulty(:,k));
                    paraB{2,j,n,k}(i)=100*sum(violation.*delay(:,j).*ipsi.*difficulty(:,k))/sum(do.*delay(:,j).*ipsi.*difficulty(:,k));
                    paraB{3,j,n,k}(i)=100*sum(miss.*delay(:,j).*ipsi.*difficulty(:,k))/sum(delay(:,j).*ipsi.*difficulty(:,k));
                    paraB{4,j,n,k}(i)=sum(indRT.*dataChoiceResult{i}(:,5).*delay(:,j).*ipsi.*difficulty(:,k))/sum(indRT.*delay(:,j).*ipsi.*difficulty(:,k));
                    paraB{5,j,n,k}(i)=sum(indLCI.*dataChoiceResult{i}(:,7).*delay(:,j).*ipsi.*difficulty(:,k))/sum(indLCI.*delay(:,j).*ipsi.*difficulty(:,k));
                    paraB{1,j,n+2,k}(i)=100*sum(correct.*delay(:,j).*contra.*difficulty(:,k))/sum(novio.*delay(:,j).*contra.*difficulty(:,k));
                    paraB{2,j,n+2,k}(i)=100*sum(violation.*delay(:,j).*contra.*difficulty(:,k))/sum(do.*delay(:,j).*contra.*difficulty(:,k));
                    paraB{3,j,n+2,k}(i)=100*sum(miss.*delay(:,j).*contra.*difficulty(:,k))/sum(delay(:,j).*contra.*difficulty(:,k));
                    paraB{4,j,n+2,k}(i)=sum(indRT.*dataChoiceResult{i}(:,5).*delay(:,j).*contra.*difficulty(:,k))/sum(indRT.*delay(:,j).*contra.*difficulty(:,k));
                    paraB{5,j,n+2,k}(i)=sum(indLCI.*dataChoiceResult{i}(:,7).*delay(:,j).*contra.*difficulty(:,k))/sum(indLCI.*delay(:,j).*contra.*difficulty(:,k));
                end
            else
                for j=1:size(delay,2)%各种长短的delay
                    paraB{1,j,n,k}(i)=100*sum(correct.*delay(:,j).*difficulty(:,k))/sum(novio.*delay(:,j).*difficulty(:,k));
                    paraB{2,j,n,k}(i)=100*sum(violation.*delay(:,j).*difficulty(:,k))/sum(do.*delay(:,j).*difficulty(:,k));
                    paraB{3,j,n,k}(i)=100*sum(miss.*delay(:,j).*difficulty(:,k))/sum(delay(:,j).*difficulty(:,k));
                    paraB{4,j,n,k}(i)=sum(indRT.*dataChoiceResult{i}(:,5).*delay(:,j).*difficulty(:,k))/sum(indRT.*delay(:,j).*difficulty(:,k));
                    paraB{5,j,n,k}(i)=sum(indLCI.*dataChoiceResult{i}(:,7).*delay(:,j).*difficulty(:,k))/sum(indLCI.*delay(:,j).*difficulty(:,k));
                end
            end
        end
    end
end

pparaB=ones(size(paraB,1),size(paraB,2),size(paraB,3)/2,size(paraB,4));%保存各个paraB中的参数CNO-saline的p值
for i1=1:size(pparaB,1)%parameter
    for i2=1:size(pparaB,2)%delay
        for i4=1:size(pparaB,4)%difficulty
            paraBdiff{i1,i2,1,i4}=paraB{i1,i2,2,i4}-paraB{i1,i2,1,i4};%CNO-saline
            if size(paraB, 3)==2
                [~,p]=ttest(paraB{i1,i2,1,i4},paraB{i1,i2,2,i4});
                %                     p=round(p,4);%保留小数点后4位
                pparaB(i1,i2,1,i4)=p;
            else%size(paraB, 3)==4
                [~,pipsi]=ttest(paraB{i1,i2,1,i4},paraB{i1,i2,2,i4});
                pparaB(i1,i2,1,i4)=pipsi;
                [~,pcontra]=ttest(paraB{i1,i2,1,i4},paraB{i1,i2,2,i4});
                pparaB(i1,i2,2,i4)=pcontra;
                paraBdiff{i1,i2,3,i4}=paraB{i1,i2,4,i4}-paraB{i1,i2,3,i4};
            end
        end
    end
end
%     compare difference CNO-saline
pparaBdiff=ones(2,size(paraB,1));%保存paraB中的参数CNO-saline的值在区分delay/difficulty时的p值,只为作图服务，索引对应fig位置,第一行不分组，不需比较，故此处第一维代表分delay，第二维代表分difficulty
for i=1:size(paraB,1)
    [~,pbydelay]=ttest(paraBdiff{i,2,1,1},paraBdiff{i,end,1,1});
    pparaBdiff(1,i)=pbydelay;
    [~,pbydifficulty]=ttest(paraBdiff{i,1,1,2},paraBdiff{i,1,1,end});
    pparaBdiff(2,i)=pbydifficulty;
end
end
function [fig]=fPlot2FactorInteraction(data,plotMethod,rowname,colname)
%plotMethod={'scatter','line','3D'}
%data is 3-D cells. 1-D treatment group(saline/CNO), 2-D 1st factor, 3-D 2nd factor
%each row for 1st factors, each colum for 2nd factors.
if data{1,1,1}>100 %if true, then is RT data, otherwise is rate data
    unit='(ms)';
else
%     cellfun(@(x) x*100,data);
    unit='(%)';
end
nrow=size(data,2);
ncol=size(data,3);
if nrow~=length(rowname)
    warning('data and row name length not match');
end
if ncol~=length(colname)
    warning('data and column name length not match');
end
fig=figure;
if ~strcmp(plotMethod, '3D')
    set(gcf,'position',[0,0,300*nrow,300*ncol]);
    for i=1:nrow
        for j=1:ncol
            subplot(nrow,ncol,ncol*i-ncol+j);
            if strcmp(plotMethod,'scatter')
                scatter(data{1,i,j},data{2,i,j},20,'k');
                hold on;
                plot(0:100,0:100,'--k');
                hold on;
                ystr=strcat('CNO',unit);
                xlabel(strcat('saline',unit));
            elseif strcmp(plotMethod,'line')
                x=ones(length(data{1,i,j}),1);
                for k=1:size(data,1)
                    scatter(k*x,data{k,i,j},20,'k');
                    hold on;
                end
                for k=1:length(data{1,i,j})
                    plot([1,2],[data{1,i,j}(k), data{2,i,j}(k)],'k','LineWidth',0.5);
                    hold on;
                end
                ystr=unit;
                set(gca,'XTick',[1,2]);
                set(gca,'XTickLabel',{'saline','CNO'});
                set(gca,'xlim',[0,3]);
                set(gca,'ylim',[40,100]);
            end
            if i==1
                title(colname{j});
            end
            if j==1
                ylabel({rowname{i};ystr});
            else
                ylabel(ystr);
            end
            set(gca,'FontName','Arial','FontSize',14);
        end
    end
else
    for i=1:nrow
        for j=1:ncol
            ndot=length(data{1,i,j});
            scatter3(i*ones(ndot,1),j*ones(ndot,1),data{2,i,j}-data{1,i,j},20,'k');
            hold on;
        end
    end
    set(gca,'xlim',[0,nrow+1]);
    set(gca,'XTick',1:1:nrow);
    set(gca,'XTickLabel',rowname);
    set(gca,'ylim',[0,ncol+1]);
    set(gca,'YTick',1:1:ncol);
    set(gca,'YTickLabel',colname);
    zlabel(unit);
    zeroplane=patch([0,nrow+1,nrow+1,0],[0,0,ncol+1,ncol+1],[0,0,0,0]);
    alpha(zeroplane,0.3);
%     for i=1:2:360%旋转3D视图
%         %view(a,b):a是角度，b是仰视角
%         view(i,20);
%         pause(0.06);
%     end
end

end
function [delaybylength,delaygroup]= fDelayGroup(delayraw,varargin)
if nargin==3
    if strcmp(varargin{1},'ngroup')
        ngroup=varargin{2};
        %divide the delay into ngroup groups according to delay length
        mindelay=min(delayraw);
        maxdelay=max(delayraw);
        part=(maxdelay-mindelay)/ngroup;
        delaybylength=mindelay:part:maxdelay;
        delaygroup=ones(length(delayraw),1);
        for i=1:ngroup
            delaycol=(mindelay+(i-1)*part<=delayraw).*(mindelay+i*part>delayraw);
            delaygroup=[delaygroup delaycol];
        end
    elseif strcmp(varargin{1},'user-defined')
        delaybylength=varargin{2};
        if length(delaybylength)~=4
            warning('Input should be a 1*4 double vector');
        else
            delaygroup=ones(length(delayraw),1);
            for i=1:2
                delaycol=(delaybylength(2*i-1)<=delayraw).*(delaybylength(2*i)>=delayraw);
                delaygroup=[delaygroup delaycol];
            end
        end
    else
        warning('Input error for function fDelayGroup');
    end
else
    warning('Input number error for function fDelayGroup');
end
end
function [str]=plabel(pvalue)
if pvalue<0.05 && pvalue>=0.01
    str='p<0.05';
elseif pvalue<0.01 && pvalue>=0.001
    str='p<0.01';
elseif pvalue<0.001
    str='p<0.001';
else
    pvalue=round(pvalue,2);%保留两位即可
    str=strcat('p=',num2str(pvalue));
end
end
function [str]=plabelsymbol(pvalue)
if pvalue<0.05 && pvalue>=0.01
    str=' *';
elseif pvalue<0.01 && pvalue>=0.001
    str=' **';
elseif pvalue<0.001
    str=' ***';
else
    pvalue=round(pvalue,2);%保留两位即可
    str=strcat('p=',num2str(pvalue));
end
end