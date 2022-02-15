%直接读beh文件，进行单个seesion 或者 paired session绘图
clear;
set_n_session=1;%how many sessions to be plotted
path='E:\xulab\behavior\pyx385';
% path='E:\2P\CD088_20190111\im_data_reg\result_save';
temp=strsplit(path,'\');
figname=temp{end};
[animal_name_cell,dataChoice]=fFindChoiceMat(path,'session');%cell array of choice matrix, each cell represent a session
n_session=min(length(dataChoice)-1,set_n_session);%choose the min of set session number and actual data session number

%获取实验组和对照组的两组数据，并画图 
 choice_last=dataChoice{end-1};  
 animal_name=animal_name_cell{end-1};
[correctProbe_last,toneOct_last,choicep_last,f_last,x_last,yfit_last]=fProbePerformance(choice_last);
A=figure; 

if length(dataChoice)>2   
    for i=length(dataChoice)-n_session:length(dataChoice)-2
        color_session=(length(dataChoice)-1-i)/n_session*[1,1,1];
        choice_previous=dataChoice{i};
        animal_name=animal_name_cell{i};
        [correctProbe_previous,toneOct_previous,choicep_previous,f_previous,x_previous,yfit_previous]=fProbePerformance(choice_previous);
        curve(i-length(dataChoice)+n_session+2)=fPlotProbePerformance(A,color_session,correctProbe_previous,toneOct_previous,choicep_previous,x_previous,yfit_previous,animal_name);
        lgdstr{i-length(dataChoice)+n_session+2}=['day',num2str(length(dataChoice)-n_session-1-i)];
    end
    curve(1)=fPlotProbePerformance(A,[0,0,0],correctProbe_last,toneOct_last,choicep_last,x_last,yfit_last,animal_name); 
    lgdstr{1}='day0';
    figure(A);
    subplot(1,4,4);
    hlegend=legend(curve,lgdstr,'Location','Best');
    set(hlegend,'Box','off');
else
    curve=fPlotProbePerformance(A,[0,0,0],correctProbe_last,toneOct_last,choicep_last,x_last,yfit_last,animal_name); 
    lgdstr='day0';
    figure(A);
    subplot(1,4,4);
    hlegend=legend(curve,lgdstr,'Location','Best');
    set(hlegend,'Box','off');
end
saveas(A,strcat(path,'\behavior of ',figname),'png');
function [curve] = fPlotProbePerformance(fig,color,correctProbe,toneOct,choicep,x,yfit,animal_name)
%x1 = linspace(min(toneOct)-0.1, max(toneOct)+0.1, 100);
    figure(fig);
    set(gcf, 'position', [100 100 1200 300]);%控制fig尺寸
    subplot(1,4,1);
    scatter(2.^toneOct*20,correctProbe(1:size(correctProbe,1)-3,2)*100,50,color);
    hold on;
    scatter([2.^(toneOct(end)+1)*20,2.^(toneOct(end)+2)*20,2.^(toneOct(end)+3)*20],correctProbe(size(correctProbe,1)-2:size(correctProbe,1),2)*100,50,color,'filled');
    scatter(1,correctProbe(end,2)*100,50,color,'filled');
    hold on;    
    ylabel('% Correct','FontName','Arial','FontSize',14);
    set(gca,'XTick',1,'LineWidth',1);
    set(gca,'XTickLabel',' ');
    %xlabel('Log_2 (click rate/20)','FontName','Arial','FontSize',14);
    xlabel('clicks/s','FontName','Arial','FontSize',14);
    title(['Performance of ', animal_name],'FontName','Arial','FontWeight','Bold','FontSize',16);
    plot([2.^(toneOct(1)-0.5)*20, 2.^(toneOct(end)+0.5)*20],[80,80],'Color',color,'LineWidth',1);
    hold on;
    set(gca,'XScale','log');
    set(gca,'XTick',[20,125],'LineWidth',1);
    set(gca,'XTickLabel',[20,125]);
    set(gca, 'XLim',[2.^(toneOct(1)-0.5)*20, 2.^(toneOct(end)+3.5)*20]);
    set(gca, 'YLim',[0 100]);
    set(gca,'FontName','Arial','FontSize',14);
    
    subplot(1,4,2);
    plot(2.^toneOct*20 ,correctProbe(1:size(correctProbe,1)-3,3)*100,'Color',color,'LineWidth',2);
    hold on;
    %xlabel('Training day','FontName','Arial','FontSize',14);
    ylabel('% Miss','FontName','Arial','FontSize',14);
    %xlabel('Log_2 (click rate/20)','FontName','Arial','FontSize',14);
    xlabel('clicks/s','FontName','Arial','FontSize',14);
    title(['Miss rate of ' ,animal_name],'FontName','Arial','FontWeight','Bold','FontSize',16);
    hold on;
    set(gca,'XScale','log');
    set(gca,'XTick',[20,125],'LineWidth',1);
    set(gca,'XTickLabel',[20,125]);
    set(gca, 'XLim',[2.^(toneOct(1)-0.5)*20, 2.^(toneOct(end)+0.5)*20]);
    
    subplot(1,4,3);
    plot(2.^toneOct*20 ,correctProbe(1:size(correctProbe,1)-3,5)*100,'Color',color,'LineWidth',2);
    hold on;
    %xlabel('Training day','FontName','Arial','FontSize',14);
    ylabel('% Violation','FontName','Arial','FontSize',14);
    xlabel('clicks/s','FontName','Arial','FontSize',14);
%     title(['Violation rate of ' ,animal_name],'FontName','Arial','FontWeight','Bold','FontSize',16);
    set(gca,'XScale','log');
    set(gca,'XTick',[20,125],'LineWidth',1);
    set(gca,'XTickLabel',[20,125]);
    set(gca, 'XLim',[2.^(toneOct(1)-0.5)*20, 2.^(toneOct(end)+0.5)*20]);
     
    subplot(1,4,4);
    box off;
    points_legend=scatter(2.^choicep(:,1)*20,choicep(:,2)*100,correctProbe(1:size(correctProbe,1)-3,6),color,'filled');
    hold on;
    pointsize=sort(correctProbe(1:size(correctProbe,1)-3,6));
    %pointsizelevel=ceil(pointsize(end-2)/10)*10;%find probe trial numbers,choose the largest and decide legend markersize.
    pointsizelevel=ceil(pointsize(1)/10)*10;%find probe trial numbers,choose the smallest and decide legend markersize.
    [lgd,icons,~,~]=legend(points_legend,strcat(num2str(pointsizelevel),' trials'),'Location','Best');
    set(lgd,'Box','off');
    icons(2).Children.MarkerSize = sqrt(pointsizelevel);%legend as example how big the marker is
    %errorbar(choicep(:,1),choicep(:,2),correctProbe(:,4),'LineWidth',0);
    %hold on;
    curve=plot(x ,yfit*100 ,'Color',color,'LineWidth',2);%标记该曲线，为legend提供名称
    hold on;
%     ylabel('% Choice of High Click Rate Side','FontName','Arial','FontSize',14);
ylabel('% Choice Right','FontName','Arial','FontSize',14);
    xlabel('clicks/s','FontName','Arial','FontSize',14);
%     title(['Psychometric curve'],'FontName','Arial','FontWeight','Bold','FontSize',16);
    set(gca,'FontName','Arial','FontSize',14);
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
    click(isnan(click))=[];
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