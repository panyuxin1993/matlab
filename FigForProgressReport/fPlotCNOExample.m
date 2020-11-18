function [ fig ] = fPlotCNOExample( fig,position,exp )
%FPLOTCNOEXAMPLE plot example of CNO inhibition
%   
if strcmp(exp,'ip')
    dir='D:\xulab\behavior\exp_data\pyx011\probe';
elseif strcmp(exp,'infusion')
    dir='D:\xulab\behavior\exp_data\pyx019\probe';
end
[animal_name,dataChoice]=fFindChoice(dir,2);
b=dataChoice;
%correctProbe=fCorrectRateProbePerSession(choice);
%��ȡʵ����Ͷ�������������ݣ�����ͼ 
 p_ctrl=find(dataChoice(1,:,1)==0);%������bug��������1000trialʱ�ǿ�ֵ
 if isempty(p_ctrl) 
     p_ctrl=1000;
 end
 choice_ctrl=squeeze(dataChoice(1,1:p_ctrl-1,:));  
[correctProbe_ctrl,toneOct_ctrl,choicep_ctrl,f_ctrl,x_ctrl,yfit_ctrl]=fProbePerformance(choice_ctrl);

ctrl_curve=fPlotProbePerformance(fig,'k',correctProbe_ctrl,toneOct_ctrl,choicep_ctrl,x_ctrl,yfit_ctrl,animal_name); 
hold on;
if size(dataChoice,1)>1    
    p_exp=find(dataChoice(2,:,1)==0);    
    choice_exp=squeeze(dataChoice(2,1:p_exp-1,:));  
    [correctProbe_exp,toneOct_exp,choicep_exp,f_exp,x_exp,yfit_exp]=fProbePerformance(choice_exp);
    exp_curve=fPlotProbePerformance(fig,'r',correctProbe_exp,toneOct_exp,choicep_exp,x_exp,yfit_exp,animal_name);
    hlegend=legend([ctrl_curve, exp_curve],'ctrl','ip CNO','Location','Best');
    set(hlegend,'Box','off');
    set(hlegend,'position',[0.95,0.4,0.1,0.1]);
end
set(fig,'PaperPosition',position);
end
function [curve] = fPlotProbePerformance(fig,color,correctProbe,toneOct,choicep,x,yfit,animal_name)
%x1 = linspace(min(toneOct)-0.1, max(toneOct)+0.1, 100);
    figure(fig);
    set(gcf, 'position', [100 100 1200 300]);%����fig�ߴ�

    box off;
    points_legend=scatter(2.^choicep(:,1)*20,choicep(:,2)*100,50,color,'filled');
    hold on;
%     pointsize=sort(correctProbe(1:size(correctProbe,1)-3,6));
    %pointsizelevel=ceil(pointsize(end-2)/10)*10;%find probe trial numbers,choose the largest and decide legend markersize.
%     pointsizelevel=ceil(pointsize(1)/10)*10;%find probe trial numbers,choose the smallest and decide legend markersize.
%     [lgd,icons,~,~]=legend(points_legend,strcat(num2str(pointsizelevel),' trials'),'Location','Best');
%     set(lgd,'Box','off');
%     icons(2).Children.MarkerSize = sqrt(pointsizelevel);%legend as example how big the marker is
    %errorbar(choicep(:,1),choicep(:,2),correctProbe(:,4),'LineWidth',0);
    %hold on;
    curve=plot(x ,yfit*100 ,color,'LineWidth',2);%��Ǹ����ߣ�Ϊlegend�ṩ����
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
    if length(click)==2
        warning('Not enough stimuli');
    end
    freq_n=length(click);
    correctProbe=zeros(freq_n+3,6);%ÿ�б�ʾһ��Ƶ�ʣ�ÿ�зֱ��ʾƵ�ʣ���ȷ�ʣ�miss�ʣ���ȷ�ʵ�std��violation�ʣ�correct+error��trial����֮�ͣ������������У��ֱ��ʾend trial,��probe trial, all trials

    pclick=zeros(size(choice,1),freq_n);%���ڼ�¼ÿ��click rate��trialλ��
    for i=1:freq_n
        pclick(:,i)=double(choice(:,1)==click(i));
    end
    pclick(:,freq_n+1)=pclick(:,1)+pclick(:,freq_n);
    pclick(:,freq_n+2)=sum(pclick(:,2:freq_n-1),2);%�������
    pclick(:,freq_n+3)=pclick(:,freq_n+1)+pclick(:,freq_n+2);
    pperformance=zeros(size(choice,1),4);%���ڼ�¼�������͵�trialλ�ã�����correct,error,miss����
    performance=[1 2 3 4];%�ֱ��ʾCORRECT,ERROR,MISS,VIOLATION
    for i=1:4
        pperformance(:,i)=double(choice(:,2)==performance(i));
    end
    correctProbe(1:freq_n,1)=click;%�������������end ��probez����ʱ����ֵ��
    for i=1:freq_n+3
        do=(pperformance(:,1)+pperformance(:,2)).*pclick(:,i);
        notMiss=(do+pperformance(:,4)).*pclick(:,i);
        correctProbe(i,2)=sum(pperformance(:,1).*pclick(:,i))/sum(do);
        correctProbe(i,3)=sum(pperformance(:,3).*pclick(:,i))/sum(pclick(:,i));
        %�������Ƶ�ʵ���ȷ�ʵ�std
        correctNoMiss=pperformance(:,1).*pclick(:,i);
        correctNoMiss(~notMiss)=nan;
        correctProbe(i,4)=nanstd(correctNoMiss);
        correctProbe(i,5)=sum(pperformance(:,4).*pclick(:,i))/sum(notMiss);
        correctProbe(i,6)=sum(do);
    end
    
    %���
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
