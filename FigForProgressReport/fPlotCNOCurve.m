function [ fig ] = fPlotCNOCurve( fig, position, exp )
%FPLOTCNOCURVE plot curves of different experiment, distinguish delay
%length
%   Input: exp decide which exp(ip/infusion)
if strcmp(exp,'ip')
    dirC='D:\xulab\behavior\CNO experiment\CNO ctrl\ip\';
    dir='D:\xulab\behavior\CNO experiment\ip\';
elseif strcmp(exp,'infusion')
    dir='D:\xulab\behavior\CNO experiment\infusion\bilateral\';
    dirC='D:\xulab\behavior\CNO experiment\CNO ctrl\infusion\bilateral\';
end
group={'saline', 'CNO'};
color={'k','r'}; 
type=  {'single case','population'};
curve=cell(length(group),1);%用于存储图像句柄  
% transform_method={'top-down'};
transform_method={'left-right'};
for n=1:length(group)
    [animal_name,dataChoiceResult]=fFindChoiceMat(strcat(dir,group{n}));
    nsample=size(dataChoiceResult,1)-1;
    
    for i=1:nsample
        [correctProbe,toneOct,choicep,f]=fProbePerformance(dataChoiceResult{i},transform_method);
        fPlotProbePerformance(fig,color{n},correctProbe,toneOct,choicep,f,animal_name{end},type{1},transform_method);
    end
    %population mean
    %%%%%%%%%%%%%%%%%%此处datachoice的数据将不同session汇总以后会丢失rule信息，左右转的情况下无法旋转；考虑重新设计datachoice
    [correctProbe,toneOct,choicep,f]=fProbePerformance(dataChoiceResult{nsample+1},transform_method);
    curve{n}=fPlotProbePerformance(fig,color{n},correctProbe,toneOct,choicep,f,animal_name{end},type{2},transform_method);
end
leg=legend([curve{1} curve{2}],'saline','CNO');
set(leg,'position',[0.95,0.5,0.1,0.1]);
legend boxoff;
set(fig,'PaperPosition',position);

end

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
    set(gcf, 'position', [0 0 900 300]);%控制fig尺寸
    titlestr={'All','Short delay','Long delay'};
    for i=1:3
        subplot(1,3,i);
        title(titlestr{i},'FontName','Arial','FontSize',14);
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
            if i==1
            ylabel('% Choice of Right','FontName','Arial','FontSize',14);
            end
            xlabel('stimulus side','FontName','Arial','FontSize',14);
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

