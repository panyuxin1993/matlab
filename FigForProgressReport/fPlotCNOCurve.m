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
curve=cell(length(group),1);%���ڴ洢ͼ����  
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
    %%%%%%%%%%%%%%%%%%�˴�datachoice�����ݽ���ͬsession�����Ժ�ᶪʧrule��Ϣ������ת��������޷���ת�������������datachoice
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
    set(gcf, 'position', [0 0 900 300]);%����fig�ߴ�
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
            curve=plot(x ,yfit*100 ,'Color',color,'LineWidth',2);%��Ǹ����ߣ�Ϊlegend�ṩ����
            hold on;
        else
            curve=plot(x ,yfit*100 ,'Color',color,'LineWidth',0.5);%��Ǹ����ߣ�Ϊlegend�ṩ����
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
     %delay��ȷ����
    maxdelay=max(choice(:,4));
    mindelay=min(choice(:,4));
    part=(maxdelay-mindelay)/3;%��delay�ֳ��ļ��ֳ���
    shortdelay=(choice(:,4)<=(mindelay+part));
    longdelay=(choice(:,4)>(maxdelay-part));
    delay=[ones(size(choice,1),1),shortdelay,longdelay];%���ܳɶ�ά���飬ÿһ�д���ͬ���ȵ�delay,% logical ��Ӿͱ��double
    correctProbe=zeros(size(delay,2),9,7);%��һά����delay���ȣ�֮��ڶ�ά��ʾһ��Ƶ�ʣ�����ά�ֱ��ʾƵ��(��������ҵ����Ļ�����������Ӧ������)����ȷ�ʣ�miss�ʣ���ȷ�ʵ�std��violation;RT,LCI�ʣ��ڶ�ά�����������У��ֱ��ʾend trial,��probe trial, all trials
    click=[20 28 40 62 90 125];
    if strcmp(transform_method,'left-right')
        choice(:,1)=choice(:,3);%���ҵ����Ļ����õ����е�������Ӧ�������滻��һ��Ƶ��
    elseif strcmp(transform_method,'top-down')
        
    elseif strcmp(transform_method,'uni')
        choice(:,1)=choice(:,3);%���ҵ����Ļ����õ����е�������Ӧ�������滻��һ��Ƶ��
        protate=logical((choice(:,6)==2)+(choice(:,6)==3));%��Щtrials��Ҫ��ת��Ҳ�������£����Ҹ��ߵ�һ��
        choice(protate,1)=7-choice(protate,1);      
    end
    pclick=zeros(size(choice,1),9);%���ڼ�¼ÿ��click rate��trialλ��
    for i=1:6
        pclick(:,i)=double(choice(:,1)==click(i));
    end
    pclick(:,7)=pclick(:,1)+pclick(:,6);
    pclick(:,8)=pclick(:,2)+pclick(:,3)+pclick(:,4)+pclick(:,5);
    pclick(:,9)=pclick(:,7)+pclick(:,8);
    pperformance=zeros(size(choice,1),4);%���ڼ�¼�������͵�trialλ�ã�����correct,error,miss��violation 4��
    performance=[1 2 3 4];%�ֱ��ʾCORRECT,ERROR,MISS,VIOLATION ��right(ע��action choice ����Ϊ1)
    for i=1:4
        pperformance(:,i)=double(choice(:,2)==performance(i));
    end
    choicep=cell(size(delay,2),1);
    fitdata=cell(size(delay,2),1);
    f=cell(size(delay,2),1);
    for j=1:size(delay,2)%���ֲ�ͬdelay
        correctProbe(j,1:6,1)=click;%�������������end ��probez����ʱ����ֵ��
        for i=1:9
            do=pperformance(:,1)+pperformance(:,2); 
            notMiss=do+pperformance(:,4);
            correctProbe(j,i,2)=sum(pperformance(:,1).*pclick(:,i).*delay(:,j))/sum(do.*pclick(:,i).*delay(:,j));
            correctProbe(j,i,3)=sum(pperformance(:,3).*pclick(:,i).*delay(:,j))/sum(pclick(:,i).*delay(:,j));
            %�������Ƶ�ʵ���ȷ�ʵ�std
            correctNoMiss=pperformance(:,1).*pclick(:,i).*delay(:,j);
            errorNoMiss=pperformance(:,2).*pclick(:,i).*delay(:,j);
            pcor=sum(correctNoMiss)/sum(correctNoMiss+ errorNoMiss);
            correctProbe(j,i,4)=sqrt(sum(correctNoMiss+ errorNoMiss)*pcor*(1-pcor)/(sum(correctNoMiss+ errorNoMiss)-1));
%             correctNoMiss(~notMiss)=nan;
%             correctProbe(j,i,4)=nanstd(correctNoMiss);
            correctProbe(j,i,5)=sum(pperformance(:,4).*pclick(:,i).*delay(:,j))/sum(notMiss.*pclick(:,i).*delay(:,j));%%%�˴�����population����������󣬻���ƫС     
            correctProbe(j,i,6)=sum(choice(:,5).*pperformance(:,1).*pclick(:,i).*delay(:,j))/sum(pperformance(:,1).*pclick(:,i).*delay(:,j));%RT only include correct trials
            correctProbe(j,i,7)=sum(choice(:,5).*pperformance(:,1).*pclick(:,i).*delay(:,j))/sum(pperformance(:,1).*pclick(:,i).*delay(:,j));%LCI also only include correct trials
%             correctProbe(j,i,7)=sum(choice(:,5).*pclick(:,i).*delay(:,j))/sum(pclick(:,i).*delay(:,j));%LCI include all trials(cor,err,miss,vio)
        end

        %����choicep����ʾѡ��ĳ������ı���������curve����Ϊɢ��
        %����toneOct����ʾ��ͼ���õĺ�����
        %Ϊ�˷�����������±�ȵ�һ���ԣ���Ԥ��Ϊ�̼���ͬ�߶�
        if strcmp(transform_method,'left-right')
            correctProbe(j,1:6,1)=[20;28;40;62;90;125];
        end
        toneOct  = log2(correctProbe(j,1:6,1)/correctProbe(j,1,1));
        choicep{j} =squeeze(correctProbe(j,1:6,1:2));%6*2����
        choicep{j}(:,1)=toneOct;
        choicep{j}(1:3,2)=1-correctProbe(j,1:3,2);
        
        %���
        ffun=fittype('a+b./(1+exp(-k.*(x-c)))','independent','x');
        %����һ������average�����ܵ���ȷ�ʼ��㣬��ÿ��stimulusһ����
%         f{j} =fit(choicep{j}(:,1),choicep{j}(:,2),ffun,'Startpoint',[0,1,1,1]);
        %�����������õ���trial����ϣ�һ��trialһ���㣬�Ա���Ϊstimulus�������0-1��
        %���ַ������end point�������� 
        fitdata{j}=choice(logical((choice(:,2)<3).*delay(:,j)),1:2);%��һ��1-6�������ִ̼����ڶ��зֱ��ʾCORRECT,ERROR,MISS,VIOLATION��ֻ����correct,error����1,2
        if strcmp(transform_method,'left-right')
            fitdata{j}((fitdata{j}(:,1)<50),2)=fitdata{j}((fitdata{j}(:,1)<50),2)-1;%ǰ����stimulus��Ӧ��ѡ��ת��Ϊ0���ԣ���-1�����ң�
            fitdata{j}((fitdata{j}(:,1)>50),2)=2-fitdata{j}((fitdata{j}(:,1)>50),2);%������stimulus��Ӧ��ѡ�ң�ת��Ϊ0������-1���ԣ��ң�    
        elseif strcmp(transform_method,'top-down')
            fitdata{j}((fitdata{j}(:,1)<50),2)=fitdata{j}((fitdata{j}(:,1)<50),2)-1;%ǰ����stimulus��Ӧ��ѡ��Ƶ�࣬ת��Ϊ0���ԣ��ͣ�-1�����ߣ�
            fitdata{j}((fitdata{j}(:,1)>50),2)=2-fitdata{j}((fitdata{j}(:,1)>50),2);%������stimulus��Ӧ��ѡ��Ƶ�࣬ת��Ϊ0�����ͣ�-1���ԣ��ߣ�
        end
        fitdata{j}(:,1)=log2(fitdata{j}(:,1)/min(fitdata{j}(:,1)));%stimulus��ת��Ϊ����
        f{j} =fit(fitdata{j}(:,1),fitdata{j}(:,2),ffun,'Startpoint',[0,1,1,1]);
    end
end

