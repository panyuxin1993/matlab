function [ fig ] = fPlotCNOSummary( fig,position,exp )
%FPLOTCNOSUMMARY plot CNO ip/infusion effect, compare between delay length,
%difficulties, both for hM4D mice and mcherry mice

if strcmp(exp,'ip')
    dirC='D:\xulab\behavior\CNO experiment\CNO ctrl\ip\';
    dirE='D:\xulab\behavior\CNO experiment\ip\';
elseif strcmp(exp,'infusion')
    dirE='D:\xulab\behavior\CNO experiment\infusion\bilateral\';
    dirC='D:\xulab\behavior\CNO experiment\CNO ctrl\infusion\bilateral\';
end
group={'saline', 'CNO'};
color={'k','r'}; 
type=  {'single case','population'};
curve=cell(length(group),1);%���ڴ洢ͼ����  
transform_method={};
ndelaygroup=2;
delaysetting=[300 900 900 1500];
ndifficultygroup=2;
[paraE,pparaE,paraEdiff,pparaEdiff,delaybylength]=fGetpara(dirE, transform_method, ndelaygroup ,group, ndifficultygroup);
[paraC,pparaC,paraCdiff,pparaCdiff,delaybylength]=fGetpara(dirC, transform_method, ndelaygroup ,group, ndifficultygroup);
fig=fCmpOverallPerformance(paraE,paraEdiff,paraC,paraCdiff,ndelaygroup);%��CNO-saline�Ĳ�ֵ���Ƚϲ�ͬ����µ�����
set(fig,'PaperPosition',position);
end

function [fig] = fCmpOverallPerformance(paraE,paraEdiff,paraC,paraCdiff,ndelaygroup)
    fig=figure;%2*4��ͼ���ֱ��correct rate, violation rate, miss rate, reaction time����һ���ܽ����  �ڶ�������delay����ֻ�ǳ��̣����Ƿֳ�200msһ���6��,����������difficulty
    set(gcf, 'position', [0 0 1400 800]);%����fig�ߴ�
    ylabelstr={'\it\Delta\rmcorrect rate(%)';'\it\Delta\rmViolation rate(%)';'\it\Delta\rmMiss rate(%)';'\it\Delta\rmReaction time(ms)'};
    titlestr={'Correct';'Violation';'Miss';'Reaction time'};
%     a=suptitle('CNO-saline');
%     set(a,'FontSize',14,'position',[0.5 0.95 0.02 0.01]);
    pararaw=zeros(length(paraE{1,1,1,1}),4);%��һά��ʾ��������ֵ��ʾCNO-saline������ڶ�ά������ָ��
    pararawC=zeros(length(paraC{1,1,1,1}),4);
    x=[0.8,1.2,1.8,2.2];%�����x����
    for i=1:4 %��������
        %���в��ϲ�
        subplot(3,4,i);%���кϲ�
        pararaw(:,i)=paraE{i,1,2,1}-paraE{i,1,1,1};
        pararawC(:,i)=paraC{i,1,2,1}-paraC{i,1,1,1};      
        scatter(ones(size(pararaw,1),1),pararaw(:,i),20,'r','filled');
        hold on;
        scatter(ones(size(pararawC,1),1)*2,pararawC(:,i),20,'k','filled');
        ytext=get(gca,'Ylim');    
        [~,p]=ttest(pararaw(:,i));
        text(1,ytext(end),plabelsymbol(p));
        [~,p]=ttest(pararawC(:,i));
        text(2,ytext(end),plabelsymbol(p));
        [~,p2]=ttest2(pararaw(:,i),pararawC(:,i));
        text(1.5,ytext(end)+0.2*(ytext(end)-ytext(1)),plabelsymbol(p2));
        plot([1,2],[ytext(end)+0.1*(ytext(end)-ytext(1)),ytext(end)+0.1*(ytext(end)-ytext(1))],'k-');
        set(gca,'YLim',[ytext(1), ytext(end)+0.1*(ytext(end)-ytext(1))]);%this way avoid ylim changing when changing fig size
        ylabel(ylabelstr{i});
        plot([0,3],[0,0],'--k','LineWidth',0.5);%0��ˮƽ��
        %     for i=1:size(pararaw,1)%����������ͬsession�Ĳ���
        %         plot(pararaw(i,:));
        %         hold on;
        %     end
        hold on;
        %     set(gca,'Ylim',[-400,600]);
        set(gca, 'XLim',[0 3]);
        set(gca,'XTick',[1,2]);
        set(gca,'XTickLabel',{'hM4D ','mCherry '});
        %�ĳ���ת���label�Խ�ʡ�ռ�
        xtb = get(gca,'XTickLabel');% ��ȡ���������ǩ���
        xt = get(gca,'XTick');% ��ȡ��������̶Ⱦ��
        yt = get(gca,'YTick'); % ��ȡ��������̶Ⱦ��
        xtextp=xt;%ÿ����ǩ����λ�õĺ����꣬�����ȻӦ�ú�ԭ����һ���ˡ�
        ytextp=yt(1)*ones(1,length(xt)); % ������ʾ��ǩ��λ�ã�д����Ψһ��������ʵ����Ϊÿ����ǩ�ҷ���λ�õ�������
        % rotation��������ת�Ƕȴ�����ʱ����ת����ת�������HorizontalAlignment�������趨��
        % ��3������ֵ��left��right��center��������Ը�������ֵ���Լ�rotation��ĽǶȣ�����д����45
        text(xtextp,ytextp,xtb,'HorizontalAlignment','right','rotation',30,'fontsize',8);
        set(gca,'xticklabel',[]);% ��ԭ�еı�ǩ��ȥ
        box off;
        
        subplot(3,4,4+i);%�ڶ���,delay
        parabydelay=zeros(length(paraEdiff{i,1,1,1}),ndelaygroup);
        parabydelayC=zeros(length(paraCdiff{i,1,1,1}),ndelaygroup);
        for j=2:size(paraEdiff,2)
            parabydelay(:,j-1)=paraEdiff{i,j,1,end}';
            parabydelayC(:,j-1)=paraCdiff{i,j,1,end}';
        end
        for j=1:size(parabydelay,1)
            plot([x(1),x(2)],parabydelay(j,:),'r','LineWidth',0.5);
            hold on;
        end
        for j=1:size(parabydelayC,1)
            plot([x(3),x(4)],parabydelayC(j,:),'k','LineWidth',0.5);
        end
        fige1=scatter(x(1)*ones(size(parabydelay,1),1),parabydelay(:,1),20,'r','filled');
        fige2=scatter(x(2)*ones(size(parabydelay,1),1),parabydelay(:,2),20,'r');
        figc1=scatter(x(3)*ones(size(parabydelayC,1),1),parabydelayC(:,1),20,'k','filled');
        figc2=scatter(x(4)*ones(size(parabydelayC,1),1),parabydelayC(:,2),20,'k');
        plot([0,3],[0,0],'--k','LineWidth',0.5);%0��ˮƽ��
        hold on;
        set(gca, 'XLim',[0 size(parabydelay,2)+1]);
        ylim=get(gca,'YLim');
        ylabel(ylabelstr{i});
        [~,p]=ttest(parabydelay(:,1),parabydelay(:,end));
        xtext=get(gca,'Xlim');
        ytext=get(gca,'Ylim');
        text(1,ytext(end),plabelsymbol(p));
        [~,p]=ttest(parabydelayC(:,1),parabydelayC(:,end));
        text(2,ytext(end),plabelsymbol(p));
        parabydelaydiff=parabydelay(:,end)-parabydelay(:,1);
        parabydelaydiffC=parabydelayC(:,end)-parabydelayC(:,1);
        [~,p2]=ttest2(parabydelaydiff,parabydelaydiffC);
        %        [~,p2]=ranksum(parabydelaydiff,parabydelaydiffC);
        text(1.5,ytext(end)+0.1*(ytext(end)-ytext(1)),plabelsymbol(p2));
        plot([1,2],[ytext(end)+0.01*(ytext(end)-ytext(1)),ytext(end)+0.01*(ytext(end)-ytext(1))],'k-');
        set(gca,'YLim',[ylim(1), ytext(end)+0.1*(ytext(end)-ytext(1))]);%this way avoid ylim changing when changing fig size
        set(gca,'XTick',[1,2]);
        set(gca,'XTickLabel',{'hM4D ','mCherry '});
        %�ĳ���ת���label�Խ�ʡ�ռ�
        xtb = get(gca,'XTickLabel');% ��ȡ���������ǩ���
        xt = get(gca,'XTick');% ��ȡ��������̶Ⱦ��
        yt = get(gca,'YTick'); % ��ȡ��������̶Ⱦ��
        xtextp=xt;%ÿ����ǩ����λ�õĺ����꣬�����ȻӦ�ú�ԭ����һ���ˡ�
        ytextp=yt(1)*ones(1,length(xt)); % ������ʾ��ǩ��λ�ã�д����Ψһ��������ʵ����Ϊÿ����ǩ�ҷ���λ�õ�������
        % rotation��������ת�Ƕȴ�����ʱ����ת����ת�������HorizontalAlignment�������趨��
        % ��3������ֵ��left��right��center��������Ը�������ֵ���Լ�rotation��ĽǶȣ�����д����45
        text(xtextp,ytextp,xtb,'HorizontalAlignment','right','rotation',30,'fontsize',8);
        set(gca,'xticklabel',[]);% ��ԭ�еı�ǩ��ȥ
        if i==4
            legend([fige1,fige2],'short delay','long delay');
            leg=legend([figc1,figc2],'short delay','long delay');
            set(leg,'position',[0.95 0.5 0.02 0.05]);
            legend boxoff;
        end
        %        set(gca,'FontName','Arial','FontSize',14);
        box off;
       
       subplot(3,4,8+i);%������,difficulty
       easy=paraEdiff{i,1,1,2};
       hard=paraEdiff{i,1,1,3};
       parabydifficulty=[easy' hard'];
       parabydifficultyC=[paraCdiff{i,1,1,2}' paraCdiff{i,1,1,3}'];
       for j=1:size(parabydifficulty,1)
           plot([x(1) x(2)],parabydifficulty(j,:),'r','LineWidth',0.5);
           hold on;
       end
       for j=1:size(parabydifficultyC,1)
           plot([x(3) x(4)],parabydifficultyC(j,:),'k','LineWidth',0.5);
           hold on;
       end

       fige1=scatter(ones(size(parabydifficulty,1),1)*x(1),parabydifficulty(:,1),20,'r','filled');
       fige2=scatter(ones(size(parabydifficulty,1),1)*x(2),parabydifficulty(:,2),20,'r');
       figc1=scatter(ones(size(parabydifficultyC,1),1)*x(3),parabydifficultyC(:,1),20,'k','filled');
       figc2=scatter(ones(size(parabydifficultyC,1),1)*x(4),parabydifficultyC(:,2),20,'k');
       plot([0,3],[0,0],'--k','LineWidth',0.5);
       hold on;
       ylabel(ylabelstr{i});
       set(gca, 'XLim',[0 size(parabydelay,2)+1]);  
       ylim=get(gca,'YLim');
%        title(titlestr{i});
       [~,p]=ttest(parabydifficulty(:,1),parabydifficulty(:,2));
%        xtext=get(gca,'Xlim');
       ytext=get(gca,'Ylim');
       text(1,ytext(end),plabelsymbol(p));
       [~,p]=ttest(parabydifficultyC(:,1),parabydifficultyC(:,2));
       text(2,ytext(end),plabelsymbol(p));
       parabydifficultydiff=parabydifficulty(:,2)-parabydifficulty(:,1);
       parabydifficultydiffC=parabydifficultyC(:,2)-parabydifficultyC(:,1);
       [~,p2]=ttest2(parabydifficultydiff,parabydifficultydiffC);
%        [~,p2]=ranksum(parabydifficultydiff,parabydifficultydiffC);
       text(1.5,ytext(end)+0.1*(ytext(end)-ytext(1)),plabelsymbol(p2));
       plot([1,2],[ytext(end)+0.01*(ytext(end)-ytext(1)),ytext(end)+0.01*(ytext(end)-ytext(1))],'k-');
       set(gca,'YLim',[ylim(1), ytext(end)+0.1*(ytext(end)-ytext(1))]);%this way avoid ylim changing when changing fig size
       set(gca,'XTick',[1,2]);
       set(gca,'XTickLabel',{'hM4D ','mCherry '});
       %�ĳ���ת���label�Խ�ʡ�ռ�
       xtb = get(gca,'XTickLabel');% ��ȡ���������ǩ���
       xt = get(gca,'XTick');% ��ȡ��������̶Ⱦ��
       yt = get(gca,'YTick'); % ��ȡ��������̶Ⱦ��
       xtextp=xt;%ÿ����ǩ����λ�õĺ����꣬�����ȻӦ�ú�ԭ����һ���ˡ�
       ytextp=yt(1)*ones(1,length(xt)); % ������ʾ��ǩ��λ�ã�д����Ψһ��������ʵ����Ϊÿ����ǩ�ҷ���λ�õ�������
       % rotation��������ת�Ƕȴ�����ʱ����ת����ת�������HorizontalAlignment�������趨��
       % ��3������ֵ��left��right��center��������Ը�������ֵ���Լ�rotation��ĽǶȣ�����д����45
       text(xtextp,ytextp,xtb,'HorizontalAlignment','right','rotation',30,'fontsize',8);
       set(gca,'xticklabel',[]);% ��ԭ�еı�ǩ��ȥ
       if i==4
%        legend([fige1,fige2],'easy','hard','Location','best');
       leg=legend([figc1,figc2],'easy','hard','Location','best');
       set(leg,'position',[0.95 0.15 0.02 0.05]);
       legend boxoff;
       end
%        set(gca,'FontName','Arial','FontSize',14);
       box off;
    end
end
function [paraB,pparaB,paraBdiff, pparaBdiff,delaybylength]=fGetpara(dir, transform_method, vardelay,treatmentgroup, ndifficultygroup)
% paraB=cell(5,ndelaygroup+1,length(group),3);
%��һά��������ָ�꣨cor rate, violation rate, miss rate, RT, licking consistency index LCI)��
%�ڶ�ά�����ܡ��̡���delay������ndelaygroupȷ������;��vardelayΪ�Զ��廮�ֱ�׼����Ĭ�Ϸֳ�����
%����ά��������ʵ��������saline,CNO)��
%����ά���Ѷȣ�end/probe)
% pparaB=ones(size(paraB,1),size(paraB,2),size(paraB,3)/2,size(paraB,4));%�������paraB�еĲ���CNO-saline��pֵ
% paraBdiff=cell(5,ndelaygroup+1,length(group)/2,3);%����CNO-saline��Ч��
% pparaBdiff=ones(2,size(paraB,1));%����paraB�еĲ���CNO-saline��ֵ������delay/difficultyʱ��pֵ,ֻΪ��ͼ����������Ӧfigλ��,��һ�в����飬����Ƚϣ��ʴ˴���һά�����delay���ڶ�ά�����difficulty
% delayPart ����ndelaygroupȷ����ÿ�����ֵ�delay����
if length(vardelay)==1
    ndelaygroup=vardelay;
elseif length(vardelay)==4
    ndelaygroup=2;
    delaysetting=vardelay;
end
if strcmp(transform_method,'uni')
    paraB=cell(5,ndelaygroup+1,2*length(treatmentgroup),ndifficultygroup+1);%����ά����ipsi, contra����
    paraBdiff=cell(5,ndelaygroup+1,length(treatmentgroup),ndifficultygroup+1);%����CNO-saline��Ч��
else
    paraB=cell(5,ndelaygroup+1,length(treatmentgroup),ndifficultygroup+1);%��һά��������ָ�꣬�ڶ�ά�����ܡ��̡���delay������ά��������ʵ��������saline,CNO)������ά���Ѷȣ�end/probe)
    paraBdiff=cell(5,ndelaygroup+1,length(treatmentgroup)/2,ndifficultygroup+1);
end
for n=1:length(treatmentgroup)
    [animal_name,dataChoiceResult]=fFindChoiceMat(strcat(dir,treatmentgroup{n}));
    nsample=size(dataChoiceResult,1)-1;
    for i=1:nsample
        %delay��ȷ����
        if length(vardelay)==1
            [delaybylength,delay]=fDelayGroup(dataChoiceResult{i}(:,4),'ngroup',ndelaygroup);%���ܳɶ�ά���飬ÿһ�д���ͬ���ȵ�delay,% logical ��Ӿͱ��double
        elseif length(vardelay)==4
            [delaybylength,delay]=fDelayGroup(dataChoiceResult{i}(:,4),'user-defined',delaysetting);
        end
        %��ȡtrial results
        correct=dataChoiceResult{i}(:,2)==1;%ȡֵ1-4�ֱ����cor,err,miss.err
        error=dataChoiceResult{i}(:,2)==2;
        miss=dataChoiceResult{i}(:,2)==3;
        violation=dataChoiceResult{i}(:,2)==4;
        novio=correct+error;
        do=novio+violation;
        %��ȡendtrial/probetrial
        click=unique(dataChoiceResult{i}(:,1));
        sort(click);
        pend=logical((dataChoiceResult{i}(:,1)==click(1))+(dataChoiceResult{i}(:,1)==click(end)));
        pprobe=logical(ones(size(dataChoiceResult{i},1),1)-pend);
        difficulty=[pend+pprobe,pend,pprobe];%���ܳɶ�ά���飬ÿһ�д���ͬ�Ѷ�
        %��ȡreaction time
        indRT=double((dataChoiceResult{i}(:,5)>0).*(dataChoiceResult{i}(:,2)==1));%ֻ����ȷ��trial
        indLCI=double(dataChoiceResult{i}(:,2)==1);%ֻ����ȷ��trial
%         indLCI=double(dataChoiceResult{i}(:,2)==2);%ֻ�������trial
%         indLCI=double(dataChoiceResult{i}(:,2)~=3);%������trial����ȥ����ˮ��miss
        for k=1:size(difficulty,2)%�����Ѷ�
            %��ȡipsi/contra
            if strcmp(transform_method,'uni')
                dataChoiceResult{i}(:,1)=dataChoiceResult{i}(:,3);%���ҵ����Ļ����õ����е�������Ӧ�������滻��һ��Ƶ��
                protate=logical((dataChoiceResult{i}(:,6)==2)+(dataChoiceResult{i}(:,6)==3));%��Щtrials��Ҫ��ת��Ҳ�������£����Ҹ��ߵ�һ��
                dataChoiceResult{i}(protate,1)=7-dataChoiceResult{i}(protate,1);
                ipsi=(dataChoiceResult{i}(:,1)<=3);
                contra=(dataChoiceResult{i}(:,1)>=3);
                for j=1:size(delay,2)%���ֳ��̵�delay
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
                for j=1:size(delay,2)%���ֳ��̵�delay
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

pparaB=ones(size(paraB,1),size(paraB,2),size(paraB,3)/2,size(paraB,4));%�������paraB�еĲ���CNO-saline��pֵ
for i1=1:size(pparaB,1)%parameter
    for i2=1:size(pparaB,2)%delay
        for i4=1:size(pparaB,4)%difficulty
            paraBdiff{i1,i2,1,i4}=paraB{i1,i2,2,i4}-paraB{i1,i2,1,i4};%CNO-saline
            if size(paraB, 3)==2
                [~,p]=ttest(paraB{i1,i2,1,i4},paraB{i1,i2,2,i4});
                %                     p=round(p,4);%����С�����4λ
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
pparaBdiff=ones(2,size(paraB,1));%����paraB�еĲ���CNO-saline��ֵ������delay/difficultyʱ��pֵ,ֻΪ��ͼ����������Ӧfigλ��,��һ�в����飬����Ƚϣ��ʴ˴���һά�����delay���ڶ�ά�����difficulty
for i=1:size(paraB,1)
    [~,pbydelay]=ttest(paraBdiff{i,2,1,1},paraBdiff{i,end,1,1});
    pparaBdiff(1,i)=pbydelay;
    [~,pbydifficulty]=ttest(paraBdiff{i,1,1,2},paraBdiff{i,1,1,end});
    pparaBdiff(2,i)=pbydifficulty;
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
    pvalue=round(pvalue,2);%������λ����
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
    str='n.s.';
%     pvalue=round(pvalue,2);%������λ����
%     str=strcat('p=',num2str(pvalue));
end
end