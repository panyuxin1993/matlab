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
curve=cell(length(group),1);%用于存储图像句柄  
transform_method={};
ndelaygroup=2;
delaysetting=[300 900 900 1500];
ndifficultygroup=2;
[paraE,pparaE,paraEdiff,pparaEdiff,delaybylength]=fGetpara(dirE, transform_method, ndelaygroup ,group, ndifficultygroup);
[paraC,pparaC,paraCdiff,pparaCdiff,delaybylength]=fGetpara(dirC, transform_method, ndelaygroup ,group, ndifficultygroup);
fig=fCmpOverallPerformance(paraE,paraEdiff,paraC,paraCdiff,ndelaygroup);%画CNO-saline的差值，比较不同情况下的区别
set(fig,'PaperPosition',position);
end

function [fig] = fCmpOverallPerformance(paraE,paraEdiff,paraC,paraCdiff,ndelaygroup)
    fig=figure;%2*4张图，分别存correct rate, violation rate, miss rate, reaction time；第一行总结果，  第二行区分delay，不只是长短，而是分成200ms一组的6组,第三行区分difficulty
    set(gcf, 'position', [0 0 1400 800]);%控制fig尺寸
    ylabelstr={'\it\Delta\rmcorrect rate(%)';'\it\Delta\rmViolation rate(%)';'\it\Delta\rmMiss rate(%)';'\it\Delta\rmReaction time(ms)'};
    titlestr={'Correct';'Violation';'Miss';'Reaction time'};
%     a=suptitle('CNO-saline');
%     set(a,'FontSize',14,'position',[0.5 0.95 0.02 0.01]);
    pararaw=zeros(length(paraE{1,1,1,1}),4);%第一维表示样本数，值表示CNO-saline情况，第二维是四种指标
    pararawC=zeros(length(paraC{1,1,1,1}),4);
    x=[0.8,1.2,1.8,2.2];%各点的x坐标
    for i=1:4 %控制四列
        %首行不合并
        subplot(3,4,i);%首行合并
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
        plot([0,3],[0,0],'--k','LineWidth',0.5);%0处水平线
        %     for i=1:size(pararaw,1)%用线连接相同session的参数
        %         plot(pararaw(i,:));
        %         hold on;
        %     end
        hold on;
        %     set(gca,'Ylim',[-400,600]);
        set(gca, 'XLim',[0 3]);
        set(gca,'XTick',[1,2]);
        set(gca,'XTickLabel',{'hM4D ','mCherry '});
        %改成旋转后的label以节省空间
        xtb = get(gca,'XTickLabel');% 获取横坐标轴标签句柄
        xt = get(gca,'XTick');% 获取横坐标轴刻度句柄
        yt = get(gca,'YTick'); % 获取纵坐标轴刻度句柄
        xtextp=xt;%每个标签放置位置的横坐标，这个自然应该和原来的一样了。
        ytextp=yt(1)*ones(1,length(xt)); % 设置显示标签的位置，写法不唯一，这里其实是在为每个标签找放置位置的纵坐标
        % rotation，正的旋转角度代表逆时针旋转，旋转轴可以由HorizontalAlignment属性来设定，
        % 有3个属性值：left，right，center，这里可以改这三个值，以及rotation后的角度，这里写的是45
        text(xtextp,ytextp,xtb,'HorizontalAlignment','right','rotation',30,'fontsize',8);
        set(gca,'xticklabel',[]);% 将原有的标签隐去
        box off;
        
        subplot(3,4,4+i);%第二行,delay
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
        plot([0,3],[0,0],'--k','LineWidth',0.5);%0处水平线
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
        %改成旋转后的label以节省空间
        xtb = get(gca,'XTickLabel');% 获取横坐标轴标签句柄
        xt = get(gca,'XTick');% 获取横坐标轴刻度句柄
        yt = get(gca,'YTick'); % 获取纵坐标轴刻度句柄
        xtextp=xt;%每个标签放置位置的横坐标，这个自然应该和原来的一样了。
        ytextp=yt(1)*ones(1,length(xt)); % 设置显示标签的位置，写法不唯一，这里其实是在为每个标签找放置位置的纵坐标
        % rotation，正的旋转角度代表逆时针旋转，旋转轴可以由HorizontalAlignment属性来设定，
        % 有3个属性值：left，right，center，这里可以改这三个值，以及rotation后的角度，这里写的是45
        text(xtextp,ytextp,xtb,'HorizontalAlignment','right','rotation',30,'fontsize',8);
        set(gca,'xticklabel',[]);% 将原有的标签隐去
        if i==4
            legend([fige1,fige2],'short delay','long delay');
            leg=legend([figc1,figc2],'short delay','long delay');
            set(leg,'position',[0.95 0.5 0.02 0.05]);
            legend boxoff;
        end
        %        set(gca,'FontName','Arial','FontSize',14);
        box off;
       
       subplot(3,4,8+i);%第三行,difficulty
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
       %改成旋转后的label以节省空间
       xtb = get(gca,'XTickLabel');% 获取横坐标轴标签句柄
       xt = get(gca,'XTick');% 获取横坐标轴刻度句柄
       yt = get(gca,'YTick'); % 获取纵坐标轴刻度句柄
       xtextp=xt;%每个标签放置位置的横坐标，这个自然应该和原来的一样了。
       ytextp=yt(1)*ones(1,length(xt)); % 设置显示标签的位置，写法不唯一，这里其实是在为每个标签找放置位置的纵坐标
       % rotation，正的旋转角度代表逆时针旋转，旋转轴可以由HorizontalAlignment属性来设定，
       % 有3个属性值：left，right，center，这里可以改这三个值，以及rotation后的角度，这里写的是45
       text(xtextp,ytextp,xtb,'HorizontalAlignment','right','rotation',30,'fontsize',8);
       set(gca,'xticklabel',[]);% 将原有的标签隐去
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
    str='n.s.';
%     pvalue=round(pvalue,2);%保留两位即可
%     str=strcat('p=',num2str(pvalue));
end
end