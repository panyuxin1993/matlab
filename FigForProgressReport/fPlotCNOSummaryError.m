function [ fig ] = fPlotCNOSummaryError( fig,position,exp )
%FPLOTCNOSUMMARYERROR Summary of this function goes here
%   Detailed explanation goes here
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
ylabelstr={'\it\Delta\rmError(%)';'\it\Delta\rmViolation(%)';'\it\Delta\rmMiss(%)';'\it\Delta\rmReaction time(ms)'};
titlestr={'Correct';'Violation';'Miss';'Reaction time'};
animal_list={'pyx011','pyx012','pyx014','pyx015','pyx017','pyx019','pyx070','pyx072'};
animal_color={[0,0,0],[0,0,1],[0,1,0],[1,0,0],[1,0,1],[0,1,1],[1,1,0],[1,0.5,0]};

[paraE,pparaE,paraEdiff,pparaEdiff,delaybylength,animal_name_E]=fGetpara(dirE, transform_method, ndelaygroup ,group, ndifficultygroup);
[paraC,pparaC,paraCdiff,pparaCdiff,delaybylength,animal_name_C]=fGetpara(dirC, transform_method, ndelaygroup ,group, ndifficultygroup);
animal_ind=zeros(length(animal_name_E)-1,1);%1-8 means pyx011-pyx072
for i=1:length(animal_name_E)-1%last one is mean
    temp=cellfun(@(x) strcmp(x,animal_name_E{i}),animal_list);
    animal_ind(i)=find(temp);
end
derr=paraE{1,1,1,1}-paraE{1,1,2,1};
derr_delay=[paraE{1,2,1,1}-paraE{1,2,2,1};paraE{1,3,1,1}-paraE{1,3,2,1}];%2row, n_ample column
derr_difficulty=[paraE{1,1,1,2}-paraE{1,1,2,2};paraE{1,1,1,3}-paraE{1,1,2,3}];
%3*3 figure,1st row-derr, 2nd row-derr_delay, 3rd row-derr_difficulty;
% 1st colum-delta err by session. 2nd colum-delta err, showing different
% animal, 3rd column-delta err by animals
figure(fig);
subplot(3,3,1);%derr
fPlotScatter_Mean(derr);
title('by session');
set(gca,'XTick',[]);
subplot(3,3,2);%derr, different color showing different animal
fPlotScatter_Mean(derr,animal_ind,animal_color);
title('by session and animal');
set(gca,'XTick',[]);
subplot(3,3,4);%derr_delay
fPlotScatter_Mean(derr_delay);
set(gca,'XTick',[1,2],'XTickLabel',{'short','long'});
subplot(3,3,5);%derr_delay, different color showing different animal
fPlotScatter_Mean(derr_delay,animal_ind,animal_color);
set(gca,'XTick',[1,2],'XTickLabel',{'short','long'});
subplot(3,3,7);%derr_difficulty
fPlotScatter_Mean(derr_difficulty);
set(gca,'XTick',[1,2],'XTickLabel',{'easy','hard'});
subplot(3,3,8);%derr_difficulty, different color showing different animal
fPlotScatter_Mean(derr_difficulty,animal_ind,animal_color);
set(gca,'XTick',[1,2],'XTickLabel',{'easy','hard'});
%re-extract choices, grouped by animal
[paraE,pparaE,paraEdiff,pparaEdiff,delaybylength,animal_name_E]=fGetpara(dirE, transform_method, ndelaygroup ,group, ndifficultygroup,'animal');
derr=paraE{1,1,1,1}-paraE{1,1,2,1};
derr_delay=[paraE{1,2,1,1}-paraE{1,2,2,1};paraE{1,3,1,1}-paraE{1,3,2,1}];%2row, n_ample column
derr_difficulty=[paraE{1,1,1,2}-paraE{1,1,2,2};paraE{1,1,1,3}-paraE{1,1,2,3}];
figure(fig);
subplot(3,3,3);%derr
fPlotScatter_Mean(derr);
set(gca,'XTick',[]);
title('by animal');
subplot(3,3,6);%derr_delay
fPlotScatter_Mean(derr_delay);
set(gca,'XTick',[1,2],'XTickLabel',{'short','long'});
subplot(3,3,9);%derr_difficulty
fPlotScatter_Mean(derr_difficulty);
set(gca,'XTick',[1,2],'XTickLabel',{'easy','hard'});
set(fig,'PaperPosition',position);
end

function []=fPlotScatter_Mean(data,varargin)
color={[0,0,0],[0,0,1],[0,1,0],[1,0,0],[1,0,1],[0,1,1],[1,1,0],[1,0.5,0]};
x_mean=[0.8,2.2];%indicate where to put the mean
if nargin ==1
    animal_ind=ones(length(color),1);
    c='k';%default is all dots are black
else
    animal_ind=varargin{1};
    color=varargin{2};
    c=zeros(length(animal_ind),3);%create color for different group of dots
    for i=1:length(animal_ind)
        c(i,:)= color{animal_ind(i)};
    end
end

if size(data,1)==1 %just one condition
    scatter(ones(1,size(data,2))-0.2,data,20,c,'filled');
    hold on;
    mean_data=mean(data);
    sem=std(data)/sqrt(size(data,2));
    plot([0.8,1.2],[mean_data, mean_data],'-k');
    plot(1,mean_data,'ok');
    plot([0.9,1.1],[mean_data-sem, mean_data-sem],'-k');
    plot([0.9,1.1],[mean_data+sem, mean_data+sem],'-k');
    plot([0,2],[0,0],'k--');
    set(gca,'Xlim',[0.5,1.5]);
    [h,p]=ttest(data);
    set(gca,'ylim',[-10,20]);
    y=get(gca,'Ylim');
    text(1,y(end)*0.9,plabelsymbol(p));
elseif size(data,1)==2%should compare different condition, e.g. delay length/ difficulty
    for i=1:2%plot each group wiht dots
        scatter(ones(1,size(data,2))*i,data(i,:),20,c,'filled');
        hold on;
        mean_data=mean(data,2);%each row a group
        sem=std(data,0,2)/sqrt(size(data,2));%
        plot([x_mean(i)-0.2,x_mean(i)+0.2],[mean_data(i), mean_data(i)],'-k');
        plot(x_mean(i),mean_data(i),'ok');
        plot([x_mean(i)-0.1,x_mean(i)+0.1],[mean_data(i)-sem(i), mean_data(i)-sem(i)],'-k');
        plot([x_mean(i)-0.1,x_mean(i)+0.1],[mean_data(i)+sem(i), mean_data(i)+sem(i)],'-k');
        [h,p]=ttest(data(i,:));
        set(gca,'ylim',[-10,30]);
        y=get(gca,'Ylim');
        text(i,y(end)*0.8,plabelsymbol(p));
    end
    for i=1:size(data,2)%plot line linking two group points
        plot([1,2],data(:,i),'-','Color',[0.5,0.5,0.5],'LineWidth',0.5);
        hold on;
    end
    plot([0,3],[0,0],'k--');
    hold on;
    set(gca,'Xlim',[0.5,2.5]);
    [h,p2]=ttest(data(1,:),data(2,:));
    text(1.5,y(end),plabelsymbol(p2));
    plot([1,2],[y(end)*0.9,y(end)*0.9],'k-');
end
ylabel('\Delta error');

end

function [paraB,pparaB,paraBdiff, pparaBdiff,delaybylength,animal_name]=fGetpara(dir, transform_method, vardelay,treatmentgroup, ndifficultygroup,varargin)
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
if nargin>5%new variables indicate other usage
    mode=varargin{1};
else
    mode='session';
end
for n=1:length(treatmentgroup)
    [animal_name,dataChoiceResult]=fFindChoiceMat(strcat(dir,treatmentgroup{n}),mode);
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
elseif pvalue<0.1 && pvalue>=0.05
%     str='n.s.';
    pvalue=round(pvalue,2);%������λ����
    str=strcat('p=',num2str(pvalue));
else
    pvalue=round(pvalue,1);%����1λ����
    str=strcat('p=',num2str(pvalue));
end
end