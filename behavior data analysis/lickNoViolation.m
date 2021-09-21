%plot lick density during delay when no violation(allow violation)
%get lick time
clear;
date='2021_09_16';
animal='pyx378';
fileFolder=strcat('E:\xulab\behavior\',animal);
% session=strcat('\',date,'_',animal,'.mat');
% dirmat=strcat(fileFolder,session);
dirmat=strcat(fileFolder,'\*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
sort(filenames);
lickTime=cell(length(filenames),1);%单个cell表示一个session的所有lick

ntotal=length(filenames);
nfig=1;
divideTrialType=1;%2则一个session区分左右trialType而画成两部分,否则画在一  起
resample='resample';%[];%'resample';%trial数量太多的话可以可考虑resample
rasterSize=2;%每个lick点的横向长度，因为分辨率为1ms小鼠不可能舔水这么密集，此值甚至可以到上百
align='delayOff';%align={'stim','go','delayOff'};
npast=1;%how many sessions will be analyzed
figure;
if nfig==1
    set(gcf, 'position', [100 100 600 600]);%控制fig尺寸
end
for n=ntotal-nfig+1:ntotal
    load(strcat(fileFolder,'\',filenames{n}));
    [lickTime{n},maxdelay]=fGetLickTimes(SessionResults,align);
    %排序
    if divideTrialType==1%一共一行fig，每个session一张fig
        subplot(2,nfig+1,n-ntotal+nfig);%这种设定能画指定数量张图，无violation的1000ms和首次加上300ms的violation
        strtitle=fPlotLickRasterOneSession(lickTime{n},[0 1],'withVio','raw',filenames{n},animal,maxdelay,resample,rasterSize,align);
        subplot(2,nfig+1,n-ntotal+nfig*2+1);%第二行画licking rate
        fPlotLickRateOneSession(lickTime{n},strtitle,maxdelay,n,fileFolder,align);
    else
        subplot(2,nfig+1,n-ntotal+nfig);%第一行画左边的trial
        fPlotLickRasterOneSession(lickTime{n},0,'withVio','sort',filenames{n},animal,maxdelay,resample,rasterSize);
        subplot(2,nfig+1,n-ntotal+nfig*2);%第二行画右边的trial
        fPlotLickRasterOneSession(lickTime{n},1,'withVio','sort',filenames{n},animal,maxdelay,resample,rasterSize);
    end
end
[pEarlyLickRateIndex,pTrialEarlyLick]=fAnalyzeEarlyLick(npast,fileFolder,align,ntotal);
subplot(2,nfig+1,n-ntotal+nfig+1);
plot(pTrialEarlyLick(:,1),'k-');hold on;
plot(pTrialEarlyLick(:,2),'g-');
plot(pTrialEarlyLick(:,3),'r-');
xlabel('session');
ylabel('percent of early lick trial');
set(gca,'FontSize',14,'FontName','Arial');
legend('all trial','correct trial','error trial');
box off;
legend boxoff;
subplot(2,nfig+1,n-ntotal+nfig*2+2);
plot(pEarlyLickRateIndex(:,1),'k-');hold on;
plot(pEarlyLickRateIndex(:,2),'g-');
plot(pEarlyLickRateIndex(:,3),'r-');
xlabel('session');
ylabel('early lick rate index');
set(gca,'FontSize',14,'FontName','Arial');
legend('all trial','correct trial','error trial');
box off;
legend boxoff;
% fSavePPT(strcat(fileFolder,'\test.ppt'),filenames{n},'-f1');

function [strtitle]=fPlotLickRasterOneSession(lickTime,trialType,resultType,sortType,filename,animal,maxdelay,nresample,rasterSize,align)
%trialType选择trial为左0右1,均选择则用[0 1];ismember(lickTime.trialType,trialType)
%resultType选择result是‘withVio'所有trial，’noVio'不画violation
%sortType选择是否sort，'sort'表示sort，'raw'表示不sort
if strcmp(nresample,'resample')
    nsample=unidrnd(length(lickTime),1,100);%trial 数量太多，抽取其中一部分
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
set(gca, 'XLim',[-maxdelay maxdelay+2000]);
xlabel(['time(ms) from ',align],'FontName','Arial','FontSize',14);
% xlabel('time(ms) from go cue','FontName','Arial','FontSize',14);
ylabel('trials','FontName','Arial','FontSize',14);
% title(strtitle,'FontName','Arial','FontSize',14);
plot([0,0],[0,1000],'k','linewidth',1);
hold on;
set(gca,'FontName','Arial','FontSize',14);
end
function [strtitle]=fPlotLickRateOneSession(lickTime,strtitle,maxdelay,n,fileFolder,align)
%input- n-index of session; 
timeMin=-maxdelay;%start from () before aligned event, consistent with licking raster
binSize=200;%100ms
binStep=50;%moving window
trialLength=2*maxdelay+2000;
[~,dataChoiceResult ] = fFindChoiceMat(fileFolder );
color={[0,1,0];[1,0,0];[1,0.8,0]};%color for cor,err,vio

corTrialInd=(dataChoiceResult{n}(:,2)==1);%select 1correct trials 2 error trials, 3 miss trials, 4 violation trials
[ ~,lickRate_cor ] = fLickRate( lickTime, binSize,trialLength,timeMin, corTrialInd, binStep);
errTrialInd=(dataChoiceResult{n}(:,2)==2);%select 1correct trials 2 error trials, 3 miss trials, 4 violation trials
[ ~,lickRate_err ] = fLickRate( lickTime, binSize,trialLength,timeMin, errTrialInd, binStep);
vioTrialInd=(dataChoiceResult{n}(:,2)==4);%select 1correct trials 2 error trials, 3 miss trials, 4 violation trials
[ ~,lickRate_vio ] = fLickRate( lickTime, binSize,trialLength,timeMin, vioTrialInd, binStep);
xts=timeMin:binStep:trialLength+timeMin-binSize+binStep;
xts=double(xts);
ind_align=find(xts>=0);
if length(xts)~=length(lickRate_cor)
    n_ts=min(length(xts),length(lickRate_cor));
    plot(xts(1:n_ts),lickRate_cor(1:n_ts,1),'Color',color{1},'LineWidth',2);%green for cor
    hold on;
    plot(xts(1:n_ts),lickRate_err(1:n_ts,1),'Color',color{2},'LineWidth',2);%red for err
    plot(xts(1:n_ts),lickRate_vio(1:n_ts,1),'Color',color{3},'LineWidth',2);%orange for vio
    box off;
else
    plot(xts,lickRate_cor(:,1),'Color',color{1},'LineWidth',2);%green for cor
    hold on;
    plot(xts,lickRate_err(:,1),'Color',color{2},'LineWidth',2);%red for err
    plot(xts,lickRate_vio(:,1),'Color',color{3},'LineWidth',2);%orange for vio
    box off;
end
y_text=get(gca,'Ylim');
plot([0,0],[y_text(1),y_text(end)],'k-');
xlabel(['time(ms) from ',align],'FontName','Arial','FontSize',14);
ylabel('licks/s','FontName','Arial','FontSize',14);
set(gca, 'XLim',[-maxdelay maxdelay+2000]);
set(gca,'FontName','Arial','FontSize',14);
dy_text=y_text(end)-y_text(1);
dxts=xts(end)-xts(1);
pLick_cor=sum(lickRate_cor(1:ind_align(1),1))/sum(lickRate_cor(:,1));
pLick_err=sum(lickRate_err(1:ind_align(1),1))/sum(lickRate_err(:,1));
pLick_vio=sum(lickRate_vio(1:ind_align(1),1))/sum(lickRate_vio(:,1));
% text(xts(end)-0.1*dxts,y_text(1)+0.9*dy_text,['p=',num2str(round(pLick_cor,2))],'Color',color{1});
% text(xts(end)-0.1*dxts,y_text(1)+0.8*dy_text,['p=',num2str(round(pLick_err,2))],'Color',color{2});
% text(xts(end)-0.1*dxts,y_text(1)+0.7*dy_text,['p=',num2str(round(pLick_vio,2))],'Color',color{3});
end
function [pEarlyLickRateIndex,pTrialEarlyLick]=fAnalyzeEarlyLick(npast,fileFolder,align,varargin)
%   input:
%   npast-set how many sessions data to show, 
%   fileFolder-set where data from, 
%   varargin-set which session to be show in detail(one subplot show changes 
%       within that session),if it is string, it decide filename of chosen
%       session, if it is num, it set index of chosen session
%output-1 early lick rate index,which is the ratio of lickrate before 0(the aligned event) versus that after0; 
%2-percent of trials that lick early
dirmat=strcat(fileFolder,'\*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
sort(filenames);
if isnumeric(varargin{1})
    i_session=varargin{1};
elseif ischar(varargin{1})
    indSession=cellfun(@(x) strcmp(x,varargin{1}),filenames);
    i_session=find(indSession);
end
if nargin==1
    i_session=length(filenames);%default is to plot the last session
end
pEarlyLickRateIndex=zeros(npast,3);
pTrialEarlyLick=zeros(npast,3);
if i_session-npast+1<=0%如果此值大于总数，则将所有的数据画出来；小于总数则画出指定数量
    istart=1;
else
    istart=i_session-npast+1;
end
[~,dataChoiceResult ] = fFindChoiceMat(fileFolder );
for isession=istart:i_session
    load(strcat(fileFolder,'\',filenames{isession}));
    [lickTime,maxdelay]=fGetLickTimes(SessionResults,align);
    earlyLickCount_L=arrayfun(@(x) sum((x.lickleft>x.stimOn).*(x.lickleft<=x.go)),lickTime);%count of lick between stim_onset and go cue
    earlyLickCount_R=arrayfun(@(x) sum((x.lickright>x.stimOn).*(x.lickright<=x.go)),lickTime);
    earlyLickCount=earlyLickCount_L+earlyLickCount_R;
    earlyLickInd=(earlyLickCount>0);
    corTrialInd=(dataChoiceResult{isession}(:,2)==1);%dataChoiceResult{n}(:,2) 1correct trials 2 error trials, 3 miss trials, 4 violation trials
    errTrialInd=(dataChoiceResult{isession}(:,2)==2);
    pTrialEarlyLick_all=sum(earlyLickInd)/length(earlyLickInd);
    pTrialEarlyLick_cor=sum(earlyLickInd.*corTrialInd')/sum(corTrialInd);
    pTrialEarlyLick_err=sum(earlyLickInd.*errTrialInd')/sum(errTrialInd);
    pTrialEarlyLick(isession-istart+1,:)=[pTrialEarlyLick_all,pTrialEarlyLick_cor,pTrialEarlyLick_err];
    
    lickRate_early_L=arrayfun(@(x) double(sum((x.lickleft>x.stimOn).*(x.lickleft<=x.go)))/double(x.go-x.stimOn)*1000,lickTime);%licks/s;early licks are from stim_onset to go cue
    lickRate_early_R=arrayfun(@(x) double(sum((x.lickright>x.stimOn).*(x.lickright<=x.go)))/double(x.go-x.stimOn)*1000,lickTime);%licks/s;late llicks are from go cue to end of trial(estimated)
    lickRate_late_L=arrayfun(@(x) double(sum(x.lickleft>x.go))/double((max(max(max(x.lickleft),max(x.lickright)),x.answer+2000))-x.go)*1000,lickTime);
    lickRate_late_R=arrayfun(@(x) double(sum(x.lickright>x.go))/double((max(max(max(x.lickleft),max(x.lickright)),x.answer+2000))-x.go)*1000,lickTime);
    pEarlyLickRate=(lickRate_early_L+lickRate_early_R)./(lickRate_early_L+lickRate_early_R+lickRate_late_L+lickRate_late_R);
    pEarlyLickRateIndex_all=nansum(pEarlyLickRate)/length(pEarlyLickRate);
    pEarlyLickRateIndex_cor=nansum(double(pEarlyLickRate).*double(corTrialInd'))/sum(corTrialInd);
    pEarlyLickRateIndex_err=nansum(double(pEarlyLickRate).*double(errTrialInd'))/sum(errTrialInd);
    pEarlyLickRateIndex(isession-istart+1,:)=[pEarlyLickRateIndex_all,pEarlyLickRateIndex_cor,pEarlyLickRateIndex_err];
end
end