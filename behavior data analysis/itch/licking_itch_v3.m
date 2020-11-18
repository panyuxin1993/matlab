%plot licking raster and compare opto/non-opto trials for single session
%get lick time
clear;
animal='';
maxdelay=0;
fileFolder='D:\xulab\behavior\practice\LTY\LTY_ChR2_01\opto_effect';
%fileFolder='D:\xulab\behavior\CNO experiment\infusion\bilateral\saline';
% session=strcat('\',date,'_',animal,'.mat');
% dirmat=strcat(fileFolder,session);
dirmat=strcat(fileFolder,'\*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
sort(filenames);
figure;
% set(gcf, 'position', [0 0 1000 400]);%控制fig尺寸 
ntotal=1;%length(filenames);
nfig=ntotal;
divideTrialType=2;%2则一个session区分左右trialType而画成两部分,否则画在一起
resample=[];%'resample';%trial数量太多的话可以可考虑resample
rasterSize=2;%每个lick点的横向长度，因为分辨率为1ms小鼠不可能舔水这么密集，此值甚至可以到上百
for n=ntotal-nfig+1:ntotal
    load(strcat(fileFolder,'\',filenames{n}));
    %[lickTime{n},maxdelay]=fGetLickTimes(SessionResults,'go');%align={'stim','go'};
    %get a 'lickTime' variable, in default, align to go cue
    lickTime=struct('trialType',[],'lickleft',[],'lickright',[],'choice',[],'go',[]);
    optoType=double(cellfun(@(x) x.Time_optoStimOnset~=0, SessionResults));%1opto,0ctrl
    for i=1:length(SessionResults)
        base=SessionResults{i}.Time_delayOffset;  
        lickTime(i).go=0;
        lickTime(i).trialType=SessionResults{i}.Trial_Type;%左0右1
        left=strsplit(SessionResults{i}.Action_lickTimeLeft,'|');
        left(:,1)=[];
        lickTime(i).lickleft=str2double(left)-double(base);
        right=strsplit(SessionResults{i}.Action_lickTimeRight,'|');
        right(:,1)=[];
        lickTime(i).lickright=str2double(right)-double(base);
        lickTime(i).choice=double(SessionResults{i}.Action_choice);
        %lickTime(i).optoType=double((lickTime(i).choice)==(SessionResults{i}.Time_optoStimOnset~=0));%0左1右
    end
    %筛选opto/non-opto trials
    baseline=lickTime(1:200);
    experiment=lickTime(201:end);
    indOpto=(optoType(201:end)==1);
    ctrl=experiment(~indOpto);
    opto=experiment(indOpto);
    %licking rate
    binSize=100;%100ms
    trialLength=3000;%3s answer period + 3s opto stimulation
    sampleNum=trialLength/binSize+1;
    timeMin=-500;%start from 1s before go cue
    lickBaseline=zeros(length(baseline),trialLength,2);%第三维分表示左/右边舔水率
    lickCtrl=zeros(length(ctrl),trialLength,2);
    lickOpto=zeros(length(opto),trialLength,2);
    lickRate=zeros(3,sampleNum,2);%第一维分别是baseline, ctrl, opto，第三维是左右
    for i=1:length(baseline)
        for j=1:length(baseline(i).lickleft)
            if ~isnan(baseline(i).lickleft(j)) && baseline(i).lickleft(j)-timeMin<size(lickBaseline,2)
                lickBaseline(i,baseline(i).lickleft(j)-timeMin,1)=1;
            end
        end
        for j=1:length(baseline(i).lickright)
            if ~isnan(baseline(i).lickright(j)) && baseline(i).lickright(j)-timeMin<size(lickBaseline,2)
                lickBaseline(i,baseline(i).lickright(j)-timeMin,2)=1;
            end
        end
    end
    for i=1:length(ctrl)
        for j=1:length(ctrl(i).lickleft)
            if ~isnan(ctrl(i).lickleft(j)) && ctrl(i).lickleft(j)-timeMin<size(lickCtrl,2)
                lickCtrl(i,ctrl(i).lickleft(j)-timeMin,1)=1;
            end
        end
        for j=1:length(ctrl(i).lickright)
            if ~isnan(ctrl(i).lickright(j)) && ctrl(i).lickright(j)-timeMin<size(lickCtrl,2)
                lickCtrl(i,ctrl(i).lickright(j)-timeMin,2)=1;
            end
        end
    end
    for i=1:length(opto)
        for j=1:length(opto(i).lickleft)
            if ~isnan(opto(i).lickleft(j)) && opto(i).lickleft(j)-timeMin<size(lickOpto,2)
                lickOpto(i,opto(i).lickleft(j)-timeMin,1)=1;
            end
        end
        for j=1:length(opto(i).lickright)
            if ~isnan(opto(i).lickright(j)) && opto(i).lickright(j)-timeMin<size(lickOpto,2)
                lickOpto(i,opto(i).lickright(j)-timeMin,2)=1;
            end
        end
    end
    for i=1:2
        for j=1:sampleNum-1
            lickRate(1,j,i)=sum(sum(lickBaseline(:,(100*j-99):(100*j),i)))*10/length(baseline);
            lickRate(2,j,i)=sum(sum(lickCtrl(:,100*j-99:100*j,i)))*10/length(ctrl);
            lickRate(3,j,i)=sum(sum(lickOpto(:,100*j-99:100*j,i)))*10/length(opto);
        end
    end
    %一共2行fig，每个session一张fig一列，第一行画ctrl,第二行画opto
    %set(gcf,'position',[0,0,200*nfig,400]);
    subplot(4,nfig,n-ntotal+nfig);%baseline
    fPlotLickRasterOneSession(baseline,[0 1 2],'withVio','raw',filenames{n},animal,maxdelay,resample,rasterSize);
    title('baseline');
    set(gca, 'XLim',[timeMin trialLength+timeMin]);
    subplot(4,nfig,n-ntotal+nfig*2);%ctrl
    fPlotLickRasterOneSession(ctrl,[0 1 2],'withVio','raw',filenames{n},animal,maxdelay,resample,rasterSize);
    title('ctrl');
    set(gca, 'XLim',[timeMin trialLength+timeMin]);
    subplot(4,nfig,n-ntotal+nfig*3);%opto
    fPlotLickRasterOneSession(opto,[0 1 2],'withVio','raw',filenames{n},animal,maxdelay,resample,rasterSize);
    title('opto');
    set(gca, 'XLim',[timeMin trialLength+timeMin]);
    subplot(4,nfig,n-ntotal+nfig*4);%licking rate, three curves in a panel-->left
    color={[0,0,0],[0.5,0.5,0.5],[1,0,0]};
    for i=1:3
        plot(timeMin:binSize:trialLength+timeMin,lickRate(i,:,1)+lickRate(i,:,2),'Color',color{i},'LineWidth',2);
        hold on;
    end
    if n==ntotal
        legend('baseline','ctrl','opto');
    end
    title('lick rate');
    xlabel('time(ms) from go cue','FontName','Arial','FontSize',14);
    ylabel('licks/s');
    set(gca,'FontName','Arial','FontSize',14);
    set(gca, 'XLim',[timeMin trialLength+timeMin]);
    

%     subplot(5,nfig,n-ntotal+nfig*5);%licking rate, three curves in a panel-->right
%     for i=1:3
%         plot(timeMin:binSize:trialLength+timeMin,lickRate(i,:,2),'Color',color{i},'LineWidth',2);
%         hold on;
%     end
%     if n==ntotal
%         legend('baseline','ctrl','opto');
%     end
%     title('lick rate of right');
%     set(gca,'FontName','Arial','FontSize',14);
%     set(gca, 'XLim',[timeMin trialLength+timeMin]);
end