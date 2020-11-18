clear;
close all;%close all figure
%  dir='D:\xulab\behavior\CNO experiment\CNO ctrl\infusion\bilateral\';
% dir='D:\xulab\behavior\CNO experiment\CNO ctrl\infusion\unilateral\';
% dir='D:\xulab\behavior\CNO experiment\CNO ctrl\ip\';
dire='D:\xulab\behavior\CNO experiment\infusion\bilateral\';
%  dir='D:\xulab\behavior\CNO experiment\infusion\unilateral\';
%dire='D:\xulab\behavior\CNO experiment\ip\';
group={'saline', 'CNO'};
color={'k','r'};
colorCase={[0.7,0.7,0.7],[1,0.7,0.7]};
timeMin=-500;%start from 0.5s before answer
binSize=200;%100ms
binStep=50;%moving window
trialLength=3000;%3s answer period + 3s opto stimulation
for n=1:length(group)
    fileFolder=strcat(dire,group{n});
    [animal_name,dataChoiceResult ] = fFindChoiceMat(fileFolder );
    dirmat=strcat(fileFolder,'\*.mat');
    dirs=dir(dirmat);
    dircell=struct2cell(dirs);
    filenames=dircell(1,:);
    sort(filenames);%将文件按照训练先后顺序排列起来
    lickTime=cell(length(filenames)+1,1);%最后一个cell存所有trial合起来
    lickRate=cell(length(filenames)+1,1);%最后一个cell存所有trial合起来
    for i=1:length(filenames)
        load(strcat(fileFolder,'\',filenames{i}));
        corTrialInd=(dataChoiceResult{i}(:,2)==1);%select 1correct trials 2 error trials, 3 miss trials, 4 violation trials
        [lickTime{i},maxdelay] =fGetLickTimes(SessionResults,'answer');
        [ lickMat,lickRate{i} ] = fLickRate( lickTime{i}, binSize,trialLength,timeMin ,corTrialInd,binStep);
    end
    %calculate licking rate mean, method 1, by pooling all trials together
    eachTrialNum=cellfun(@(x) length(x),lickTime);
    for i=1:length(filenames)
        trialInd(1,i)=sum(eachTrialNum(1:i));
    end
    trialInd=[0 trialInd];
    for i=1:length(trialInd)-2
        lickTime{end}(trialInd(i)+1:trialInd(i+1))=lickTime{i};
    end
    corTrialInd=(dataChoiceResult{end}(:,2)==1);%%select 1correct trials 2 error trials, 3 miss trials, 4 violation trials
    [ lickMat,lickRate{end} ] = fLickRate( lickTime{end}, binSize,trialLength,timeMin,corTrialInd,binStep);
    %calculate licking rate mean and ci, here is the mean of each sessions
    lickRateCell=cellfun(@(x) x(:,1)',lickRate,'UniformOutput',false);
    lickRateMat=cell2mat(lickRateCell);
    [lickRateMean, lickRateCI]=fMean_CI(lickRateMat,0.05);
    %plot licking rate
%     for i=1:length(filenames)
%         plot(timeMin:binStep:trialLength+timeMin-binSize+binStep,lickRate{i}(:,1),'Color',colorCase{n},'LineWidth',0.3);
%         hold on;
%     end
    %plot(timeMin:binStep:trialLength+timeMin-binSize+binStep,lickRate{end}(:,1),'Color',color{n},'LineWidth',2);
    plot(timeMin:binStep:trialLength+timeMin-binSize+binStep,lickRateMean,'Color',color{n},'LineWidth',3);
    hold on;
    plot(timeMin:binStep:trialLength+timeMin-binSize+binStep,lickRateCI(1,:),'Color',color{n},'LineWidth',0.5);
    hold on;
    plot(timeMin:binStep:trialLength+timeMin-binSize+binStep,lickRateCI(2,:),'Color',color{n},'LineWidth',0.5);
    hold on;
    xlabel('time(ms) from answer','FontName','Arial','FontSize',14);
    ylabel('lick rate(/s)','FontName','Arial','FontSize',14);
    set(gca,'FontName','Arial','FontSize',14);
end