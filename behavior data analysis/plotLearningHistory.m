%from a xlsx file get animal names, then get learning history data from mat
%files, and write results back to the xlsx
rootpath='D:\xulab\behavior';
tableFile=strcat(rootpath,'\learning_history.xlsx');
[num,txt,raw] = xlsread(tableFile,1);
ncol=size(txt,2);
n_animal=ncol-1;
figLH=figure;
set(gcf,'position',[10,10,700,250]);
for i=1:n_animal
    wd=strcat(rootpath,'\',txt{1,i+1});
    [animal_name,dataChoiceResult] = fFindChoiceMat(wd);
    [animal_name_check,learningHistory]=fGetLearningHistory(wd,'txt');
    if strcmp(animal_name{1},animal_name_check)
        xlRange=strcat(char(65+i),num2str(2));%write to the xls from 2th row, char(65+i)th column;char(65) is'A', which is first column, used as training day
        xlswrite(strcat(rootpath,'\learning_history.xlsx'),learningHistory,1,xlRange);
        ind_2AFC=find(cellfun(@(x) strcmp(x,'2AFC'),learningHistory));
        ind_2AFC=fExtractConsecutiveVector(ind_2AFC','first');
        ind_retract=find(cellfun(@(x) strcmp(x,'retract'),learningHistory));
        ind_retract=fExtractConsecutiveVector(ind_retract','first');
        if ~isempty(ind_2AFC) && ~isempty(ind_retract) &&ind_2AFC(1)<ind_retract(1)%this animal undergo regular 2AFC training 
            color=[0,0,0];
            ind_chosen=ind_2AFC;
        elseif isempty(ind_retract)
            color=[0,0,0];
            ind_chosen=ind_2AFC;
        else%this animal learn 2AFC with lickport retracted
            color=[1,0,0];
            ind_chosen=ind_retract;
        end
        figure(figLH);
        subplot(1,3,1);%plot performace as a function of training day
        cor_rate=cellfun(@(x) sum(x(:,2)==1)/(sum(x(:,2)==1)+sum(x(:,2)==2)),dataChoiceResult);
        plot(cor_rate(ind_chosen(1):ind_chosen(end)),'color',color);hold on;
        ylabel('correct rate');
        xlabel('training day');
        set(gca,'FontName','Arial','FontSize',14);
        subplot(1,3,2);%plot LCI as a function of training day
        LCI = fGetLCI(length(ind_chosen),[wd,'\'],ind_chosen(end));
        plot(LCI(:,1),'color',color);hold on;
        ylabel('licking consistency index');
        xlabel('training day');
        set(gca,'FontName','Arial','FontSize',14);
        subplot(1,3,3);%plot licking rate as a function of training day
        dirmat=strcat(wd,'\*.mat');
        dirs=dir(dirmat);
        dircell=struct2cell(dirs);
        filenames=dircell(1,:);
        sort(filenames);%将文件按照训练先后顺序排列起来
        load(strcat(wd,'\',filenames{1}));%plot licking rate of first day
        [lickTime,maxdelay] =fGetLickTimes(SessionResults,'answer');
        timeMin=-500;%start from 0.5s before answer
        binSize=200;%100ms
        binStep=50;%moving window
        trialLength=4000;%3s answer period + 3s opto stimulation
        trialInd=ones(length(lickTime),1);
        [ lickMat,lickRate ] = fLickRate( lickTime, binSize,trialLength,timeMin, trialInd, binStep);
        plot(timeMin:binStep:trialLength+timeMin-binSize+binStep,lickRate(:,1),'color',color);hold on;
        ylabel('lick rate on first day');
        xlabel('time from answer(ms)');
        set(gca,'FontName','Arial','FontSize',14);
    else
        warning(strcat('Inconsistent animal name: ', animal_name{1},'&',animal_name_check));
    end
end