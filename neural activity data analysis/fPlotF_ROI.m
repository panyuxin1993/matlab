function [ figOut ] = fPlotF_ROI( path, fdata,datatype,ind_ROI, ind_trial )
%FPLOTF_ROI Plot F raw or dF/F for one ROI, distinguish trial
%type (stimuli,choice, etc.). Colorplot
%Input-
%   fdata='raw'(default)|'dff', indicate which data to show
%   datatype= '2P'(default)|'FP', indicate which data structure
%   ind_ROI= a vector indicating the roi number to present
%   ind_trial= a vector indicating the trial index to present, only support
%   continues trials, e.g. 1:10, otherwise can't garantee bug free
dirmat=strcat(path,'\*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
if ~exist('Data_extract','var')
    file_beh=cellfun(@(x) contains(x,'_Virables'), filenames);
    if sum(file_beh)>0
        load([path,'\',filenames{1,file_beh}]);%load behavior data
    else
        %if no file then not plot behavior event
    end
end
if strcmp(datatype,'2P')
    file_behmat=cellfun(@(x) contains(x,'imaging.mat'),filenames);
    load([path,filesep,filenames{1,file_behmat}]);%load .mat file that transform form .beh
    if  ~exist('SavedCaTrials','var')
        file_imaging=cellfun(@(x) contains(x,'Ca'), filenames);
        load([path,filesep,filenames{1,file_imaging}]);%load imaging data
    end
    %plot f_raw and label properly
    frT = SavedCaTrials.FrameTime;
    if isempty(ind_trial)
        ind_trial=1:length(SavedCaTrials.f_raw);%if not enter range, just plot whole session
    end
    if strcmp(fdata,'raw')
        ylabelstr='f raw';
        f_data=[];
        for i=1:length(ind_trial)
            f_data=[f_data,SavedCaTrials.f_raw{ind_trial(i)}];
        end
        trialLengthTime=size(f_data,2)*frT/1000; 
    elseif strcmp(fdata,'dff')
        ylabelstr='\it\DeltaF/F';
        if exist([path,filesep,'dff.mat'],'file')
            temp=load([path,filesep,'dff.mat']);
            f_data=temp.dff;
        else
            f_data=fnx_getDff( path,path);
        end
        trialLengthTime=0;
        for i=1:length(ind_trial)
            trialLengthFrame_temp=size(SavedCaTrials.f_raw{ind_trial(i)},2);
            trialLengthTime=trialLengthTime+trialLengthFrame_temp*frT/1000;
        end
        ind_tr_1=1;
        ntr=length(SavedCaTrials.f_raw);
        nFrameEachTrial=cellfun(@(x) size(x,2),SavedCaTrials.f_raw);
        ind_1stFrame=zeros(1,length(nFrameEachTrial));
        ind_1stFrame(1)=1;
        ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
        indFrame_trial_start=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials
        if isempty(ind_trial)
            ind_trial=1:length(indFrame_trial_start);%if not enter range, just plot whole session
        end
%         trialLengthTime=indFrame_trial_start(min(max(ind_trial)+1,length(indFrame_trial_start)))-indFrame_trial_start(max(min(ind_trial),1));
        if max(ind_trial)>length(indFrame_trial_start)
            f_data=f_data(:,indFrame_trial_start(min(ind_trial)):end);
        else
            f_data=f_data(:,indFrame_trial_start(min(ind_trial)):(indFrame_trial_start(max(ind_trial)+1)));%baseline correction, fit 470 and 410, dff
        end
    end
elseif strcmp(datatype,'FP')
    if strcmp(fdata,'dff')
        ylabelstr='\it\DeltaF/F';
    end
    file_behmat=cellfun(@(x) contains(x,'FP.mat'),filenames);
    load([path,filesep,filenames{1,file_behmat}]);%load .mat file that transform form .beh
    file_FP=cellfun(@(x) contains(x,'dff_temp'), filenames);
    load([path,filesep,filenames{1,file_FP}]);%load imaging data

    FrameInfo = dir([path,filesep,'*.log']);
    fileID = fopen([path,filesep,FrameInfo.name]);
    C=textscan(fileID,'%d %d','HeaderLines',16);
    fclose(fileID);
    TrialCount=C{1,1};
    TrialStart_FrameCount=C{1,2};
    indFrame_trial_start=round(double(TrialStart_FrameCount')/2);
    if isempty(ind_trial)
        ind_trial=1:length(indFrame_trial_start);%if not enter range, just plot whole session
    end
    frT=50;%for FP, sample rate=20Hz
    trialLengthTime=(indFrame_trial_start(min(max(ind_trial)+1,length(indFrame_trial_start)))-indFrame_trial_start(max(min(ind_trial),1)))*frT/1000;
    if max(ind_trial)==length(indFrame_trial_start)
        f_data=dff{1}(:,indFrame_trial_start(min(ind_trial)):end);
    else
        f_data=dff{1}(:,indFrame_trial_start(min(ind_trial)):(indFrame_trial_start(max(ind_trial)+1)));%baseline correction, fit 470 and 410, dff
    end
end
if max(ind_ROI)>size(f_data,1) || min(ind_ROI)<1
    warning('index of ROI out of range');
end
%fmax=max(SavedCaTrials.f_raw{ind_trial}(ind_ROI,:));
figOut=figure;
clf;
if length(ind_ROI)==1 
    plot(f_data(ind_ROI,:),'k-');%plot f_raw
    hold on;
    if length(ind_trial)==1
        titlestr=strcat('ROI',num2str(ind_ROI),'-trial',num2str(ind_trial));
    else
        titlestr=strcat('ROI',num2str(ind_ROI),'-trial',num2str(ind_trial(1)),'-',num2str(ind_trial(end)));
    end
else%plot multiple neural traces, with no y tick labels, just pile up one by one
    range=[];
    for i=1:length(ind_ROI)
        current_f=f_data(ind_ROI(i),:);
        current_f=current_f-min(current_f)+sum(range);
        plot(current_f,'k-');
        hold on;
        xlim=get(gca,'Xlim');
        range=[range,max(current_f)-min(current_f)];
        text(xlim(end),sum(range),strcat('ROI-',num2str(ind_ROI(i))));
    end
    if length(ind_trial)==1
        titlestr=strcat(num2str(length(ind_ROI)),'ROIs','-trial',num2str(ind_trial));
    else
        titlestr=strcat(num2str(length(ind_ROI)),'ROIs','-trial',num2str(ind_trial(1)),'-',num2str(ind_trial(end)));
    end
    set(gca,'ytick',[],'yticklabel',[]);
    %plot scale bar and text
    if length(range)>=10
        scalebar_lengthy=fRound(prctile(range,10));
    else
        scalebar_lengthy=fRound(min(range)/2);
    end
    scalebar_lengthx=fRound(trialLengthTime/15);
    %n_digit=fix(log10(scalebar_lengthy));
    %scalebar_lengthy=round(scalebar_lengthy/10^n_digit)*10^n_digit;
    xlim=get(gca,'Xlim');
    plot([xlim(end),xlim(end)],[0,scalebar_lengthy],'k-','LineWidth',2);
    text(xlim(end),scalebar_lengthy/2,['\it\DeltaF/F ','\rm',num2str(scalebar_lengthy)],'FontName','Arial','FontSize',10);
    plot(xlim(end)-[0,scalebar_lengthx*round(1000/frT)],zeros(2,1),'k-','LineWidth',2);
    text(xlim(end)-scalebar_lengthx*round(1000/frT)/2,0,[num2str(scalebar_lengthx),'s'],'FontName','Arial','FontSize',10);
%     ylabelstr=strcat('f raw(',num2str(scalebar_length),')');
end       
    
ylim=get(gca,'Ylim');
%set(gca,'ytick',ylim(1):ylim(2),'yticklabel',sum(selectedTrialInd));
n_tick=10;
binTime=floor(trialLengthTime/n_tick);
set(gca,'xtick',[1:round(1000/frT)*binTime:trialLengthTime*1000/frT],'xticklabel',[0:binTime:trialLengthTime]);
%set(gca,'xtick',[],'xticklabel',[]);%for simplicity
xlabel('Time from trial start(s)');
ylabel(ylabelstr);
title(titlestr);
set(gca,'FontName','Arial','FontSize',14);
box off;

%plot behavior event
i_selectivity=4;%*********variable**************
selectivitystr={'stimuli','sensory difficulty','sensory','choice'};%sensory means grouping difficulties;
trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
if strcmp(datatype,'FP')
    fiberSide={'left','right'};
    fiberSide=fiberSide{ind_ROI};
elseif strcmp(datatype,'2P')
    fiberSide='left';
end
if exist('Data_extract','var') 
%     trialType=cellfun(@(x) double(x.Trial_Type), SessionResults);
    [trialType,behrule,trialTypeStr] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix',fiberSide);%NOTE HERE£¬ only one ROI for FP data,since they from different hemisphere
    trialType=squeeze(sum(trialType(1:2,:,:)));%only correct and error trials
    trialType=trialType(1,:)+2*trialType(2,:);%change to 1/2 values
    indFrame_trial_start=cellfun(@fTrialStartStr2double, SessionResults);
    relativeTrialStart=indFrame_trial_start(ind_trial)-indFrame_trial_start(ind_trial(1));%relative to the ind_trial(1) time
    relativeTrialStart=[relativeTrialStart, indFrame_trial_start(max(ind_trial)+1)-indFrame_trial_start(ind_trial(1))];
    for j=1:length(ind_trial)
        %plot bar indicating the trial number
        color_trial={[0,0,0],[0,0,1],[1,0,0]};%other, ipsi, contra
        plot([relativeTrialStart(j),relativeTrialStart(j+1)-2]*1000/frT,[ylim(2),ylim(2)],'color',color_trial{trialType(ind_trial(j))+1},'LineWidth',2);
%         %plot stim 
%         xpatch=relativeTrialStart(j)*1000/frT+[Data_extract.Stim_onset_time(1,ind_trial(j))/frT,Data_extract.Stim_offset_time(1,ind_trial(j))/frT,Data_extract.Stim_offset_time(1,ind_trial(j))/frT,Data_extract.Stim_onset_time(1,ind_trial(j))/frT];
%         ypatch=[ylim(1),ylim(1),ylim(2),ylim(2)];
%         p=patch(xpatch,ypatch,[0.5,0.5,0.5]);
%         p.FaceAlpha=0.3;
%         hold on;
%         if isfield(Data_extract,'Opto_trial_index') && Data_extract.Opto_trial_index(1,ind_trial(j))==1% for opto trials, plot opto-stimuli duration
%             xopto=relativeTrialStart(j)*1000/frT+[Data_extract.Opto_Onset_Time(1,ind_trial(j))/frT,Data_extract.Opto_Off_Time(1,ind_trial(j))/frT,Data_extract.Opto_Off_Time(1,ind_trial(j))/frT,Data_extract.Opto_Onset_Time(1,ind_trial(j))/frT];
%             p_opto=patch(xopto,ypatch,[0.5,0.5,1]);
%             p_opto.FaceAlpha=0.3;
%             hold on;
%         end
        % plot([Data_extract.Stim_onset_time(1,ind_trial(j))/frT,Data_extract.Stim_onset_time(1,ind_trial(j))/frT],[ylim(1) ylim(2)],'k');%stim onset time
        % hold on;
        % plot([Data_extract.Stim_offset_time(1,ind_trial(j))/frT,Data_extract.Stim_offset_time(1,ind_trial(j))/frT],[ylim(1) ylim(2)],'k');%stim offset time
        % hold on;
%         plot(relativeTrialStart(j)*1000/frT+[Data_extract.Go_time(1,ind_trial(j))/frT,Data_extract.Go_time(1,ind_trial(j))/frT],[ylim(1) ylim(2)],'Color',[0.5,0.5,0.5]);%go cue time
%         hold on;
%         for i=1:length(Data_extract.Left_lick_time{1,ind_trial(j)})
%             plot(relativeTrialStart(j)*1000/frT+[Data_extract.Left_lick_time{1,ind_trial(j)}(i)/frT,Data_extract.Left_lick_time{1,ind_trial(j)}(i)/frT],[ylim(1) ylim(2)],'b');%left lick time
%             hold on;
%         end
%         for i=1:length(Data_extract.Right_lick_time{1,ind_trial(j)})
%             plot(relativeTrialStart(j)*1000/frT+[Data_extract.Right_lick_time{1,ind_trial(j)}(i)/frT,Data_extract.Right_lick_time{1,ind_trial(j)}(i)/frT],[ylim(1) ylim(2)],'r');%right lick time
%             hold on;
%         end
%         Data_extract.Reward_time=double(Data_extract.Reward_time);
%         Data_extract.Reward_time(Data_extract.Reward_time==0)=nan;
%         plot(relativeTrialStart(j)*1000/frT+[Data_extract.Reward_time(1,ind_trial(j))/frT,Data_extract.Reward_time(1,ind_trial(j))/frT],[ylim(1) ylim(2)],'c');%reward time
%         hold on;
    end
end
end

%fuction for getting a round number, example: 166-->100, 25-->10
function [out]=fRound(in)
    n_digit=floor(log10(in));
    out=round(in/10^n_digit)*10^n_digit;
end

%transform the trial start time from string form to double form; but this
%method introduce errors since time of arduino is different from 2P
function [timeTrialStart]=fTrialStartStr2double(strTrialStartcell)
t=strsplit(strTrialStartcell.Time_trialStart,'_');
t=str2double(t);
timeTrialStart=t(4)*3600+t(5)*60+t(6)+t(7)/1000; %time(s) relative to trial1
end