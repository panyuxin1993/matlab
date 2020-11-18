dbstop if error;

close all;
clear;
filepath='F:\video tracking\M2 imaging video';
savepath=[filepath,filesep,'summary figures'];
summaryFile=[filepath,filesep,'imaging_video_data_summary.xlsx'];
[~,~,temp]=xlsread(summaryFile,1);
dataSummaryT=cell2table(temp(2:end,:));
dataSummaryT.Properties.VariableNames =temp(1,:);
nSession=size(dataSummaryT,1);
DLCiteration=2;

%% global settings
%coordinates of which body part
%for coordinates
% bodyparts={'Tongue','Tongue','LeftHandFingerTip','LeftHandFingerTip','RightHandFingerTip','RightHandFingerTip','Nose','Nose','LeftWhiskerTip','LeftWhiskerTip','RightWhiskerTip','RightWhiskerTip','LeftLickPort','LeftLickPort','RightLickPort','RightLickPort'};%{Tongue,LeftHandFingerTip,LeftHandFingerRoot,LeftLickPort,Nose,LeftWhiskerTip,LeftWhiskerRoot, etc.}
% coordinates={'x','y','x','y','x','y','x','y','x','y','x','y','x','y','x','y'};%{x,y,likelihood};
bodyparts={'Tongue','Tongue', 'LeftHandFingerTip','RightHandFingerTip','Nose','LeftLickPort', 'RightLickPort'};
coordinates={'x','y','x','x', 'x','x', 'x'};
% %for likelihood
% bodyparts={'Tongue','LeftHandFingerTip','RightHandFingerTip','LeftHandFingerRoot','RightHandFingerRoot','Nose','LeftWhiskerTip','RightWhiskerTip','LeftLickPort','RightLickPort'};
% coordinates={'likelihood','likelihood','likelihood','likelihood','likelihood','likelihood','likelihood','likelihood','likelihood','likelihood'};
%plot part from the pool
indBodyPartsChosen=1:length(bodyparts);
bodyparts=bodyparts(indBodyPartsChosen);
coordinates=coordinates(indBodyPartsChosen);
nbodyparts=length(bodyparts);
treshold4likelihood=0.5;
%align to behavior event
behEventAlign='stimOnset';%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
masklick='yes';
i_selectivity=4;%*********variable**************
selectivitystr={'stimuli','difficulty','sensory','choice'};%sensory means grouping difficulties;
trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
outcomeType='divideCorErr';%'combineCorErr','divideCorErr',
%plotting settings
if strcmp(behEventAlign,'stimOnset')
    frameNumTime=[0.5,2];%from 5s before align point to 5s after align point
else
    frameNumTime=[1,2];%from 5s before align point to 5s after align point
end
binsizeMovingTtest=3;%bin size for moving ttest
binstepMovingTtest=1;%bin step for moving ttest
pSigTtest=0.01;%significance level for moving ttest

meanEachSession=cell(nSession,nbodyparts);%each cell represent one session mean
meanDiffEachSession=cell(nSession,nbodyparts);%each cell represent one session mean
for iSession=1:nSession%can be a loop
led_file=[filepath,filesep,dataSummaryT.OLEDFileName{iSession},'.csv'];
behdata=[filepath,filesep,dataSummaryT.session{iSession},'_Virables.mat'];
switch DLCiteration
    case 1
        file_trace=[filepath,filesep,'iteration1',filesep,dataSummaryT.DLCFileName1{iSession},'.csv'];
    case 2
        file_trace=[filepath,filesep,'iteration2',filesep,dataSummaryT.DLCFileName2{iSession},'.csv'];
end

%% 导入数据
name_file_trace=strsplit(file_trace,'.');
if ~exist([name_file_trace{1},'.mat'],'file')
    [dcnum,dctxt,~] = xlsread(file_trace,1);%this process is time consuming, so not do very often
    save([name_file_trace{1},'.mat'],'dcnum','dctxt');
else
    load([name_file_trace{1},'.mat'])
end

%% from .beh file calculated a vector of trial start, for trial start from video use
fr=24;
load(behdata);
timeTrialStartBeh=cellfun(@fTrialStartStr2double,Data_extract.Time_trialStart);% time(s) relative to start of session day
name_file_led=strsplit(led_file,filesep);
name_file_led2=strsplit(name_file_led{end},'.');
if ~exist([filepath,filesep,name_file_led2{1},'.mat'],'file')
    T=readtable(led_file);%different sessions may need mannually set parameters and save result separately
    frameTrialStartVideo=fTrialStartByOLEDrefArduino( T, fr ,timeTrialStartBeh);
    save([filepath,filesep,name_file_led2{1},'.mat'],'frameTrialStartVideo');
else
    load([filepath,filesep,name_file_led2{1},'.mat'])
end
frameTrialStartVideo=reshape(frameTrialStartVideo,1,[]);
% figure;
% plot(frameTrialStartVideo);


%% get behavior event from .beh
ind_tr_1=1;%using data from trial 1
frameNum=double(round(frameNumTime*fr));
[behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( Data_extract, frameTrialStartVideo, 1000/fr ,ind_tr_1);%get behavior event time
[trialType,rule,trialTypeStr] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix','left',outcomeType);%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials

for iBodyPart=1:nbodyparts
    %% get coordinates of body parts, analogy to dff
    indcol1=cellfun(@(x) strcmp(bodyparts{iBodyPart},x),dctxt(2,:));
    indcol2=cellfun(@(x) strcmp(coordinates{iBodyPart},x),dctxt(3,:));
    indcol_likelihood=cellfun(@(x) strcmp('likelihood',x),dctxt(3,:));
    indcol_co=find(indcol1.*indcol2);
    indcol_like=find(indcol1.*indcol_likelihood);
    bodyco=dcnum(:,indcol_co);
    bodycoli=dcnum(:,indcol_like);           
    bodyco(bodycoli<treshold4likelihood)=nan;%rule out those low likelihood data
    if ~strcmp(coordinates{iBodyPart},'likelihood')
        bodyco=fBaselineCorrection(bodyco,round(20*fr));%40s as span
        bodyco=bodyco-nanmean(bodyco);%calculate pixel shift, if it is likelihood, no need for normalization
        % rule out those body parts with low likelihood
        if mean(bodycoli)<treshold4likelihood
            bodyco(:)=nan;
        end
    end
    bodyco=reshape(bodyco,1,[]);
    %% aligned the coordinates change of body parts to behavior events
    if strcmp(behEventAlign,'stimOnset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
        [ bodyco_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( bodyco, behEventFrameIndex,  frameNum );
    else
        [ bodyco_aligned, behEvent_aligned, licking_aligned ] = fAlignSigBehEvent( bodyco, behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
    end
    nfigcol=size(trialType,1)-2;%4 column(correct/error/miss/violation),companied with 4 lick raster
    meanEachSession{iSession,iBodyPart}=cell(size(trialType,2),nfigcol);
    meanDiffEachSession{iSession,iBodyPart}=cell(1,nfigcol);%diff between two groups mean
    for nStim=1:size(trialType,2) %for each stimulus
        for  nResult=1:nfigcol
            selectedTrialInd=trialType(nResult,nStim,:);
            selectedTrialInd=logical(squeeze(selectedTrialInd))';
            neuralActivity=bodyco_aligned(selectedTrialInd,:);
%             indPnan=sum(isnan(neuralActivity),2)/size(neuralActivity,2);
%             ind2nan=(indPnan>=0.5);%here threshold can be adjusted
%             neuralActivity(ind2nan,:)=nan;
            [neuralActivityMean, neuralActivityCI]=fMean_CI(neuralActivity,0.05);
            meanEachSession{iSession,iBodyPart}{nStim,nResult}=neuralActivityMean;
        end
    end
    for  nResult=1:nfigcol
        tempmeanEachSession1=cell2mat(meanEachSession{iSession,iBodyPart}(1:size(trialType,2)/2,nResult));
        tempmeanEachSession2=cell2mat(meanEachSession{iSession,iBodyPart}(size(trialType,2)/2+1:end,nResult));
        meanDiffEachSession{iSession,iBodyPart}{1,nResult}=nanmean(tempmeanEachSession1,1)-nanmean(tempmeanEachSession2,1);
    end
end
end

%% settings for plot
if strcmp(coordinates{iBodyPart},'likelihood')
    yrange=[ 0 , 1.1 ];
end
%plot color plot,4 column(correct/error/miss/violation)
nfigcol=size(trialType,1)-2;%4 column(correct/error/miss/violation),companied with 4 lick raster
figMeanTrace=figure;%plot mean trace
% set(gcf, 'position', [0 0 250*nfigcol 200*nbodyparts]);
set(gcf, 'PaperPosition', [0 0 2.5*nfigcol 2*nbodyparts]);
figMeanDiff=figure;%plot difference of mean trace for different choice etc. in one session
set(gcf, 'PaperPosition', [0 0 2.5*nfigcol 2*nbodyparts]);
for iBodyPart=1:nbodyparts
    %% plot coordinates change, ref plotDffPSTH
    if strcmp(outcomeType,'divideCorErr')
        titlestr={'Correct','Error','Miss','Violation'};
    else
        titlestr={'Do','Miss','Violation'};
    end
    titlestr=strcat(bodyparts{iBodyPart},'-',coordinates{iBodyPart},'-',titlestr);
    if size(trialType,2)==6
        color_mean_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};
    elseif size(trialType,2)==2
        color_mean_trace={[0 0 1],[1 0 0]};
        color_cases={[0.5 0.5 1],[1 0.5 0.5]};
    end
    if strcmp(rule,'low click rate-right')
        color_mean_trace=fliplr(color_mean_trace);
    end
    for nResult=1:nfigcol
        matCasesStimResult=cell(size(meanEachSession,1),size(trialType,2));
        for  nStim=1:size(trialType,2) %for each stimulus
            figure(figMeanTrace);%save mean trace
            subplot(nbodyparts,nfigcol,nResult+(iBodyPart-1)*nfigcol);
            for iSession=1:nSession %plot individual traces                
                meanEachSessionCase=meanEachSession{iSession,iBodyPart}{nStim,nResult};
                matCasesStimResult{iSession,nStim}=meanEachSessionCase;
                curve_case(nStim)=plot(1:size(meanEachSessionCase,2),meanEachSessionCase,'Color',color_cases{nStim},'linewidth',1);
                hold on;
            end       
            %plot mean trace
            figure(figMeanTrace);%save mean trace
            %subplot(size(trialType,2)+3,2*nfigcol,2*nResult-1+2*nfigcol*(size(trialType,2)+2));
            subplot(nbodyparts,nfigcol,nResult+(iBodyPart-1)*nfigcol);
            matCasesStimResultNow=cell2mat(matCasesStimResult(:,nStim));
            [neuralActivityMean,neuralActivityCI]=fMean_CI(matCasesStimResultNow,0.05);
            %plot CI
            if size(trialType,2)==2%if only two stimuli, plot CI, otherwise, do not plot for the reason of simplicity
                xpatch=[1:size(neuralActivity,2), fliplr(1:size(neuralActivity,2))];
                ypatch=[neuralActivityCI(1,:),fliplr(neuralActivityCI(2,:))];
                p=patch(xpatch,ypatch,color_mean_trace{nStim});%plot confidence interval
                p.FaceAlpha=0.1;
                p.EdgeColor=color_mean_trace{nStim};%'none';
                hold on;
            end
            %plot mean
            curve_meanTrace(nStim)=plot(1:size(neuralActivity,2),neuralActivityMean,'Color',color_mean_trace{nStim},'linewidth',2);
            hold on;
            title(titlestr{nResult});
        end
        if ~exist('yrange','var')
            y_lim=get(gca,'Ylim');
        else
            y_lim=yrange;
        end
        pTtest = fMovingTtest(cell2mat(matCasesStimResult(:,1)),cell2mat(matCasesStimResult(:,2)),binsizeMovingTtest,binstepMovingTtest,'ttest');
        yTtest=ones(size(pTtest))*y_lim(2)*0.9;
        xTtest=(1:length(pTtest));
        xTtest(pTtest>=pSigTtest)=nan;
        yTtest(pTtest>=pSigTtest)=nan;
        plot(xTtest,yTtest,'k-','LineWidth',1);
        %plot meandiff
        matMeanDiffEachSession=cell(size(meanDiffEachSession));
        for iSession=1:nSession
            figure(figMeanDiff);
            subplot(nbodyparts,nfigcol,nResult+(iBodyPart-1)*nfigcol);
            matMeanDiffEachSession{iSession,iBodyPart}=meanDiffEachSession{iSession,iBodyPart}{1,nResult};
            curve_Diffcase(nStim)=plot(1:size(matMeanDiffEachSession{iSession,iBodyPart},2),matMeanDiffEachSession{iSession,iBodyPart},'Color',[0.5,0.5,0.5],'linewidth',1);         
            hold on;
        end
        figure(figMeanDiff);
        subplot(nbodyparts,nfigcol,nResult+(iBodyPart-1)*nfigcol);
        matCasesStimResult=cell2mat(matMeanDiffEachSession(:,iBodyPart));
        [neuralActivityMean,neuralActivityCI]=fMean_CI(matCasesStimResult,0.05);
        %plot CI
        if size(trialType,2)==2%if only two stimuli, plot CI, otherwise, do not plot for the reason of simplicity
            xpatch=[1:size(neuralActivity,2), fliplr(1:size(neuralActivity,2))];
            ypatch=[neuralActivityCI(1,:),fliplr(neuralActivityCI(2,:))];
            p=patch(xpatch,ypatch,[0,0,0]);%plot confidence interval
            p.FaceAlpha=0.1;
            p.EdgeColor=[0,0,0];%'none';
            hold on;
        end
        %plot mean
        curve_meanDiff(nStim)=plot(1:size(neuralActivity,2),neuralActivityMean,'k-','linewidth',2);
        hold on;
        plot([1,size(neuralActivity,2)],[0,0],'k-','linewidth',1);
        title(titlestr{nResult});
        plot(xTtest,yTtest,'k-','LineWidth',1);
    end

end
for figid=1:2
for nResult=1:nfigcol
    for iBodyPart=1:nbodyparts
        figure(figid);
        subplot(nbodyparts,nfigcol,nResult+(iBodyPart-1)*nfigcol);
        if ~exist('yrange','var')
            y_lim=get(gca,'Ylim');
        else
            y_lim=yrange;
        end
        xlim([0,frameNumTime(2)+frameNumTime(1)]*fr);
        plot([round(frameNumTime(1)*fr),round(frameNumTime(1)*fr)]+1,[y_lim(1),y_lim(2)],'k-');%align to a behavior event
        hold on;
        if strcmp(behEventAlign,'stimOnset')
            plot(((round(frameNumTime(1)+0.5)*fr)+1)*ones(1,2),[y_lim(1),y_lim(2)],'k-');%align to a behavior event
        end
        if iBodyPart==nbodyparts
            xlabel(['time(s) from ',behEventAlign]);
        end
        if nResult==1
            if strcmp(coordinates{iBodyPart},'likelihood')
                ylabel('likelihood');
            else
                ylabel('\it\DeltaC');
            end
        end
        set(gca,'xtick',[round(fr*(frameNumTime(1)-floor(frameNumTime(1))))+1:round(fr):size(neuralActivity,2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
        %         set(gca,'xtick',[round(fr*(frameNumTime(1)-floor(frameNumTime(1))))+1:round(fr):size(neuralActivity,2)],'xticklabel',0:1:frameNumTime(2));
        set(gca,'FontName','Arial','FontSize',12);
        set(gca,'Ylim',y_lim);
    end
end
end

figure(figMeanTrace);
subplot(nbodyparts,nfigcol,1);
if contains(trialTypeStr,'stimuli')
    h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','','contra hard','','contra easy'},'Location','best');
    %legend(curve_meanTrace(:),num2str(Data_extract.Stimuli(:)),'Location','best');
elseif contains(trialTypeStr,'difficulty')
    h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','contra hard','contra easy'},'Location','best');
elseif contains(trialTypeStr,'first lick')
    h=legend(curve_meanTrace(:),{'ipsi lick first','contra lick first'},'Location','best');
elseif contains(trialTypeStr,'sensory')
    h=legend(curve_meanTrace(:),{'low click rate','high click rate'},'Location','best');
end
set(h,'box','off');
if strcmp(coordinates{iBodyPart},'likelihood')
    saveas(figMeanTrace,[savepath,filesep,'likelihood-algin to ',behEventAlign,'-sort ',behEventSort,'mean trace.pdf'],'pdf');
    saveas(figMeanDiff,[savepath,filesep,'Diff-likelihood-algin to ',behEventAlign,'-sort ',behEventSort,'mean trace.pdf'],'pdf');
else 
    saveas(figMeanTrace,[savepath,filesep,'coordinates-algin to ',behEventAlign,'-sort ',behEventSort,'mean trace.pdf'],'pdf');
    saveas(figMeanDiff,[savepath,filesep,'Diff-coordinates-algin to ',behEventAlign,'-sort ',behEventSort,'mean trace.pdf'],'pdf');
end

%% auxiliary fuction
function [timeTrialStart]=fTrialStartStr2double(strTrialStart)
%transform the trial start time from string form to double form
t=strsplit(strTrialStart,'_');
t=str2double(t);
timeTrialStart=t(4)*3600+t(5)*60+t(6)+t(7)/1000; %time(s) relative to trial1
end
function [out]=fBaselineCorrection(in,span)
%span = round(20*1000/frT/2); %every 40s, get a mean and substrate to compensate for baseline change
    x = [];
    for i = 1:length(in)
        %     ind_x = ind_x + 1;
        ind1 = max(i- span,1);
        ind2 = min(i+ span,length(in));
        x(i) = prctile(in(ind1:ind2),5);
    end
    out = in - x + nanmean(x);
end