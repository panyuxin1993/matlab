function [meanActivityByTrialType,cellActivityByTrialType] = fMeanPSTHbyROI(activities_data,ind_tr_1,Data_extract,SavedCaTrials,frameNumTime,behEventAlign,masklick,i_selectivity,trial2include,trial2exclude)
%FMEANPSTHBYROI get a cell matrix of mean activities grouped by ROI
%   Detailed explanation goes here
%Input-
% activity_form={'dff','spkr','zscored-spkr'},indicating the activity form


ntr = length(SavedCaTrials.f_raw); %use all data, or  just set the number of trials to use
frT = SavedCaTrials.FrameTime;
% align to behavior event
nFrameEachTrial=cellfun(@(x) size(x,2),SavedCaTrials.f_raw);
ind_1stFrame=zeros(1,length(nFrameEachTrial));
ind_1stFrame(1)=1;
ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
ind_1stFrame=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials
if isnan(trial2include)
    indTrial2include=fExcludeTrials(trial2exclude,ind_1stFrame,'logical');
    trial2include=[1,length(ind_1stFrame)];
elseif strcmp(trial2include,'all')
    trial2include=[1,length(ind_1stFrame)];
    indTrial2include=fIncludeTrials(trial2include,ind_1stFrame,'logical');
else
    indTrial2include=fIncludeTrials(trial2include,ind_1stFrame,'logical');
end
frameNum=double(round(frameNumTime*1000/frT));
[behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( Data_extract, ind_1stFrame, SavedCaTrials.FrameTime ,ind_tr_1);%get behavior event time

%   align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
if strcmp(behEventAlign,'delay onset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
    [ activities_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( activities_data, behEventFrameIndex,  frameNum );
else
    [ activities_aligned, behEvent_aligned, licking_aligned ] = fAlignSigBehEvent( activities_data, behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
end

[trialType,rule] = fGetTrialType( Data_extract,[],i_selectivity,'matrix','left','divideCorErr');%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
meanActivityByTrialType=cell(size(trialType,2),size(trialType,1));%1d-stimuli, 2d-(correct/error/miss/violation)
cellActivityByTrialType=cell(size(trialType,2),size(trialType,1));%1d-stimuli, 2d-(correct/error/miss/violation)
for nStim=1:size(trialType,2) %for each stimulus
    for nResult=1:size(trialType,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
        selectedTrialInd=trialType(nResult,nStim,:);
        selectedTrialInd=reshape(logical(squeeze(selectedTrialInd)),[],1);
        indTrial2include=reshape(indTrial2include,[],1);
        selectedTrialInd=logical(selectedTrialInd.*indTrial2include);
        neuralActivity=activities_aligned(selectedTrialInd,:);
        cellActivityByTrialType{nStim,nResult}=cell(1,1);
        cellActivityByTrialType{nStim,nResult}{1}=neuralActivity;
        meanActivityByTrialType{nStim,nResult}=nanmean(neuralActivity,1);
    end
end
end