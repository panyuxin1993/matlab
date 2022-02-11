function [ behEventFrameIndex, lickingFrameIndex ] = fGetBehEventTime( behData, ind_1stFrame, frameTime ,varargin)
%FGETBEHEVENTTIME Get behavior event time,transform to frame index in a
% whole-session-long trial
%   Detailed explanation 
%   Input is behavior data(struct) and 1st_frame index of each trials(1*n
%   matrix) in the whole session, frameTime(a number);
%   Output is the frame number index of different behavior event, such as
%   stimuli onset, stimuli offset, go cue onset, first lick time, answer
%   time, reward time,
if isempty(varargin) 
    ind_tr_1=1;
else
    ind_tr_1=varargin{1};
end

n_trial_beh=length(behData.Trial_index);
if n_trial_beh~=length(ind_1stFrame)
    warning('dff data and behavior data trial number not match');
end
ind_1stFrame=reshape(ind_1stFrame,1,[]);%ensure that array are of same size
n_trial_output=min(n_trial_beh,length(ind_1stFrame));%if trial number not consistent, choose the smaller one
ind_trial=ind_tr_1:ind_tr_1+n_trial_output-1;%the selected trials for further analysis
ind_1stFrame=ind_1stFrame(ind_trial);
if isfield(behData,'Opto_trial_index') % for opto sessions, get the opto onset and offset
    optoOnset=behData.Opto_Onset_Time(ind_trial);
    optoOffset=behData.Opto_Off_Time(ind_trial);
    optoOnset(behData.Opto_trial_index==0)=nan;
    optoOffset(behData.Opto_trial_index==0)=nan;
    behEventFrameIndex.optoOnset=double(round(optoOnset/frameTime))+ind_1stFrame;
    behEventFrameIndex.optoOffset=double(round(optoOffset/frameTime))+ind_1stFrame;
end
stimOnset=behData.Stim_onset_time(ind_tr_1:ind_tr_1+n_trial_output-1);
stimOffset=behData.Stim_offset_time(ind_tr_1:ind_tr_1+n_trial_output-1);
goTime=behData.Go_time(ind_tr_1:ind_tr_1+n_trial_output-1);
delayOffset=behData.Time_delayOffset(ind_tr_1:ind_tr_1+n_trial_output-1);
ansTime=double(behData.Answer_time(ind_tr_1:ind_tr_1+n_trial_output-1)); %zero when no answer, so put to nan
rewTime=double(behData.Reward_time(ind_tr_1:ind_tr_1+n_trial_output-1)); %zero when no reward, so put to nan
ansTime(ansTime==0)=nan;
rewTime(rewTime==0)=nan;
leftLick=behData.Left_lick_time(ind_tr_1:ind_tr_1+n_trial_output-1);
rightLick=behData.Right_lick_time(ind_tr_1:ind_tr_1+n_trial_output-1);
lickFirst_left=zeros(1,length(leftLick));
lickFirst_right=zeros(1,length(rightLick));
lickLast_left=zeros(1,length(leftLick));
lickLast_right=zeros(1,length(rightLick));
lickFirst=lickFirst_left;
lickLast=lickLast_left;
for i=1:length(leftLick)
    leftLick{i}=double(leftLick{i});
    rightLick{i}=double(rightLick{i});
    leftLick{i}(leftLick{i}<stimOnset(i))=[];
    rightLick{i}(rightLick{i}<stimOnset(i))=[];
    if isempty(leftLick{i})
        lickFirst_left(i)=nan;
        lickLast_left(i)=nan;
    else
        lickFirst_left(i)=leftLick{i}(1);
        lickLast_left(i)=leftLick{i}(end);
    end
    if isempty(rightLick{i})
        lickFirst_right(i)=nan;
        lickLast_right(i)=nan;
    else
        lickFirst_right(i)=rightLick{i}(1);
        lickLast_right(i)=rightLick{i}(end);
    end
    if isempty(rightLick{i}) && isempty(leftLick{i})
        lickFirst(i)=nan;
        lickLast(i)=nan;
    else
        lickFirst(i)=nanmin(lickFirst_left(i),lickFirst_right(i));
        lickLast(i)=nanmax(lickLast_left(i),lickLast_right(i));
    end
end
behEventFrameIndex.stimOnset=double(round(stimOnset/frameTime))+ind_1stFrame;%Integers can only be combined with integers of the same class, or scalar doubles.
behEventFrameIndex.stimOffset=double(round(stimOffset/frameTime))+ind_1stFrame;
behEventFrameIndex.ansTime=double(round(ansTime/frameTime))+ind_1stFrame;
behEventFrameIndex.rewTime=double(round(rewTime/frameTime))+ind_1stFrame;
behEventFrameIndex.lickFirst_left=double(round(lickFirst_left/frameTime))+ind_1stFrame;
behEventFrameIndex.lickFirst_right=double(round(lickFirst_right/frameTime))+ind_1stFrame;
behEventFrameIndex.lickFirst=double(round(lickFirst/frameTime))+ind_1stFrame;
behEventFrameIndex.lickLast=double(round(lickLast/frameTime))+ind_1stFrame;
behEventFrameIndex.go=double(round(goTime/frameTime))+ind_1stFrame;
behEventFrameIndex.delayOffset=double(round(delayOffset/frameTime))+ind_1stFrame;
behEventFrameIndex.start=ind_1stFrame;
leftLick=cellfun(@(x) double(round(x/frameTime)),behData.Left_lick_time,'UniformOutput',false);
rightLick=cellfun(@(x) double(round(x/frameTime)),behData.Right_lick_time,'UniformOutput',false);
for i=1:length(ind_1stFrame)
    lickingFrameIndex.leftLick{i}=leftLick{i}+ind_1stFrame(i);
    lickingFrameIndex.rightLick{i}=rightLick{i}+ind_1stFrame(i);
end
% save(behEventFrameIndex,behEventFrameIndex);
end

