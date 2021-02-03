function [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,dff,frameTime,str_nFrames)
%FGETSIGBEHEPOCH use dff and behavior event frame index to extract neural
%activity for further analysis AUC or for other population analyses
%   input: behEventFrameIndex-index of behavior event across whole session
%       dff-dff across whole session
%       frameTime-?ms/frame, use to decide how many frames to include fro
%       further analysis
%   output: T_SigbyEpoch-table including ITI,sound,delay,response,lick,
%       double for correct and error trials, nan for miss and violation trials
ntrial=length(behEventFrameIndex.stimOnset);
varTypes = {'double','double','double','double','double'};
[ITI,sound,delay,response,lick]=deal(nan(ntrial,1));
nFrames1k=floor(1000/frameTime);%default, each epoch use activities from 1s
nFrames500=floor(500/frameTime);

if strcmp(str_nFrames,'500ms')
    for i=1:ntrial
        if isnan(behEventFrameIndex.ansTime(i))%miss or violation, then set the row of table as nan
            continue;
        else
            ITI(i)=nanmean(dff(max(behEventFrameIndex.stimOnset(i)-nFrames500,1):behEventFrameIndex.stimOnset(i)-1));
            baseline=ITI(i);
            ITI(i)=ITI(i)-baseline;
            sound(i)=nanmean(dff(behEventFrameIndex.stimOnset(i):behEventFrameIndex.stimOffset(i)))-baseline;
            delay(i)=nanmean(dff(max(behEventFrameIndex.go(i)-nFrames500+1,behEventFrameIndex.stimOffset(i)+1):behEventFrameIndex.go(i)-1))-baseline;
            response(i)=nanmean(dff(behEventFrameIndex.go(i):behEventFrameIndex.ansTime(i)-1))-baseline;
            lick(i)=nanmean(dff(behEventFrameIndex.ansTime(i):behEventFrameIndex.ansTime(i)+2*nFrames500))-baseline;
        end
    end
elseif strcmp(str_nFrames,'1s')
    for i=1:ntrial
        if isnan(behEventFrameIndex.ansTime(i))%miss or violation, then set the row of table as nan
            continue;
        else
            ITI(i)=nanmean(dff(max(behEventFrameIndex.stimOnset(i)-nFrames1k,1):behEventFrameIndex.stimOnset(i)-1));
            baseline=nanmean(dff(max(behEventFrameIndex.stimOnset(i)-nFrames500,1):behEventFrameIndex.stimOnset(i)-1));
            ITI(i)=ITI(i)-baseline;
            sound(i)=nanmean(dff(behEventFrameIndex.stimOnset(i):behEventFrameIndex.stimOffset(i)))-baseline;
            delay(i)=nanmean(dff(max(behEventFrameIndex.go(i)-nFrames500+1,behEventFrameIndex.stimOffset(i)+1):behEventFrameIndex.go(i)-1))-baseline;
            response(i)=nanmean(dff(behEventFrameIndex.go(i):behEventFrameIndex.ansTime(i)-1))-baseline;
            lick(i)=nanmean(dff(behEventFrameIndex.ansTime(i):behEventFrameIndex.ansTime(i)+2*nFrames1k))-baseline;
        end
    end
end

T_SigbyEpoch = table(ITI,sound,delay,response,lick,...
    'VariableNames',{'ITI','sound','delay','response','lick'});
end
