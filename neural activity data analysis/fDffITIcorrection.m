function [ dffoutput ] = fDffITIcorrection( dffraw, behEventFrameIndex,fr,varargin )
%FDFFITICORRECTION subtract ITI activiities before stim onset 1s as
%baseline to get a dff with same size of input
%Input-
%   dffraw- dff to be corrected;
%   behEventFrameIndex- frame index of beh event, used to extract baseline
%   epoch
%   fr- frame rate of sampling
%   varargin- can be 1) baselineLength- time of baseline duration, default
%   0.5
if isempty(varargin)
    baselineLength=0.5;
else
    baselineLength=varargin{1};
end
for i=1:length(behEventFrameIndex.stimOnset)
    ind_trialstart=behEventFrameIndex.start(i);
    if i<length(behEventFrameIndex.stimOnset)
        ind_trialend=behEventFrameIndex.start(i+1);
    else
        ind_trialend=size(dffraw,2);
    end
    baseline_activity=dffraw(:,behEventFrameIndex.stimOnset(i)-baselineLength*fr:behEventFrameIndex.stimOnset(i));
    baseline=nanmean(baseline_activity,2);
    dffraw(:,ind_trialstart:ind_trialend)=dffraw(:,ind_trialstart:ind_trialend)-baseline;
end
dffoutput=dffraw;
end

