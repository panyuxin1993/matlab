function [ meanout, ciout, prctileout,activityout ] = fMovingBootstrpMeanCIPrctile( activity,nboot,binsize,binstep ,varargin)
%FMOVINGBOOTSTRPMEANCIPCTIL using a temporal series to calculate moving
%Bootstrp mean, ci and prctile of 0 in the resampled distribution
%Input-
%   activity-neural activities for bootstrp analysis, n-by-m matrix, where 
%   n trials each has m frames
%   n_shuffle-times of shuffle used to calculate distribution of data
%   binsize- size of bin to calculate mean activities used for bootstrp
%   calculation, is in fact smooth
%   binstep- moving step of bin to calculate bootstrp stats
%Output-
%   meanout- 1-by-k double, moving mean of resampled data
%   ciout- 2-by-k double, moving 95% confidence interval
%   prctileout- 1-by-k double, moving percentile of 0 in the distribution 
%   activityout- same size of activity, useful for smoothed process. Thus
%   garantee the activity and stats are corresponding.

if ~isempty(varargin)
    baseline=varargin{1};
else
    baseline=0;
end
nframe=size(activity,2);
meanout=nan(1,nframe);
prctileout=nan(1,nframe);
ciout=nan(2,nframe);
activityout=nan(size(activity));
for i=1:nframe
    ind1=max(1,(i-1)*binstep-floor(binsize/2)+1);
    ind2=min(size(activity,2),(i-1)*binstep+floor(binsize/2)+1);
    binActivity=activity(:,ind1:ind2);
    binMeanActivity=nanmean(binActivity,2);%n-by-1 vector
    activityout(:,i)=binMeanActivity;
    [meanout(i),ciout(:,i),prctileout(i)] = fBootstrpMeanCIPrctile(binActivity,nboot,baseline);
end
end

function [meanout,ciout,prctileout] = fBootstrpMeanCIPrctile(activity,nboot,baseline)
bootmean=bootstrp(nboot,@nanmean,activity);
meanout=nanmean(nanmean(bootmean));
ciout=[prctile(bootmean,0.5);prctile(bootmean,99.5)];
prctileout=sum(bootmean<baseline)/sum(~isnan(bootmean));
end
