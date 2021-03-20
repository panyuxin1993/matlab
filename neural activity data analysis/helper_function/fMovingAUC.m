function [AUCout,Pout,activityout] = fMovingAUC(label,activity,poslabel,n_shuffle,binsize,binstep)
%FMOVINGAUC using a temporal series to calculate moving AUC
%Input-
%   label-true label,n-by-1 vector
%   activity-neural activities for ROC analysis, n-by-m matrix, where n
%   trials each has m frames
%   poslabel-label of interest
%   n_shuffle-times of shuffle used to calculate p value of AUC
%   binsize- size of bin to calculate mean activities used for AUC
%   calculation
%   binstep- moving step of bin to calculate AUC
%Output-
%   AUCout- 1-by-k double, moving AUC
%   Pout- 1-by-k double, moving P of AUC
%   activityout- 2-by-1 cell, each n-by-k double, smoothed activity, note only support two
%   category
% nframe=floor((size(activity,2)-binsize)/binstep)+1;
nframe=size(activity,2);
AUCout=zeros(1,nframe);
Pout=ones(1,nframe);
activityout=cell(2,1);
tempactivityout=nan(size(activity));
for i=1:nframe
    ind1=max(1,(i-1)*binstep-floor(binsize/2)+1);
    ind2=min(size(activity,2),(i-1)*binstep+floor(binsize/2)+1);
    binActivity=activity(:,ind1:ind2);
    binMeanActivity=nanmean(binActivity,2);%n-by-1 vector
    tempactivityout(:,i)=binMeanActivity;
    ind_nan=isnan(binMeanActivity);
    ind_category1=(label==1);%ind_category2=(~ind_category1);
    if sum((~ind_nan).*ind_category1)<3 ||sum((~ind_nan).*(~ind_category1))<3 %if either category's non-nan is less than 3, then can not compute AUC
        AUCout(i)=nan;
        Pout(i)=nan;
    else
        [AUCout(i),Pout(i)] = fAUC(label,binMeanActivity,poslabel,n_shuffle);
    end
end
activityout{1}=tempactivityout(label==1,:);
activityout{2}=tempactivityout(label==2,:);
if binsize~=1
    warning('bin size is not 1, so check if mean trace also be smoothed');
end
end

