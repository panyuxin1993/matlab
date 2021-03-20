function [Pout,activityout] = fMovingANOVA1(activity,group,binsize,binstep)
%FMOVINGANOVA2 using a temporal series to calculate moving ANOVA1 p-value
%Input-
%   group-true label,n-by-1 vector
%   activity-neural activities for ROC analysis, n-by-m matrix, where n
%   trials each has m frames
%   binsize- size of bin to calculate mean activities used for ANOVA1
%   calculation
%   binstep- moving step of bin to calculate ANOVA1
%Output-
%   Pout- 1-by-k double, moving P of ANOVA1
nframe=size(activity,2);
Pout=ones(1,nframe);
labels=unique(group);
for i=1:nframe
    ind1=max(1,(i-1)*binstep-floor(binsize/2)+1);
    ind2=min(size(activity,2),(i-1)*binstep+floor(binsize/2)+1);
    binActivity=activity(:,ind1:ind2);
    binMeanActivity=nanmean(binActivity,2);%n-by-1 vector
    tempactivityout(:,i)=binMeanActivity;
    %test for normal distribution
    for ilabel=1:length(labels)
        flagNor=1;
        [h,p] = chi2gof(binMeanActivity(group==labels(ilabel)),'Alpha',0.05);
        if ~h
            flagNor=0;
            break;
        end
    end
    if flagNor==1
        Pout(i) = anova1(binMeanActivity,group,'off');
    else
        warning('not normal distribution, replace anova1 with kruskalwallis test');
        Pout(i) = kruskalwallis(binMeanActivity,group,'off');
    end
end
for i=1:length(labels)
    activityout{i}=tempactivityout(group==labels(i),:);
end
if binsize~=1
    warning('bin size is not 1, so check if mean trace also be smoothed');
end
end

