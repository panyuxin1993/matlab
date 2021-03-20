function [Pout,activityout] = fANOVA1(activity,group)
%FANOVA1 calculate ANOVA1 p-value
%Input-
%   group-true label,n-by-1 vector
%   activity-neural activities for ROC analysis, n-by-1 matrix, where n
%   trials each 1 value
%Output-
%   Pout- 1-by-1 double

labels=unique(group);
binMeanActivity=activity;%n-by-1 vector
tempactivityout=binMeanActivity;
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
    Pout = anova1(binMeanActivity,group,'off');
else
    warning('not normal distribution, replace anova1 with kruskalwallis test');
    Pout = kruskalwallis(binMeanActivity,group,'off');
end
for i=1:length(labels)
    activityout{i}=tempactivityout(group==labels(i));
end
end

