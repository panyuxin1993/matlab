function [pTtest] = fMovingTtest(label,activity,binsize,binstep,varargin)
%FMOVINGTTEST using a temporal series to calculate moving t-test
%Input-
%   label-true label,n-by-1 vector
%   activity-neural activities for ROC analysis, n-by-m matrix, where n
%   trials each has m frames
%   binsize- size of bin to calculate mean activities used for t-test
%   binstep- moving step of bin to calculate t-test
%Output-
%   pTtest- 1-by-k double, moving P of ttest
if isempty(varargin)
    testName='ttest2';
else
    testName=varargin{1};
end
% nframe=floor((size(activity,2)-binsize)/binstep)+1;
nframe = size(activity,2);
pTtest=ones(1,nframe);
if strcmp(testName,'ttest2')
    unilabel=unique(label);%usually should be 2
    if length(unilabel)~=2
        warning('t-test need 2 category');
    else
        activityx=activity((label==unilabel(1)),:);
        activityy=activity((label==unilabel(2)),:);
    end
elseif strcmp(testName,'ttest')%input should be change
    activityx=label;
    activityy=activity;
end
if exist('activityx','var')%in miss trials, maybe only one category, can not make ttest
    for i=1:nframe
        ind1=max(1,(i-1)*binstep-floor(binsize/2)+1);
        ind2=min(size(activity,2),(i-1)*binstep+floor(binsize/2)+1);
        binActivityx=activityx(:,ind1:ind2);
        binMeanActivityx=nanmean(binActivityx,2);%n-by-1 vector
        binActivityy=activityy(:,ind1:ind2);
        binMeanActivityy=nanmean(binActivityy,2);%n-by-1 vector
        if sum(~isnan(binMeanActivityx))<3 ||sum((~isnan(binMeanActivityy)))<3 %if either category is all nan, then can not compute AUC
            pTtest(i)=nan;
        else
            switch testName
                case 'ttest'
                    [~,pTtest(i)]=ttest(binMeanActivityx,binMeanActivityy);
                case 'ttest2'
                    [~,pTtest(i)]=ttest2(binMeanActivityx,binMeanActivityy);
            end
        end
    end
end
end

