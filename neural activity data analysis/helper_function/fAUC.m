function [AUC,pvalue] = fAUC(label,activity,poslabel,n_shuffle)
%FAUC calculate AUC and corresponding p value of one ROI, 
%   input: label-true label,n-by-1 vector
%   activity-neural activities for ROC analysis, n-by-1 vector
%   poslabel-label of interest
%   n_shuffle-times of shuffle used to calculate p value of AUC
if nargin<4
    n_shuffle=1000; %default shuffle times is 100
end

unilabel=unique(label);
if length(unilabel)~=2
    warning('not 2 category for AUC calculation');
    AUC=nan;
    pvalue=nan;
    return;
end
ind_nan=isnan(activity);
label1=(unilabel(1)==label);
if sum((~ind_nan).*label1)<3 ||sum((~ind_nan).*(~label1))<3
    AUC=nan;
    pvalue=nan;
    return;
else
    [X,Y,T,AUC]=perfcurve(label,activity,poslabel);
end
% figure;
% subplot(1,2,1);
% plot(X,Y);hold on;
% plot([0,1],[0,1],'k');%diagonal
% subplot(1,2,2);
% histogram(activity(label1),'BinWidth',0.1);hold on;
% histogram(activity(~label1),'BinWidth',0.1);
% pause(1);close;
%method 1-using for loop
AUCs=zeros(n_shuffle,1);
parfor i=1:n_shuffle%during shuffle, maybe nan will happen for one category, so before computing AUC, check that
    ind=randperm(length(label));
    label_i=label(ind);
    ind_category1=(label_i==unilabel(1));
    if sum((~ind_nan).*ind_category1)<3 ||sum((~ind_nan).*(~ind_category1))<3 %if either category's non-nan is less than 3, then can not compute AUC
        AUCs(i)=nan;
    else
        [X,Y,T,AUCs(i)]=perfcurve(label_i,activity,poslabel);
%         plot(X,Y);%pause(1);close;
    end     
end
%method 2-using cellfun
% [labelcell{1:n_shuffle}]=deal(label);
% [activitycell{1:n_shuffle}]=deal(activity);
% [poslabelcell{1:n_shuffle}]=deal(poslabel);
% [unilabelcell{1:n_shuffle}]=deal(unilabel);
% [ind_nancell{1:n_shuffle}]=deal(ind_nan);
% AUCs=cellfun(@fAUCshuffle,labelcell,activitycell,poslabelcell,unilabelcell,ind_nancell);

%calculate p-value
nAUC=sum(~isnan(AUCs));
if AUC==0.5
    pvalue=0.5;
else
    pvalue=nansum(AUC<AUCs)/nAUC;%one-tail
    %pvalue=2*(0.5-abs(0.5-nansum(AUC<AUCs)/nAUC));%two-tail
end
%visualize the AUC in distribution for confirmation of p-value
%{
histogram(AUCs);
hold on;
plot(AUC,0,'ko');
close;
%}
end

function [AUC]=fAUCshuffle(label,activity,poslabel,unilabel,ind_nan)
ind=randperm(length(label));
label_i=label(ind);
ind_category1=(label_i==unilabel(1));
if sum((~ind_nan).*ind_category1)<3 ||sum((~ind_nan).*(~ind_category1))<3 %if either category's non-nan is less than 3, then can not compute AUC
    AUC=nan;
else
    [~,~,~,AUC]=perfcurve(label_i,activity,poslabel);
end
end