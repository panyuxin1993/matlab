function [label] = fTrialType2Label(trialType,id)
%FTRIALTYPE2LABEL in order to perform AUC analysis, get a label vector from
%trialType parameter 
%Input-
%   trialType=3d/4d matrix made of 0/1, indicating trial type,usually the
%   2nd(method of grouping) or 4th(if exist, whether opto or not) dimension
%   is used, i.e. transformed to vector along that dimension
%   id=2(default)|4, which dimension to transform, 
%   {labelstr- 1-by-n string vector, n match the size of dimention 'id',default is using 1-n number}
%Output-
%   label- m-by-1 string vector, m matches size of 3rd dimention of
%   trialType
n_d=length(size(trialType));
if n_d<id
    warning('exceed dimension of trialType');
end
ntrial=size(trialType,3);
label=zeros(ntrial,1);
if n_d==4 && id==4%label by the 4-d
    temp=sum(trialType,[1,2]);
elseif n_d==4 && id==2 
    temp=sum(trialType,[1,4]);
elseif n_d==3 && id==2 
    temp=sum(trialType,1);
elseif n_d==3 && id==1
    temp=sum(trialType,2);
else
    warning('unreasonable parameter for which dimension to be transform');
end
label2d=squeeze(temp);%now label2d is m-by-n matrix, n is trial number, m is number of category
for i=1:ntrial
    for ind_category=1:size(label2d,1)
        if label2d(ind_category,i)>0
            label(i)=ind_category;
        end
    end
end
label(label==0)=nan; %if not signed value, maybe some unwanted category, so set to nan
