function [Tout]=fOrthogonalSubtraction(Tin,ind_trial,label_SVM,label_SVM_orthogonal)
%FORTHOGONALSUBTRACTION subtract unrelated activities, e.g. when
%calculating choice AUC/SVM, subtract sensory activities
%Input- 
%   Tin- table or array of data to be manipulated
%   ind_trial- choosen data to be used (may exclude trials, such as miss,
%       violation, movement etc.
%   label_SVM- label to be used to calculate AUC/SVM, 
%   label_SVM_orthogonal- label to avoid confound when calculating AUC/SVM

if istable(Tin)
    Tout=table2array(Tin(ind_trial,:));
else
    Tout=Tin(ind_trial,:);
end
n_condition=length(unique(label_SVM));
n_condition_orth=length(unique(label_SVM_orthogonal));
for i=1:n_condition_orth
    mean_orth=zeros(n_condition,size(Tout,2));
    for j=1:n_condition
        indTrial=logical((label_SVM==j).*(label_SVM_orthogonal==i));
        mean_orth(j,:)=nanmean(Tout(indTrial,:));   
    end
    indTrial_orth=(label_SVM_orthogonal==i);
    Tout(indTrial_orth,:)=Tout(indTrial_orth,:)-nanmean(mean_orth);
end
if istable(Tin)
    Tout=array2table(Tout,'VariableNames',Tin.Properties.VariableNames);
end
end

% orignal code to perform when calculating AUC
%{
function [Tout]=fOrthogonalSubtraction(Tin,ind_trial,label_AUC,label_AUC_orthogonal)
if istable(Tin)
    Tout=table2array(Tin(ind_trial,:));
else
    Tout=Tin(ind_trial,:);
end
n_condition=length(unique(label_AUC));
n_condition_orth=length(unique(label_AUC_orthogonal));
for i=1:n_condition_orth
    mean_orth=zeros(n_condition,size(Tout,2));
    for j=1:n_condition
        indTrial=logical((label_AUC==j).*(label_AUC_orthogonal==i));
        mean_orth(j,:)=nanmean(Tout(indTrial,:));   
    end
    indTrial_orth=(label_AUC_orthogonal==i);
    Tout(indTrial_orth,:)=Tout(indTrial_orth,:)-nanmean(mean_orth);
end
if istable(Tin)
    Tout=array2table(Tout,'VariableNames',Tin.Properties.VariableNames);
end
end
%} 