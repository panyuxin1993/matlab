close all;
[num,txt,raw] =xlsread('C:\Users\PYX\Documents\DataSummary\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
savepath='E:\2P\summary';
clear TAUC_combine Tmean_combine;
trialTypeStr='cor';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
AUCtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
AUCCorrectedMethod='SensoryChoiceOrthogonalSubtraction';%'SensoryChoiceOrthogonalSubtraction';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
manipulation='control';

trialTypeStr='cor and err';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
AUCtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
load([savepath,filesep,'trialType',trialTypeStr,'-TepochMeanActivity.mat']);

%% case, see different cell type
celltype={'M2','syn','vglut2','vgat'};
manipulation='control';
%T_AUC=TAUC_combine(strcmp(TAUC_combine.celltype,'syn')|contains(TAUC_combine.celltype,'+other'),:);%syn, ..+other
%T_AUC=TAUC_combine(strcmp(TAUC_combine.celltype,celltype),:);
% ind_trial_chosen=strcmp(Tmean_combine.animal,'pyx358') & (strcmp(Tmean_combine.date,'2021/6/10') | strcmp(Tmean_combine.date,'2021/5/30'));

ind_trial_chosen=contains(Tmean_combine.ROItype,'soma');
% ind_trial_chosen=Tmean_combine.usd_as_data
T_mean=Tmean_combine(ind_trial_chosen,:);

psig=0.05;
ind_sound=logical(T_mean.sound_pselectivity<psig);
ind_delay=T_mean.delay_pselectivity<psig;
ind_response=T_mean.response_pselectivity<psig;
ind_lick=T_mean.lick_pselectivity<psig;
prob_sig=zeros(length(celltype),4);%2d, behavioral epochs
for i_celltype=1:length(celltype)
    ind_temp=cellfun(@(x) strcmp(x,celltype{i_celltype}),T_mean.celltype);
    if sum(ind_temp)==0
        warning('no data included');
    end
    prob_sig(i_celltype,1)=sum(ind_sound.*ind_temp)/sum(ind_temp);
    prob_sig(i_celltype,2)=sum(ind_delay.*ind_temp)/sum(ind_temp);
    prob_sig(i_celltype,3)=sum(ind_response.*ind_temp)/sum(ind_temp);
    prob_sig(i_celltype,4)=sum(ind_lick.*ind_temp)/sum(ind_temp);
end
bar(prob_sig');
set(gca,'XTick',1:4,'XTickLabel',{'sound','delay','response','lick'},'FontSize',10);
ylabel('Proportion of significant selectivity');
legend(celltype,'Location','northwest');
datastr=strrep(T_mean.date{1},'/','');
%}

%% case, see different ROI type
%{
% ind_trial_chosen=strcmp(Tmean_combine.animal,'pyx358') & strcmp(Tmean_combine.date,'2021/6/10') ;
ROItype_all={'soma','dendrite','spine'};
psig=0.05;
ind_sound=logical(T_mean.sound_pselectivity<psig);
ind_delay=T_mean.delay_pselectivity<psig;
ind_response=T_mean.response_pselectivity<psig;
ind_lick=T_mean.lick_pselectivity<psig;
prob_sig=zeros(length(ROItype_all),4);%2d, behavioral epochs
for i_ROItype=1:length(ROItype_all)
    ind_temp=cellfun(@(x) strcmp(x,ROItype_all{i_ROItype}),T_mean.ROItype);
    prob_sig(i_ROItype,1)=sum(ind_sound.*ind_temp)/sum(ind_temp);
    prob_sig(i_ROItype,2)=sum(ind_delay.*ind_temp)/sum(ind_temp);
    prob_sig(i_ROItype,3)=sum(ind_response.*ind_temp)/sum(ind_temp);
    prob_sig(i_ROItype,4)=sum(ind_lick.*ind_temp)/sum(ind_temp);
end
bar(prob_sig');
set(gca,'XTick',1:4,'XTickLabel',{'sound','delay','response','lick'},'FontSize',10);
ylabel('Proportion of significant selectivity');
legend(ROItype_all,'Location','northwest');
datastr=strrep(T_mean.date{1},'/','');
title([T_mean.animal{1},'\_',datastr]);
set(gcf,'PaperPosition',[0,0,4,4]);
saveas(gcf,['E:\2P\example\dentrites_vs_spines\',T_mean.animal{1},'_',datastr,'bar_prob_sig.pdf'],'pdf');
%}