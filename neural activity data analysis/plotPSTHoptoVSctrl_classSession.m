close all;
clear;
[num,txt,raw] =xlsread('C:\Users\PYX\Documents\DataSummary\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
savepath='E:\2P\summary';
clear TAUC_combine Tmean_combine;
% trialTypeStr='cor and err';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
% AUCtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
% AUCCorrectedMethod='SensoryChoiceOrthogonalSubtraction';%'None';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
trialTypeStr='cor';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
AUCtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
AUCCorrectedMethod='None';%'None';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
manipulation='control';
%% calculate AUC
%
T=cell2table(raw(2:end,1:15));
T.Properties.VariableNames=strrep(raw(1,1:15),' ','_');%table variable name can't have ' ',so replace them
% ind_session1=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*(~contains(T.behavior_performance,'probe'));%not probe session
% ind_session1=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'opto');
ind_session1=strcmp(T.used_as_data,'?').*strcmp(T.manipulation,'opto');
%  ind_session2=logical(contains(T.cell_type,'syn').*contains(T.ROI_type,'soma'));
ind_session2=logical(strcmp(T.cell_type,'vglut2')+strcmp(T.cell_type,'vgat')+strcmp(T.cell_type,'M2')+strcmp(T.cell_type,'vglut2-flpo')+strcmp(T.cell_type,'vgat-flpo'));
% ind_session2=logical(strcmp(T.cell_type,'vglut2')+strcmp(T.cell_type,'vgat')+strcmp(T.cell_type,'vglut2-flpo')+strcmp(T.cell_type,'vgat-flpo'));
% ind_session2=logical(contains(T.field,'soma')+contains(T.field,'dendrite'));
% ind_session2=logical(strcmp(T.cell_type,'M2'));
% ind_session2=logical(~contains(T.field,'soma').*contains(T.ROI_type,'soma'));%include all soma data, but exclude zoomin fields
ind_session= ind_session1 & ind_session2;
ind_trial_chosen4ScatterHist=logical(contains(T.ROI_type,'soma').*(~contains(T.field,'soma')));
ind_session= ind_session & ind_trial_chosen4ScatterHist;
ind_session=find(ind_session);
n_session=length(ind_session);
animal_unique=unique(T.animal(ind_session));
n_animal=length(animal_unique);

%decide some global variable
behEventAlignPool={'delay onset'};%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
behEventSortPool={'go cue'};% string can be in{'first lick','reward','go cue'};
masklickPool={'no'};
frameNumTime=[1,3.5];%from 5s before align point to 5s after align point
% yrange=[ -0.5 , 1 ];
i_selectivity=3;%3-choice,4-sensory*********variable**************
activity_type='dff';
ind_ROI=[];%plot all ROIs
for i_session=n_session-1
    indrow=ind_session(i_session);
    trial2include=T.used_trial{indrow};
    trial2exclude=[];
    session=[T.session{indrow},'_',T.field{indrow}];
    file_path=T.file_path{indrow};
    savepath=[file_path,filesep,session];
    objsession=SessionOpto2P(session,file_path,trial2include,trial2exclude);
    [figRasterMean,figBehavior] = objsession.mPlotRasterPSTHoptoVSctrl(ind_ROI,activity_type,behEventAlignPool, masklickPool, behEventSortPool,i_selectivity);
    close all;
end