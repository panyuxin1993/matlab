close all;
behEventAlign='delay onset';
masklick='no';
behEventSort='go cue';
bodypartsPool={'Tongue'};
coordinatesPool={'x'};
titlestrPool={'Tongue (x)'};
yrange={[-40,40]};
selectivity='choice';
treshold4likelihood=0.1;
pSigTtest=0.01;
nonOverlappingBin=1;
baseline='';%'LickPort';

[num,txt,raw] =xlsread('C:\Users\PYX\Documents\DataSummary\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
T=cell2table(raw(2:end,1:15));
T.Properties.VariableNames=strrep(raw(1,1:15),' ','_');%table variable name can't have ' ',so replace them
ind_session=logical(strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*(~contains(T.behavior_performance,'probe')).*contains(T.ROI_type,'soma'));%not probe session
Tchoose=T(ind_session,:);
% file_path='H:\2P\pyx349_20210418\im_data_reg\result_save';
% session='pyx349_20210418';

behEventAlignPool={'delay onset','go cue','first lick'};
masklickPool={'yes','no','no'};
i_selectivity=3;%3-choice,4-sensory
behEventSortPool={'go cue','first lick','go cue'};
trial2include='all';
trial2exclude=[];
activity_form='dff';%{'dff','spkr'};
chosen_result={'Correct','Error'};%{'Correct','Error','Miss','Violation'};
chosen_stim={'ipsi','contra'};%{'ipsi','contra'};
n_trial_show=50;
ind_ROI=4;

%plot
%method 1, define the index
%{
i=1;
file_path=Tchoose.file_path{i};
session=[Tchoose.session{i},'_',Tchoose.field{i}];
%}

%method 2, directly give path and session name
%
file_path='E:\2P\CD058_20180201\im_data_reg\result_save';
session='CD058_20180201';
%}


savepath=[file_path,filesep,session];
objsession=Session2P(session,file_path,trial2include,trial2exclude);

% close all;
[fig_A,fig_B]=objsession.mInspectActivity(ind_ROI,activity_form,behEventAlign,masklick,behEventSort,chosen_result,chosen_stim,n_trial_show);
