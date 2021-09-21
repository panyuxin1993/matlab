behEventAlign='delay onset';
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
ind_ROI=[];
behEventAlignPool={'delay onset','go cue','first lick'};
masklickPool={'yes','no','no'};
i_selectivity=3;%3-choice,4-sensory
behEventSortPool={'go cue','first lick','go cue'};
trial2include='all';
trial2exclude=nan;
activityType='dff';%{'dff','spkr'};
parfor i=1:size(Tchoose,1)
    file_path=Tchoose.file_path{i};
    session=[Tchoose.session{i},'_',Tchoose.field{i}];
    savepath=[file_path,filesep,session];
    objsession=Session2P(session,file_path,trial2include,trial2exclude);
    [ind_NSDelayMovement,DelayMovement_criteria]=objsession.mGetIndNSDdelayMovement(baseline);
    
    fig_rasterPSTH=objsession.mPlotActivityRasterPSTH(activityType,ind_ROI,...
        behEventAlignPool,masklickPool,i_selectivity,behEventSortPool);
        
end