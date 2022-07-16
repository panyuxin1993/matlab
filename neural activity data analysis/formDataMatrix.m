close all;
clear;
[num,txt,raw] =xlsread('C:\Users\PYX\Documents\DataSummary\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
savepath='E:\2P\summary\format_data_for_population';

T=cell2table(raw(2:end,1:15));
T.Properties.VariableNames=strrep(raw(1,1:15),' ','_');%table variable name can't have ' ',so replace them
ind_session1=logical(strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control'));
% ind_session1=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*(~contains(T.behavior_performance,'probe'));%not probe session
% ind_session2=logical(contains(T.cell_type,'syn').*contains(T.ROI_type,'soma'));

celltypePool={'M2','vglut2','vgat'};
for i_celltype=1:length(celltypePool)
    ind_session2=logical(strcmp(T.cell_type,celltypePool{i_celltype})+strcmp(T.cell_type,[celltypePool{i_celltype},'-flpo']));
    ind_session= ind_session1 & ind_session2;
    ind_trial_chosen4ScatterHist=logical(contains(T.ROI_type,'soma').*(~contains(T.field,'soma')));
    ind_session= ind_session & ind_trial_chosen4ScatterHist;
    ind_session=find(ind_session);
    n_session=length(ind_session);
    animal_unique=unique(T.animal(ind_session));
    n_animal=length(animal_unique);
    
    %decide some global variable
    behEventAlignPool={'stim onset','go cue','first lick'};%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
    frameNumTime_cell={[0.5,1],[1,1],[0.5,1]};
    
    activity_type='dff';
    frT=100;
    for i_session=1:n_session%14%
        indrow=ind_session(i_session);
        trial2include=T.used_trial{indrow};
        trial2exclude=[];
        session=[T.session{indrow},'_',T.field{indrow}];
        file_path=T.file_path{indrow};
        objsession=Session2P(session,file_path,trial2include,trial2exclude);
        frT=min(frT,objsession.metadata.frT);
        disp(['reading session-',objsession.name]);
    end
    %frT =    33.3296;
    nROI=0;
    for i_session=1:n_session
        indrow=ind_session(i_session);
        trial2include=T.used_trial{indrow};
        trial2exclude=[];
        session=[T.session{indrow},'_',T.field{indrow}];
        file_path=T.file_path{indrow};
        objsession=Session2P(session,file_path,trial2include,trial2exclude);
        disp(['processin session-',objsession.name]);
        [format_dff, trial_label, ts] = objsession.mGetFormatDff(activity_type,behEventAlignPool,frameNumTime_cell,frT);
        save([savepath,filesep,celltypePool{i_celltype},filesep,objsession.name,'.mat'],'format_dff', 'trial_label', 'ts');
        nROI=nROI+objsession.metadata.nROI;
    end
    metadata.behEventAlign=behEventAlignPool;
    metadata.frameNumTime=frameNumTime_cell;
    metadata.frameTime=frT;
    metadata.sampleRate=1000/frT;
    metadata.cellNumber=nROI;
    save([savepath,filesep,celltypePool{i_celltype},filesep,'metadata.mat'],'metadata');
end