%first extract data to a formated matrix (ref formDataMatrix) and perform analyses on that.
close all;
clear;
dbstop if error;
[num,txt,raw] =xlsread('C:\Users\PYX\Documents\DataSummary\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
savepath='E:\2P\summary\population_activities';

T=cell2table(raw(2:end,1:15));
T.Properties.VariableNames=strrep(raw(1,1:15),' ','_');%table variable name can't have ' ',so replace them
ind_session1=logical(strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control'));
% ind_session1=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*(~contains(T.behavior_performance,'probe'));%not probe session
% ind_session2=logical(contains(T.cell_type,'syn').*contains(T.ROI_type,'soma'));

%extract data to a matrix
celltypePool={'M2','vglut2','vgat'};
%decide some global variable
behEventAlignPool={'stim onset','go cue','first lick'};%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
frameNumTime_cell={[0.5,1],[1,0.5],[0.5,1]};%corresponding time window around behavior event aligned
behEventMaskPool={'go cue','first lick',[]};%mask activities after events to rule out confound
%{
for i_celltype=1:length(celltypePool)
    ind_session2=logical(strcmp(T.cell_type,celltypePool{i_celltype})+strcmp(T.cell_type,[celltypePool{i_celltype},'-flpo']));
    ind_session= ind_session1 & ind_session2;
    ind_trial_chosen4ScatterHist=logical(contains(T.ROI_type,'soma').*(~contains(T.field,'soma')));
    ind_session= ind_session & ind_trial_chosen4ScatterHist;
    ind_session=find(ind_session);
    n_session=length(ind_session);
    animal_unique=unique(T.animal(ind_session));
    n_animal=length(animal_unique);
    
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
%     frT =    33.3296;
    nROI=0;
    for i_session=1:n_session
        indrow=ind_session(i_session);
        trial2include=T.used_trial{indrow};
        trial2exclude=[];
        session=[T.session{indrow},'_',T.field{indrow}];
        file_path=T.file_path{indrow};
        objsession=Session2P(session,file_path,trial2include,trial2exclude);
        disp(['processin session-',objsession.name]);
        [format_dff, trial_label, ts] = objsession.mGetFormatDff(activity_type,behEventAlignPool,frameNumTime_cell,frT,'maskAfter',behEventMaskPool);
        save([savepath,filesep,celltypePool{i_celltype},filesep,objsession.name,'-withMask.mat'],'format_dff', 'trial_label', 'ts');
        nROI=nROI+objsession.metadata.nROI;
    end
    metadata.behEventAlign=behEventAlignPool;
    metadata.frameNumTime=frameNumTime_cell;
    metadata.frameTime=frT;
    metadata.sampleRate=1000/frT;
    metadata.cellNumber=nROI;
    save([savepath,filesep,celltypePool{i_celltype},filesep,'withMask-metadata.mat'],'metadata');
end
%}
%calcualte mean trace for each trial types and group cells together as
%pseudo-population
ppdff=[];%pseudo-population
for i_celltype=1:length(celltypePool)
    load_path=[savepath,filesep,celltypePool{i_celltype}];
    dirmat=strcat(load_path,filesep,'*.mat');
    dirs=dir(dirmat);
    dircell=struct2cell(dirs);
    filenames=dircell(1,:);
    files_dff=cellfun(@(x) ~contains(x,'metadata'), filenames);
    ind_files_dff=find(files_dff);
    n_session=sum(files_dff);
    for i=1:n_session
        load([load_path,filesep,filenames{ind_files_dff(i)}]);
        disp(['loading',filenames{files_dff(i)},' nROI=',num2str(size(format_dff,1))])
        pdff=zeros(size(format_dff,1),size(format_dff,2),2);%population dff mean trace
        pdff(:,:,1)=nanmean(format_dff(:,:,trial_label==-1),3);
        pdff(:,:,2)=nanmean(format_dff(:,:,trial_label==1),3);
        ppdff=cat(1,ppdff,pdff);
    end
end
%PCA using data across epochs
ppdff2d=reshape(ppdff,[size(ppdff,1),size(ppdff,2)*2,1]);%used for further PCA etc.
[coeff_PCA,score,latent,tsquared,explained,mu] = pca(ppdff2d');

corr_mat=corr(ppdff2d');





