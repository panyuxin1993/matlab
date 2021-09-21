close all;
[num,txt,raw] =xlsread('D:\xulab\project\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
savepath='H:\2P\summary\summary_DLCfiltered';
trialTypeStr='cor and err';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
manipulation='control';
celltype_Str='syn';
T=cell2table(raw(2:end,1:17));
T.Properties.VariableNames=strrep(raw(1,1:17),' ','_');%table variable name can't have ' ',so replace them
ind_session=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control')...
    .*strcmp(T.behavior_video,'DLC tracked');%.*strcmp(T.cell_type,celltype_Str);
ind_session=find(ind_session);
n_session=length(ind_session);
animal_unique=unique(T.animal(ind_session));
n_animal=length(animal_unique);
epoch='delay';
str_nFrames='500ms';
bodyparts='Tongue';
coordinates='x';
behEventAlign='delay onset';
bodypartsPool={'Tongue'};
coordinatesPool={'x'};
behEventAlignPool={'delay onset'};
titlestrPool={'Tongue (x)'};
yrange={[-40,40]};
selectivity='choice';
treshold4likelihood=0.1;
pSigTtest=0.01;
nonOverlappingBin=1;
baseline='';%'LickPort';
SVMtype=selectivity;
nRepeat=100;
pTraining=0.9;
for i_session=1:n_session%%%%%%%%%%%
    indrow=ind_session(i_session);
    objsession=Session2P(T.session{indrow},T.file_path{indrow},T.used_trial{indrow},nan);
    savename=[T.session{indrow},'-',epoch];
    [ind_NSDelayMovement,DelayMovement_criteria]=objsession.mGetIndNSDdelayMovement(baseline);
%     [nfig,struct_rho,struct_p] = objsession.mPlotDLCvsDff(epoch,str_nFrames,trialTypeStr,...
%         bodyparts,coordinates,treshold4likelihood,savepath,savename);
    figDLC_PSTH=objsession.mPlotDLC_PSTH(bodypartsPool,coordinatesPool,ind_NSDelayMovement,...
        selectivity,pSigTtest,treshold4likelihood,behEventAlignPool,...
        yrange,titlestrPool,nonOverlappingBin,savepath,DelayMovement_criteria);
%     TSVM_currentSession=objsession.mSVMscore(T.cell_type{indrow},trialTypeStr,ind_NSDelayMovement,SVMtype,nRepeat, pTraining);
%     if exist('TSVM_combine','var')
%         TSVM_combine=vertcat(TSVM_combine,TSVM_currentSession);
%     else
%         TSVM_combine=TSVM_currentSession;
%     end
end
save([savepath,filesep,'trialType',trialTypeStr,'-',SVMtype,'pTraining',num2str(pTraining),'-TepochSVM.mat'],'TSVM_combine');
