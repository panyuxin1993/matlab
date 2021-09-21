function [TAUC] = fGetEpochAUCtableASession(filepath,animal,date,field,celltype,trialTypeStr,AUCtype)
%FGETEPOCHAUCTABLEASESSION calculate epoch AUC from mulitiple sessions,
%with variables- ITI,sound,delay,response,lick, etc. only for control
%session, not opto session
%Input-filepath, path of one imaging session where beh, dff exist
%Output- TAUC with varibles{animal, date, field, ROI, ITI, sound,
%delay, response, lick, pITI, psound, pdelay, presponse, plick, celltype};
cd(filepath);
dirmat=strcat(filepath,filesep,'*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
file_imaging=cellfun(@(x) contains(x,'Ca'), filenames);
load(filenames{file_imaging});%load imaging data
load('dff.mat');%load dff
file_beh=cellfun(@(x) contains(x,'Virables'), filenames);
load(filenames{file_beh});

datestr=strrep(date,'/','-');
if strcmp(trialTypeStr,'cor')
    combineCorErr='divideCorErr';%and only use correct trial
elseif strcmp(trialTypeStr,'cor and err')
    combineCorErr='combineCorErr';%{'combineCorErr','divideCorErr'}
end
selectivitystr={'stimuli','sensory difficulty','sensory','choice'};%sensory means grouping difficulties;
trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
i_selectivity=cellfun(@(x) strcmp(x,AUCtype),selectivitystr);%*********variable**************
i_selectivity=find(i_selectivity);
str_nFrames='1s';%'500ms';%'1s'
fileNameT=[filepath,filesep,animal,'-',datestr,'trialType',trialTypeStr,'-',AUCtype,'timeBin',str_nFrames,'-EpochAUC.mat'];
if exist(fileNameT,'file')
    load(fileNameT);
    n_ROI=size(TAUC,1);
    disp(['Table exist, use ',animal,date,';nROI=',num2str(n_ROI)]);
    
else
    nROI=size(SavedCaTrials.f_raw{1},1);
    [delayMovingAUC,pdelayMovingAUC]=deal(cell(nROI,1));
    [varanimal{1:nROI}]=deal(animal);
    varanimal=reshape(varanimal,[],1);
    [vardate{1:nROI}]=deal(date);
    vardate=reshape(vardate,[],1);
    [varfield{1:nROI}]=deal(field);
    varfield=reshape(varfield,[],1);
    [varcelltype{1:nROI}]=deal(celltype);
    varcelltype=reshape(varcelltype,[],1);
    ind_tr_1=1;
    ntr=length(SavedCaTrials.f_raw);
    frT = SavedCaTrials.FrameTime;
    % align to behavior event
    nFrameEachTrial=cellfun(@(x) size(x,2),SavedCaTrials.f_raw);
    ind_1stFrame=zeros(1,length(nFrameEachTrial));
    ind_1stFrame(1)=1;
    ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
    ind_1stFrame=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials
    frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
    frameNum=double(round(frameNumTime*1000/frT));
    [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( Data_extract, ind_1stFrame, SavedCaTrials.FrameTime );%get behavior event time
    [trialType,~,~] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix','left',combineCorErr);
    
    if strcmp(AUCtype,'choice')|| strcmp(AUCtype,'sensory')
        [ITI, sound,delay, response, lick, pITI, psound, pdelay, presponse, plick]= deal(nan(nROI,1));
        label_choice = fTrialType2Label(trialType,2);
        if contains(trialTypeStr,'cor') % 'cor'| 'cor and err'
            ind_trial=logical(reshape(sum(trialType(1,:,:),2),[],1));%only correct trials
        end
        for roiNo = 1:size(SavedCaTrials.f_raw{1},1) %SavedCaTrials.nROIs may be not true
            disp([animal,date,'rioNo',num2str(roiNo)]);
            [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,dff(roiNo,:),frT,str_nFrames);
            poslabel=2;
            nshuffle=1000;
            [ITI(roiNo),pITI(roiNo)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.ITI(ind_trial),poslabel,nshuffle);
            [sound(roiNo),psound(roiNo)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.sound(ind_trial),poslabel,nshuffle);
            [delay(roiNo),pdelay(roiNo)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.delay(ind_trial),poslabel,nshuffle);
            [response(roiNo),presponse(roiNo)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.response(ind_trial),poslabel,nshuffle);
            [lick(roiNo),plick(roiNo)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.lick(ind_trial),poslabel,nshuffle);
            %calculate moving AUC
            [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff(roiNo,:), behEventFrameIndex,  frameNum );
            nshuffle=1000;
            binsize=1;
            binstep=1;
            label = label_choice;
            for nResult=1:size(trialType,1)-2
                indTrial=trialType(nResult,:,:);
                indTrial=sum(squeeze(indTrial),1);
                indTrial=logical(squeeze(indTrial));
                if nResult==1
                    if strcmp(trialTypeStr,'cor and err')
                        [delayMovingAUC{roiNo}.do,pdelayMovingAUC{roiNo}.do] = fMovingAUC(label(indTrial),dff_aligned(indTrial,:),2,nshuffle,binsize,binstep);
                    else
                        [delayMovingAUC{roiNo}.cor,pdelayMovingAUC{roiNo}.cor] = fMovingAUC(label(indTrial),dff_aligned(indTrial,:),2,nshuffle,binsize,binstep);
                    end
                elseif nResult==2
                    [delayMovingAUC{roiNo}.err,pdelayMovingAUC{roiNo}.err] = fMovingAUC(label(indTrial),dff_aligned(indTrial,:),2,nshuffle,binsize,binstep);
                end
            end
        end
    elseif strcmp(AUCtype,'stimuli') %here, compare auc of cor/err for each stimuli
        [ITI, sound,delay, response, lick, pITI, psound, pdelay, presponse, plick]= deal(nan(nROI,size(trialType,2)));
        for roiNo = 1:size(SavedCaTrials.f_raw{1},1)
            disp([animal,date,'rioNo',num2str(roiNo)]);
            [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,dff(roiNo,:),frT,str_nFrames);
            label_choice = fTrialType2Label(trialType(1:2,:,:),1);%only include cor and err trials
            poslabel=1;
            nshuffle=1000;
            for nStim=1:size(trialType,2)
                ind_trial=logical(reshape(sum(trialType(1:2,nStim,:),1),[],1));
                [ITI(roiNo,nStim),pITI(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.ITI(ind_trial),poslabel,nshuffle);
                [sound(roiNo,nStim),psound(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.sound(ind_trial),poslabel,nshuffle);
                [delay(roiNo,nStim),pdelay(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.delay(ind_trial),poslabel,nshuffle);
                [response(roiNo,nStim),presponse(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.response(ind_trial),poslabel,nshuffle);
                [lick(roiNo,nStim),plick(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.lick(ind_trial),poslabel,nshuffle);
            end
            %calculate moving AUC
            [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff(roiNo,:), behEventFrameIndex,  frameNum );
            binsize=1;
            binstep=1;
            label = label_choice;
            for nStim=1:size(trialType,2)
                indTrial=sum(trialType(1:2,nStim,:),1);%only for cor and err trials(1st d)
                indTrial=logical(squeeze(indTrial));
                [delayMovingAUC{roiNo}.stim{nStim},pdelayMovingAUC{roiNo}.stim{nStim}] = fMovingAUC(label(indTrial),dff_aligned(indTrial,:),poslabel,nshuffle,binsize,binstep);
            end
        end
    end
    varROI=(1:size(SavedCaTrials.f_raw{1},1));
    varROI=reshape(varROI,[],1);
    TAUC=table(varanimal,vardate,varfield,varcelltype,varROI,ITI, sound,delay, response, ...
        lick, pITI, psound, pdelay, presponse, plick,delayMovingAUC,pdelayMovingAUC,'VariableNames',...
        {'animal','date','field','celltype','nROI','ITI','sound','delay','response',...
        'lick','pITI','psound','pdelay','presponse','plick','delayMovingAUC','pdelayMovingAUC'});
end
save(fileNameT,'TAUC');
end

