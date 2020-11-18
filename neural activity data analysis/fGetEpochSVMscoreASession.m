function [TSVM] = fGetEpochSVMscoreASession(filepath,animal,date,field,celltype,trialTypeStr,SVMtype,nRepeat, pTraining)
%FGETEPOCHSVMSCOREASESSION calculate epoch SVM from mulitiple sessions,
%with variables- ITI,sound,delay,response,lick, etc. only for control
%session, not opto session
%Input-filepath, path of one imaging session where beh, dff exist
%Output- TSVM with varibles{animal, date, field, ROI, ITI, sound,
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
i_selectivity=cellfun(@(x) strcmp(x,SVMtype),selectivitystr);%*********variable**************
i_selectivity=find(i_selectivity);
str_nFrames='500ms';%'500ms';%'1s'
fileNameT=[filepath,filesep,animal,'-',datestr,'trialType',trialTypeStr,'-',SVMtype,'-timeBin',str_nFrames,'-pTraining',num2str(pTraining),'-EpochSVM.mat'];
if exist(fileNameT,'file')
    load(fileNameT);
    disp(['Table exist, use ',animal,date]);
    
else
    disp(['analyzing',fileNameT]);
    nROI=1;%each session was viewed as one ROI
    [delayMovingSVM,pdelayMovingSVM]=deal([]);
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
    
    if strcmp(SVMtype,'choice')|| strcmp(SVMtype,'sensory')
        [ITI, sound,delay, response,lick,pITI, psound, pdelay, presponse, plick] = deal(cell(ones(1,1)));
        label_choice = fTrialType2Label(trialType,2);
        if contains(trialTypeStr,'cor') % 'cor'| 'cor and err'
            ind_trial=logical(reshape(sum(trialType(1,:,:),2),[],1));%only correct trials
        end
        sigbyEpoch=[];
        sigMoving=[];
        for roiNo = 1:size(SavedCaTrials.f_raw{1},1)
            %data for epoch activity
            [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,dff(roiNo,:),frT,str_nFrames);
            sigbyEpoch=cat(3,sigbyEpoch,table2array(T_SigbyEpoch));
            %data for moving activity
            [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff(roiNo,:), behEventFrameIndex,  frameNum );
            sigMoving=cat(3,sigMoving,dff_aligned);
        end
        sigbyEpoch=permute(sigbyEpoch,[1,3,2]);
        sigMoving=permute(sigMoving,[1,3,2]);
        
        score_shuffle=cell(1,1);
        [score,score_shuffle{1}] = fSVM(sigbyEpoch(ind_trial,:,:),label_choice(ind_trial), nRepeat, pTraining,1,1); 
        [meanout,pout] = fBootstrpMeanP(score,1000,0.5);
        [ITI{1}, sound{1},delay{1}, response{1},lick{1}] = deal(score(:,1),score(:,2),score(:,3),score(:,4),score(:,5));
        [pITI{1}, psound{1},pdelay{1}, presponse{1},plick{1}] = deal(pout(:,1),pout(:,2),pout(:,3),pout(:,4),pout(:,5));
        %calculate moving SVM
        binsize=3;
        binstep=1;
        label = label_choice;
        for nResult=1:size(trialType,1)-2
            indTrial=trialType(nResult,:,:);
            indTrial=sum(squeeze(indTrial),1);
            indTrial=logical(squeeze(indTrial));
            if nResult==1
                if strcmp(trialTypeStr,'cor and err')
                    [delayMovingSVM.do,delayMovingSVM.shuffle_do]=fSVM(sigMoving(indTrial,:,:),label_choice(indTrial), nRepeat, pTraining,binsize,binstep); 
                    [~,pdelayMovingSVM.do]=fBootstrpMeanP(delayMovingSVM.do,1000,0.5);
                else
                    [delayMovingSVM.cor,delayMovingSVM.shuffle_cor]=fSVM(sigMoving(indTrial,:,:),label_choice(indTrial), nRepeat, pTraining,binsize,binstep); 
                    [~,pdelayMovingSVM.cor]=fBootstrpMeanP(delayMovingSVM.cor,1000,0.5);
                end
            elseif nResult==2
                [delayMovingSVM.err,delayMovingSVM.shuffle_err]=fSVM(sigMoving(indTrial,:,:),label_choice(indTrial), nRepeat, pTraining,binsize,binstep); 
                [~,pdelayMovingSVM.err]=fBootstrpMeanP(delayMovingSVM.err,1000,0.5);
            end
        end
        
    elseif strcmp(SVMtype,'stimuli') %here, compare auc of cor/err for each stimuli
        [ITI, sound,delay, response, lick, pITI, psound, pdelay, presponse, plick]= deal(nan(nROI,size(trialType,2)));
%         for roiNo = 1
%             [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,dff(roiNo,:),frT,str_nFrames);
%             label_choice = fTrialType2Label(trialType(1:2,:,:),1);%only include cor and err trials
%             poslabel=1;
%             nshuffle=1000;
%             for nStim=1:size(trialType,2)
%                 ind_trial=logical(reshape(sum(trialType(1:2,nStim,:),1),[],1));
%                 [ITI(roiNo,nStim),pITI(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.ITI(ind_trial),poslabel,nshuffle);
%                 [sound(roiNo,nStim),psound(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.sound(ind_trial),poslabel,nshuffle);
%                 [delay(roiNo,nStim),pdelay(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.delay(ind_trial),poslabel,nshuffle);
%                 [response(roiNo,nStim),presponse(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.response(ind_trial),poslabel,nshuffle);
%                 [lick(roiNo,nStim),plick(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.lick(ind_trial),poslabel,nshuffle);
%             end
%             %calculate moving AUC
%             [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff(roiNo,:), behEventFrameIndex,  frameNum );
%             binsize=1;
%             binstep=1;
%             label = label_choice;
%             for nStim=1:size(trialType,2)
%                 indTrial=sum(trialType(1:2,nStim,:),1);%only for cor and err trials(1st d)
%                 indTrial=logical(squeeze(indTrial));
%                 [delayMovingSVM.stim{nStim},pdelayMovingSVM.stim{nStim}] = fMovingAUC(label(indTrial),dff_aligned(indTrial,:),poslabel,nshuffle,binsize,binstep);
%             end
%         end
    end
    
    TSVM=table(varanimal,vardate,varfield,varcelltype,ITI, sound,delay, response, ...
        lick, pITI, psound, pdelay, presponse, plick,delayMovingSVM,pdelayMovingSVM,score_shuffle,'VariableNames',...
        {'animal','date','field','celltype','ITI','sound','delay','response',...
        'lick','pITI','psound','pdelay','presponse','plick','delayMovingSVM','pdelayMovingSVM','shuffleEpochSVMscore'});
end
save(fileNameT,'TSVM');
end

function [meanout,pout] = fBootstrpMeanP(activity,nboot,baseline)
bootmean=bootstrp(nboot,@nanmean,activity);
meanout=nanmean(bootmean);
pout=sum(bootmean<baseline)/sum(~isnan(bootmean));
pout(pout<0.5)=pout(pout<0.5)*2;
pout(pout>=0.5)=(1-pout(pout>=0.5))*2;
end

