function [TSVM] = fGetEpochSVMscoreASession(filepath,animal,date,field,celltype,trialTypeStr,SVMtype,nRepeat, pTraining,CorrectedMethod)
%FGETEPOCHSVMSCOREASESSION calculate epoch SVM from mulitiple sessions,
%with variables- ITI,sound,delay,response,lick, etc. only for control
%session, not opto session
%Input-filepath, path of one imaging session where beh, dff exist
%   CorrectedMethod-method to correct AUC value if correct and error 
%   trials are combined to calculate sensory/choice
%   AUC£¬{'balencedCorErrTrialNum'£¬'SensoryChoiceOrthogonalSubtraction','raw'},
%   ref fGetEpochAUCtableASession
%Output- TSVM with varibles{animal, date, field, ROI, ITI, sound,
%delay, response, lick, pITI, psound, pdelay, presponse, plick, celltype};
cd(filepath);
dirmat=strcat(filepath,filesep,'*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
file_imaging=cellfun(@(x) contains(x,'CaTrialsSIM_'), filenames);
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

%naming the temp file by method of AUC calculation correction
if strcmp(trialTypeStr,'cor')
    combineCorErr='divideCorErr';%and only use correct trial
    corErrTrialNumber='raw';
    fileNameT=[filepath,filesep,animal,'-',datestr,'trialType',trialTypeStr,'-',SVMtype,'-timeBin',str_nFrames,'-pTraining',num2str(pTraining),'-EpochSVM.mat'];
elseif strcmp(trialTypeStr,'cor and err')%used to test whether it is sensory or choice AUC,in order to be more fair, keep trial numbers balenced for correct and error trials
    combineCorErr='combineCorErr';%{'combineCorErr','divideCorErr'}
    fileNameT=[filepath,filesep,animal,'-',datestr,'trialType',trialTypeStr,'-',SVMtype,'-',CorrectedMethod,'-timeBin',str_nFrames,'-pTraining',num2str(pTraining),'-EpochSVM.mat'];
    switch CorrectedMethod
        case 'balencedCorErrTrialNum'%probably problematic, since the trial are resampled and may let training and testing dataset be same
            corErrTrialNumber='balence';
        case 'SensoryChoiceOrthogonalSubtraction'
            corErrTrialNumber='raw';
            if (~strcmp(SVMtype,'choice'))&& (~strcmp(SVMtype,'sensory'))
                warning('error input combination of fGetEpochAUCtableASession function, SensoryChoiceOrthogonalSubtraction method');
            end
        otherwise
            corErrTrialNumber='raw';
    end     
end

if exist(fileNameT,'file')
    load(fileNameT);
    disp(['Table exist, use ',animal,date]);
    
else
    disp(['analyzing',fileNameT]);
    nROI=1;%each session was viewed as one ROI
    varNROI=size(SavedCaTrials.f_raw{1},1);
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
    
    if strcmp(SVMtype,'choice')|| strcmp(SVMtype,'sensory')
        [ITI, sound,delay, response,lick,pITI, psound, pdelay, presponse, plick] = deal(cell(ones(1,1)));
        switch corErrTrialNumber
            case 'balence'
                [trialType,~,~] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix','left','divideCorErr');
                trialTypeIndCell=cell(2,2);%{ipsi cor, contra cor;ipsi err, contra err};
                for i=1:2
                    for j=1:2
                        trialTypeIndCell{i,j}=find(logical(reshape(trialType(i,j,:),[],1)));
                    end
                end
                temp=cellfun(@length,trialTypeIndCell);
                trialNumEachType=floor(sum(temp,'all')/4);%keep total trial number stable and redistributed to each trial types
                s = RandStream('mlfg6331_64');%for reproducibility
                trialTypeIndCellBalenced=cellfun(@(x) reshape(datasample(s,x,trialNumEachType,'Replace',true),[],1), trialTypeIndCell,'UniformOutput',false);
                trialTypeIndFinal=cell2mat(trialTypeIndCellBalenced);%1st col-cor; 2nd col-err
                ind_trial=reshape(trialTypeIndFinal,[],1);
                label_SVM=[ones(trialNumEachType*2,1);ones(trialNumEachType*2,1)*2];
            case 'raw'
                [trialType,~,~] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix','left',combineCorErr);
                label_choice = fTrialType2Label(trialType,2);
                if contains(trialTypeStr,'cor') % 'cor'| 'cor and err'
                    ind_trial=logical(reshape(sum(trialType(1,:,:),2),[],1));%only correct trials
                end
                label_SVM=label_choice(ind_trial);
                %used for orthogonal subtraction
                [trialType_orthogonal,~,~] = fGetTrialType( Data_extract,[],7-trialTypeVar(i_selectivity),'matrix','left',combineCorErr);
                label_choice_orthogonal = fTrialType2Label(trialType_orthogonal,2);
                label_SVM_orthogonal=label_choice_orthogonal(ind_trial);
        end

        sigbyEpoch=[];
        sigMoving=[];
        sigMovingGo=[];
        for roiNo = 1:size(SavedCaTrials.f_raw{1},1)
            %data for epoch activity
            [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,dff(roiNo,:),frT,str_nFrames);
            if strcmp(CorrectedMethod,'SensoryChoiceOrthogonalSubtraction')
                T_SigbyEpoch=fOrthogonalSubtraction(T_SigbyEpoch,ind_trial,label_SVM,label_SVM_orthogonal);
            else
                T_SigbyEpoch=T_SigbyEpoch(ind_trial,:);
            end
            sigbyEpoch=cat(3,sigbyEpoch,table2array(T_SigbyEpoch));
            %data for moving activity
            [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff(roiNo,:), behEventFrameIndex,  frameNum );
            if strcmp(CorrectedMethod,'SensoryChoiceOrthogonalSubtraction')
                dff_aligned_ortho_corrected=fOrthogonalSubtraction(dff_aligned,ind_trial,label_SVM,label_SVM_orthogonal);
            else
                dff_aligned_ortho_corrected=dff_aligned(ind_trial,:);
            end
            sigMoving=cat(3,sigMoving,dff_aligned_ortho_corrected);
            %data for moving activity around go cue
            behEventAlign='go cue';
            [ dff_aligned_go, ~, ~ ] = fAlignSigBehEvent( dff(roiNo,:), behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );
            if strcmp(CorrectedMethod,'SensoryChoiceOrthogonalSubtraction')
                dff_aligned_ortho_corrected=fOrthogonalSubtraction(dff_aligned_go,ind_trial,label_SVM,label_SVM_orthogonal);
            else
                dff_aligned_ortho_corrected=dff_aligned_go(ind_trial,:);
            end
            sigMovingGo=cat(3,sigMovingGo,dff_aligned_ortho_corrected);
        end
        sigbyEpoch=permute(sigbyEpoch,[1,3,2]);
        sigMoving=permute(sigMoving,[1,3,2]);
        sigMovingGo=permute(sigMovingGo,[1,3,2]);
        
        score_shuffle=cell(1,1);
        [score,score_shuffle{1}] = fSVM(sigbyEpoch,label_SVM, nRepeat, pTraining,1,1); 
        [meanout,pout] = fBootstrpMeanP(score,1000,0.5);
        [ITI{1}, sound{1},delay{1}, response{1},lick{1}] = deal(score(:,1),score(:,2),score(:,3),score(:,4),score(:,5));
        [pITI{1}, psound{1},pdelay{1}, presponse{1},plick{1}] = deal(pout(:,1),pout(:,2),pout(:,3),pout(:,4),pout(:,5));
        %calculate moving SVM
        binsize=3;
        binstep=3;
        for nResult=1%:size(trialType,1)-2
            if nResult==1
                if strcmp(trialTypeStr,'cor and err')
                    [delayMovingSVM.do,delayMovingSVM.shuffle_do]=fSVM(sigMoving,label_SVM, nRepeat, pTraining,binsize,binstep); 
                    [~,pdelayMovingSVM.do]=fBootstrpMeanP(delayMovingSVM.do,1000,0.5);
                    [goMovingSVM.do,goMovingSVM.shuffle_do]=fSVM(sigMovingGo,label_SVM, nRepeat, pTraining,binsize,binstep); 
                    [~,pgoMovingSVM.do]=fBootstrpMeanP(goMovingSVM.do,1000,0.5);
                else
                    [delayMovingSVM.cor,delayMovingSVM.shuffle_cor]=fSVM(sigMoving,label_SVM, nRepeat, pTraining,binsize,binstep); 
                    [~,pdelayMovingSVM.cor]=fBootstrpMeanP(delayMovingSVM.cor,1000,0.5);
                    [goMovingSVM.cor,goMovingSVM.shuffle_cor]=fSVM(sigMovingGo,label_SVM, nRepeat, pTraining,binsize,binstep);
                    [~,pgoMovingSVM.cor]=fBootstrpMeanP(goMovingSVM.cor,1000,0.5);
                end
%             elseif nResult==2
%                 [delayMovingSVM.err,delayMovingSVM.shuffle_err]=fSVM(sigMoving,label_SVM, nRepeat, pTraining,binsize,binstep); 
%                 [~,pdelayMovingSVM.err]=fBootstrpMeanP(delayMovingSVM.err,1000,0.5);
%                 [goMovingSVM.err,goMovingSVM.shuffle_err]=fSVM(sigMovingGo,label_SVM, nRepeat, pTraining,binsize,binstep);
%                 [~,pgoMovingSVM.err]=fBootstrpMeanP(goMovingSVM.err,1000,0.5);
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
    
    TSVM=table(varanimal,vardate,varfield,varcelltype,varNROI,ITI, sound,delay, response, ...
        lick, pITI, psound, pdelay, presponse, plick,delayMovingSVM,pdelayMovingSVM,...
        goMovingSVM,pgoMovingSVM,score_shuffle,'VariableNames',...
        {'animal','date','field','celltype','nROI','ITI','sound','delay','response',...
        'lick','pITI','psound','pdelay','presponse','plick','delayMovingSVM','pdelayMovingSVM',...
        'goMovingSVM','pgoMovingSVM','shuffleEpochSVMscore'});
end
save(fileNameT,'TSVM');
end

function [meanout,pout] = fBootstrpMeanP(activity,nboot,baseline)
bootmean=bootstrp(nboot,@nanmean,activity);
meanout=nanmean(bootmean);
pout=sum(bootmean<baseline)./sum(~isnan(bootmean));
pout(pout<0.5)=pout(pout<0.5)*2;
pout(pout>=0.5)=(1-pout(pout>=0.5))*2;
end

function [Tout]=fOrthogonalSubtraction(Tin,ind_trial,label_SVM,label_SVM_orthogonal)
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
