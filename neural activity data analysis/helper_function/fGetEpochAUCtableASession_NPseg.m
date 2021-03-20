function [TAUC, Tmean] = fGetEpochAUCtableASession_NPseg(filepath,animal,date,field,celltype,trialTypeStr,AUCtype,varargin)
%FGETEPOCHAUCTABLEASESSION_NPSEG calculate epoch AUC from mulitiple sessions,
%with variables- ITI,sound,delay,response,lick, etc. only for control
%session, not opto session; each session one row, the values is neural
%peel, calculated as middle segments of whole field(6,7,10,11 of 16 segment)
%Input-filepath, path of one imaging session where beh, dff exist
%   varargin-method to correct AUC value if correct and error trials are
%   combined to calculate sensory/choice AUC
%Output- 
%   TAUC with varibles{animal, date, field, ROI, ITI, sound,
%   delay, response, lick, pITI, psound, pdelay, presponse, plick, celltype};
%   Tmean with similar variables, but ITI, sound etc. means mean
%   activities, while psound, etc. means p of significant mean difference
cd(filepath);
dirmat=strcat(filepath,filesep,'*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
file_imaging=cellfun(@(x) contains(x,'Ca'), filenames);
load(filenames{file_imaging});%load imaging data
datestr=strrep(date,'/','');
CurrFolder=pwd;
if exist('dff_NPseg.mat','file')
    load('dff_NPseg.mat');%load dff
else
    dff=fnx_getDff(CurrFolder,[animal,'_',datestr],'save figure','field_NPseg');
    if isempty(dff)
        warning(['dff is empty, so skip session-',animal,'_',datestr]);
        TAUC=[];
        Tmean=[];
        return;
    end       
end
file_beh=cellfun(@(x) contains(x,'Virables'), filenames);
load(filenames{file_beh});


selectivitystr={'stimuli','sensory difficulty','sensory','choice'};%sensory means grouping difficulties;
trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
i_selectivity=cellfun(@(x) strcmp(x,AUCtype),selectivitystr);%*********variable**************
i_selectivity=find(i_selectivity);
str_nFrames='500ms';%'500ms';%'1s'
if strcmp(trialTypeStr,'cor')
    combineCorErr='divideCorErr';%and only use correct trial
    corErrTrialNumber='raw';
    fileNameT=[filepath,filesep,animal,'-',datestr,'trialType',trialTypeStr,'-',AUCtype,'-trialNum',corErrTrialNumber,'-timeBin',str_nFrames,'-EpochAUC_NPseg.mat'];
elseif strcmp(trialTypeStr,'cor and err')%used to test whether it is sensory or choice AUC,in order to be more fair, keep trial numbers balenced for correct and error trials
    combineCorErr='combineCorErr';%{'combineCorErr','divideCorErr'}
    if ~isempty(varargin) && strcmp(varargin{1},'balencedCorErrTrialNum')  
        corErrTrialNumber='balence';
        fileNameT=[filepath,filesep,animal,'-',datestr,'trialType',trialTypeStr,'-',AUCtype,'-trialNum',corErrTrialNumber,'-timeBin',str_nFrames,'-EpochAUC_NPseg.mat'];
    elseif ~isempty(varargin) && strcmp(varargin{1},'SensoryChoiceOrthogonalSubtraction')  
        corErrTrialNumber='raw';
        if strcmp(AUCtype,'choice')|| strcmp(AUCtype,'sensory')
            fileNameT=[filepath,filesep,animal,'-',datestr,'trialType',trialTypeStr,'-',AUCtype,'-AUCCorrectedMethod',varargin{1},'-timeBin',str_nFrames,'-EpochAUC_NPseg.mat'];
        else
            warning('error input combination of fGetEpochAUCtableASession function, SensoryChoiceOrthogonalSubtraction method');
        end
    else
        corErrTrialNumber='raw';
        fileNameT=[filepath,filesep,animal,'-',datestr,'trialType',trialTypeStr,'-',AUCtype,'-trialNum',corErrTrialNumber,'-timeBin',str_nFrames,'-EpochAUC_NPseg.mat'];
    end
    
end

if exist(fileNameT,'file')
    load(fileNameT);
    n_ROI=size(TAUC,1);
    disp(['Table exist, use ',animal,date,';nROI=',num2str(n_ROI)]);   
else
    nROI=1;%each field one value
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
    ntr=length(SavedCaTrials.SegNPdataAll);
    frT = SavedCaTrials.FrameTime;
    % align to behavior event
    nFrameEachTrial=cellfun(@(x) size(x,2),SavedCaTrials.SegNPdataAll);
    ind_1stFrame=zeros(1,length(nFrameEachTrial));
    ind_1stFrame(1)=1;
    ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
    ind_1stFrame=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials
    frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
    frameNum=double(round(frameNumTime*1000/frT));
    [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( Data_extract, ind_1stFrame, SavedCaTrials.FrameTime );%get behavior event time
    
    if strcmp(AUCtype,'choice')|| strcmp(AUCtype,'sensory')
        [ITI, sound,delay, response, lick, pITI, psound, pdelay, presponse, plick]= deal(nan(nROI,1));
        if strcmp(corErrTrialNumber,'balence')
            [trialType_raw,~,~] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix','left','divideCorErr');
            trialTypeIndCell=cell(2,2);%{ipsi cor, contra cor;ipsi err, contra err};
            for i=1:2
                for j=1:2
                    trialTypeIndCell{i,j}=find(logical(reshape(trialType_raw(i,j,:),[],1)));
                end
            end
            temp=cellfun(@length,trialTypeIndCell);
            trialNumEachType=floor(sum(temp,'all')/4);%keep total trial number stable and redistributed to each trial types
            s = RandStream('mlfg6331_64');%for reproducibility
            trialTypeIndCellBalenced=cellfun(@(x) reshape(datasample(s,x,trialNumEachType,'Replace',true),[],1), trialTypeIndCell,'UniformOutput',false);
            trialTypeIndFinal=cell2mat(trialTypeIndCellBalenced);%1st col-cor; 2nd col-err
            ind_trial=reshape(trialTypeIndFinal,[],1);
            label_AUC=[ones(trialNumEachType*2,1);ones(trialNumEachType*2,1)*2];
        elseif strcmp(corErrTrialNumber,'raw')
            [trialType,~,~] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix','left',combineCorErr);
            label_choice = fTrialType2Label(trialType,2);
            if contains(trialTypeStr,'cor') % 'cor'| 'cor and err'
                ind_trial=logical(reshape(sum(trialType(1,:,:),2),[],1));%only correct trials
            end
            label_AUC=label_choice(ind_trial);
            %used for orthogonal subtraction
            [trialType_orthogonal,~,~] = fGetTrialType( Data_extract,[],7-trialTypeVar(i_selectivity),'matrix','left',combineCorErr);
            label_choice_orthogonal = fTrialType2Label(trialType_orthogonal,2);
            label_AUC_orthogonal=label_choice_orthogonal(ind_trial);
        end
        for roiNo = 1
            disp([animal,date,'neuro peel']);
            [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,dff(roiNo,:),frT,str_nFrames);
            if ~isempty(varargin) && strcmp(varargin{1},'SensoryChoiceOrthogonalSubtraction')
                T_SigbyEpoch=fOrthogonalSubtraction(T_SigbyEpoch,ind_trial,label_AUC,label_AUC_orthogonal);
            else
                T_SigbyEpoch=T_SigbyEpoch(ind_trial,:);
            end
            poslabel=2;
            nshuffle=1000;%%%%check
            [ITI(roiNo),pITI(roiNo)]=fAUC(label_AUC,T_SigbyEpoch.ITI,poslabel,nshuffle);
            [sound(roiNo),psound(roiNo)]=fAUC(label_AUC,T_SigbyEpoch.sound,poslabel,nshuffle);
            [delay(roiNo),pdelay(roiNo)]=fAUC(label_AUC,T_SigbyEpoch.delay,poslabel,nshuffle);
            [response(roiNo),presponse(roiNo)]=fAUC(label_AUC,T_SigbyEpoch.response,poslabel,nshuffle);
            [lick(roiNo),plick(roiNo)]=fAUC(label_AUC,T_SigbyEpoch.lick,poslabel,nshuffle);
            %calculate moving AUC
            [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff(roiNo,:), behEventFrameIndex,  frameNum );
            if ~isempty(varargin) && strcmp(varargin{1},'SensoryChoiceOrthogonalSubtraction')
                dff_aligned_ortho_corrected=fOrthogonalSubtraction(dff_aligned,ind_trial,label_AUC,label_AUC_orthogonal);
            end
            binsize=1;
            binstep=1;
            for nResult=1:size(trialType,1)-2
                if strcmp(trialTypeStr,'cor and err')%combine cor and err together
                    indTrial=ind_trial;
                else %calculate AUC for cor and err respectively
                    indTrial=trialType(nResult,:,:);
                    indTrial=sum(squeeze(indTrial),1);
                    indTrial=logical(squeeze(indTrial));
                end
                if nResult==1
                    if strcmp(trialTypeStr,'cor and err')
                        [delayMovingAUC{roiNo}.do,pdelayMovingAUC{roiNo}.do] = fMovingAUC(label_AUC,dff_aligned_ortho_corrected,2,nshuffle,binsize,binstep);
                    else
                        [delayMovingAUC{roiNo}.cor,pdelayMovingAUC{roiNo}.cor] = fMovingAUC(label_choice,dff_aligned,2,nshuffle,binsize,binstep);
                    end
                elseif nResult==2
                    [delayMovingAUC{roiNo}.err,pdelayMovingAUC{roiNo}.err] = fMovingAUC(label_choice,dff_aligned,2,nshuffle,binsize,binstep);
                end
            end
        end
       
    elseif strcmp(AUCtype,'stimuli') %here, compare auc of cor/err for each stimuli
        [trialType,~,~] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix','left','divideCorErr');
        [ITI, sound,delay, response, lick, pITI, psound, pdelay, presponse, plick]= deal(nan(nROI,size(trialType,2)));
        for roiNo = 1
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
    

    TAUC=table(varanimal,vardate,varfield,varcelltype,ITI, sound,delay, response, ...
        lick, pITI, psound, pdelay, presponse, plick,delayMovingAUC,pdelayMovingAUC,'VariableNames',...
        {'animal','date','field','celltype','ITI','sound','delay','response',...
        'lick','pITI','psound','pdelay','presponse','plick','delayMovingAUC','pdelayMovingAUC'});
end
%calculate mean activities for each time epoch and whether they
%different from ITI as baseline
if ~exist('Tmean','var')
    nROI=1;
%     [ITI_m, sound_m,delay_m, response_m, lick_m, pITI_m, psound_m, pdelay_m, presponse_m, plick_m]= deal(nan(nROI,1));
    [trialType,~,~] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix','left',combineCorErr);
    [trialType_byChoice,~,~] = fGetTrialType( Data_extract,[],3,'matrix','left','combineCorErr');
    if contains(trialTypeStr,'cor') % 'cor'| 'cor and err'
        ind_trial=logical(reshape(sum(trialType(1,:,:),2),[],1));%only correct trials
    end
    ind_ipsi=logical(reshape(sum(trialType_byChoice(1,1,:),2),[],1));
    ind_contra=logical(reshape(sum(trialType_byChoice(1,2,:),2),[],1));
    ind_tr_1=1;
    ntr=length(SavedCaTrials.SegNPdataAll);
    frT = SavedCaTrials.FrameTime;
     % align to behavior event
    nFrameEachTrial=cellfun(@(x) size(x,2),SavedCaTrials.SegNPdataAll);
    ind_1stFrame=zeros(1,length(nFrameEachTrial));
    ind_1stFrame(1)=1;
    ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
    ind_1stFrame=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials
    frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
    frameNum=double(round(frameNumTime*1000/frT));
    [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( Data_extract, ind_1stFrame, SavedCaTrials.FrameTime );%get behavior event time    
    for roiNo = 1 %SavedCaTrials.nROIs may be not true
        [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,dff(roiNo,:),frT,str_nFrames);
        T_SigbyEpoch_mean=T_SigbyEpoch(ind_trial,:);
        [ITI_m(roiNo).mean,pITI_m(roiNo).mean]=fMeanPdiff(T_SigbyEpoch_mean.ITI,T_SigbyEpoch_mean.ITI);
        [sound_m(roiNo).mean, psound_m(roiNo).mean]=fMeanPdiff(T_SigbyEpoch_mean.ITI,T_SigbyEpoch_mean.sound);
        [delay_m(roiNo).mean,pdelay_m(roiNo).mean]=fMeanPdiff(T_SigbyEpoch_mean.ITI,T_SigbyEpoch_mean.delay);
        [response_m(roiNo).mean,presponse_m(roiNo).mean]=fMeanPdiff(T_SigbyEpoch_mean.ITI,T_SigbyEpoch_mean.response);
        [lick_m(roiNo).mean,plick_m(roiNo).mean]=fMeanPdiff(T_SigbyEpoch_mean.ITI,T_SigbyEpoch_mean.lick);
        T_SigbyEpoch_ipsi=T_SigbyEpoch(logical(ind_trial.*ind_ipsi),:);
        [ITI_m(roiNo).ipsi,pITI_m(roiNo).ipsi]=fMeanPdiff(T_SigbyEpoch_ipsi.ITI,T_SigbyEpoch_ipsi.ITI);
        [sound_m(roiNo).ipsi, psound_m(roiNo).ipsi]=fMeanPdiff(T_SigbyEpoch_ipsi.ITI,T_SigbyEpoch_ipsi.sound);
        [delay_m(roiNo).ipsi,pdelay_m(roiNo).ipsi]=fMeanPdiff(T_SigbyEpoch_ipsi.ITI,T_SigbyEpoch_ipsi.delay);
        [response_m(roiNo).ipsi,presponse_m(roiNo).ipsi]=fMeanPdiff(T_SigbyEpoch_ipsi.ITI,T_SigbyEpoch_ipsi.response);
        [lick_m(roiNo).ipsi,plick_m(roiNo).ipsi]=fMeanPdiff(T_SigbyEpoch_ipsi.ITI,T_SigbyEpoch_ipsi.lick);
        T_SigbyEpoch_contra=T_SigbyEpoch(logical(ind_trial.*ind_contra),:);
        [ITI_m(roiNo).contra,pITI_m(roiNo).contra]=fMeanPdiff(T_SigbyEpoch_contra.ITI,T_SigbyEpoch_contra.ITI);
        [sound_m(roiNo).contra, psound_m(roiNo).contra]=fMeanPdiff(T_SigbyEpoch_contra.ITI,T_SigbyEpoch_contra.sound);
        [delay_m(roiNo).contra,pdelay_m(roiNo).contra]=fMeanPdiff(T_SigbyEpoch_contra.ITI,T_SigbyEpoch_contra.delay);
        [response_m(roiNo).contra,presponse_m(roiNo).contra]=fMeanPdiff(T_SigbyEpoch_contra.ITI,T_SigbyEpoch_contra.response);
        [lick_m(roiNo).contra,plick_m(roiNo).contra]=fMeanPdiff(T_SigbyEpoch_contra.ITI,T_SigbyEpoch_contra.lick);
    end
    [varanimal{1:nROI}]=deal(animal);
    varanimal=reshape(varanimal,[],1);
    [vardate{1:nROI}]=deal(date);
    vardate=reshape(vardate,[],1);
    [varfield{1:nROI}]=deal(field);
    varfield=reshape(varfield,[],1);
    [varcelltype{1:nROI}]=deal(celltype);
    varcelltype=reshape(varcelltype,[],1);

    Tmean=table(varanimal,vardate,varfield,varcelltype,...
        ITI_m', sound_m',delay_m', response_m', lick_m', ...
        pITI_m', psound_m', pdelay_m', presponse_m', plick_m','VariableNames',...
        {'animal','date','field','celltype',...
        'ITI','sound','delay','response','lick',...
        'pITI','psound','pdelay','presponse','plick'});
end
save(fileNameT,'TAUC','Tmean');
end


function [Tout]=fOrthogonalSubtraction(Tin,ind_trial,label_AUC,label_AUC_orthogonal)
if istable(Tin)
    Tout=table2array(Tin(ind_trial,:));
else
    Tout=Tin(ind_trial,:);
end
n_condition=length(unique(label_AUC));
n_condition_orth=length(unique(label_AUC_orthogonal));
for i=1:n_condition_orth
    mean_orth=zeros(n_condition,size(Tout,2));
    for j=1:n_condition
        indTrial=logical((label_AUC==j).*(label_AUC_orthogonal==i));
        mean_orth(j,:)=nanmean(Tout(indTrial,:));   
    end
    indTrial_orth=(label_AUC_orthogonal==i);
    Tout(indTrial_orth,:)=Tout(indTrial_orth,:)-nanmean(mean_orth);
end
if istable(Tin)
    Tout=array2table(Tout,'VariableNames',Tin.Properties.VariableNames);
end
end

function [mean, p] = fMeanPdiff(baseline,data,varargin)
%Input- matrix of baseline/data, n-by-m (n repeats, m time points)
%   method- 'ranksum'(default)|'ttest', etc
%Output- vector of mean, p, 1-by-m
if isempty(varargin)
    method ='ranksum';
else
    method =varargin{1};
end
mean=nanmean(data,1);
switch method
    case 'ranksum'
        if size(data,2)~=size(baseline,2)
            warnning('input data not same column number');
            return;
        end
        for i=1:size(data,2)
            p(i)=ranksum(baseline(:,i), data(:,i));
        end
    case 'ttest'
        [~,p]=ttest(baseline,data);
end
end
