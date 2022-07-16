classdef SessionOpto2P < Session2P
    %SESSIONOPTO2P that subclass the Session2P class
    %   Detailed explanation goes here
    
    properties
        TSVM_ctrl   %table storing the cross epoch and time SVM, default is control
        TAUC_ctrl   %table storing the cross epoch and time AUC
        Tmean_ctrl  %table storing the cross epoch and time mean activities
        TSVM_opto   %table storing the cross epoch and time SVM for opto trials
        TAUC_opto   %table storing the cross epoch and time AUC for opto trials
        Tmean_opto  %table storing the cross epoch and time mean activities for opto trials
    end
    
    methods
        function obj = SessionOpto2P(sessionname,rootpath,trial2include,trial2exclude,varargin)
            %SESSIONOPTO2P Construct an instance of this class
            %   Detailed explanation goes here
            if ~isempty(varargin)
                ind_tr_1=varargin{1};
            else
                ind_tr_1=1;
            end
            obj = obj@Session2P(sessionname,rootpath,trial2include,trial2exclude,ind_tr_1);
            indOpto=reshape(obj.Data_extract.Opto_trial_index,[],1);
            obj.metadata.indOpto=logical(indOpto.*obj.metadata.indTrial2use);
            obj.metadata.indCtrl=logical((~indOpto).*obj.metadata.indTrial2use);
            
        end
        
        function [TAUC, Tmean,obj] = mGetEpochAUC(obj,celltype,ROItype,trialTypeStr,...
                AUCtype,AUCCorrectedMethod,varargin)
            if ~isempty(varargin)
                str_trialtype=varargin{1};
                if strcmp(str_trialtype,'ctrl')
                    ind_chosen=obj.metadata.indCtrl;
                elseif strcmp(str_trialtype,'opto')
                    ind_chosen=obj.metadata.indOpto;
                end
            else
                ind_chosen=obj.metadata.indTrial2use;
                str_trialtype='ctrl';
            end
            date_str=datestr(obj.metadata.datetime,'yyyy-mm-dd');
            date_str=strrep(date_str,'-0','-');
            selectivitystr={'stimuli','sensory difficulty','sensory','choice'};%sensory means grouping difficulties;
            trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
            i_selectivity=cellfun(@(x) strcmp(x,AUCtype),selectivitystr);%*********variable**************
            i_selectivity=find(i_selectivity);
            str_nFrames='1s';%'500ms';%'1s'
            %naming the temp file by method of AUC calculation correction
            filepath_inSummary='E:\2P\summary\AUC\cases';%backup another copy in the summary folder
            fileName_datapars=[obj.metadata.animal,'-',date_str,'trialType',trialTypeStr,'-',str_trialtype];
            if strcmp(trialTypeStr,'cor') || strcmp(trialTypeStr,'err')
                combineCorErr='divideCorErr';%and only use correct trial
                corErrTrialNumber='raw';
                fileName_processpars=[AUCtype,'-trialNum',corErrTrialNumber,'-timeBin',str_nFrames];
            elseif strcmp(trialTypeStr,'cor and err')%used to test whether it is sensory or choice AUC,in order to be more fair, keep trial numbers balenced for correct and error trials
                combineCorErr='combineCorErr';%{'combineCorErr','divideCorErr'}
                switch AUCCorrectedMethod
                    case 'balencedCorErrTrialNum'
                        corErrTrialNumber='balence';
                        fileName_processpars=[AUCtype,'-trialNum',corErrTrialNumber,'-timeBin',str_nFrames];
                    case 'SensoryChoiceOrthogonalSubtraction'
                        corErrTrialNumber='raw';
                        if strcmp(AUCtype,'choice')|| strcmp(AUCtype,'sensory')
                            fileName_processpars=[AUCtype,'-AUCCorrectedMethod',AUCCorrectedMethod,'-timeBin',str_nFrames];
                        else
                            warning('error input combination of fGetEpochAUCtableASession function, SensoryChoiceOrthogonalSubtraction method');
                        end
                    otherwise
                        corErrTrialNumber='raw';
                        fileName_processpars=[AUCtype,'-trialNum',corErrTrialNumber,'-timeBin',str_nFrames];
                end
            end
            fileNameT=[obj.path.root,filesep,fileName_datapars,'-',fileName_processpars,'-EpochAUC.mat'];
            fileNameT_inSummary=[filepath_inSummary,filesep,fileName_datapars,'-',fileName_processpars,'-EpochAUC.mat'];
            
            if exist(fileNameT,'file')%CD058-2018-1-27trialTypecor and err-choice-AUCCorrectedMethodSensoryChoiceOrthogonalSubtraction-timeBin1s-EpochAUC
                load(fileNameT);
                n_ROI=size(TAUC,1);
                disp(['Table exist, use ',obj.metadata.animal,date_str,';nROI=',num2str(n_ROI)]);
                %delete duplicated variables
                standard_varNames={'animal','date','field','celltype','ROItype','nROI',...
                    'ipsi_performance', 'contra_performance','overall_performance',...
                    'ITI','sound','delay','response','lick','mid_delay','late_delay',...
                    'pITI','psound','pdelay','presponse','plick','pmid_delay','plate_delay','delayMovingAUC','pdelayMovingAUC'};
                T_names=TAUC.Properties.VariableNames;
                for i=1:length(T_names)
                    temp=cellfun(@(x) strcmp(x,T_names(i)), standard_varNames);
                    if sum(temp)==0
                        TAUC(:,T_names{i})=[];
                    end
                end
            else
                nROI=size(obj.SavedCaTrials.f_raw{1},1);
                [delayMovingAUC,pdelayMovingAUC]=deal(cell(nROI,1));
                [varanimal{1:nROI}]=deal(obj.metadata.animal);
                varanimal=reshape(varanimal,[],1);
                [vardate{1:nROI}]=deal(datestr(obj.metadata.datetime,'yyyy/mm/dd'));
                vardate=reshape(vardate,[],1);
                [varfield{1:nROI}]=deal(obj.metadata.field);
                varfield=reshape(varfield,[],1);
                [varcelltype{1:nROI}]=deal(celltype);
                varcelltype=reshape(varcelltype,[],1);
                [varROItype{1:nROI}]=deal(ROItype);
                varROItype=reshape(varROItype,[],1);
                ind_tr_1=1;
                ntr=length(obj.SavedCaTrials.f_raw);
                frT = obj.SavedCaTrials.FrameTime;
                % align to behavior event
                nFrameEachTrial=cellfun(@(x) size(x,2),obj.SavedCaTrials.f_raw);
                ind_1stFrame=zeros(1,length(nFrameEachTrial));
                ind_1stFrame(1)=1;
                ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
                ind_1stFrame=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials
                frameNumTime=[1,1.5];%from 1s before delay onset to 1.5s after delay onset
                frameNum=double(round(frameNumTime*1000/frT));
                [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( obj.Data_extract, ind_1stFrame, obj.SavedCaTrials.FrameTime );%get behavior event time
                
                if strcmp(AUCtype,'choice')|| strcmp(AUCtype,'sensory')
                    [ITI, sound,delay, response, lick,mid_delay,late_delay, pITI, psound, pdelay, presponse, plick,pmid_delay,plate_delay]= deal(nan(nROI,1));
                    switch corErrTrialNumber
                        case 'balence'
                            [trialType,~,~] = fGetTrialType( obj.Data_extract,[],trialTypeVar(i_selectivity),'matrix','left','divideCorErr');
                            trialTypeIndCell=cell(2,2);%{ipsi cor, contra cor;ipsi err, contra err};
                            for i=1:2
                                for j=1:2
                                    trialTypeIndCell{i,j}=find(logical(reshape(trialType(i,j,:),[],1).*ind_chosen));
                                end
                            end
                            temp=cellfun(@length,trialTypeIndCell);
                            trialNumEachType=floor(sum(temp,'all')/4);%keep total trial number stable and redistributed to each trial types
                            s = RandStream('mlfg6331_64');%for reproducibility
                            trialTypeIndCellBalenced=cellfun(@(x) reshape(datasample(s,x,trialNumEachType,'Replace',true),[],1), trialTypeIndCell,'UniformOutput',false);
                            trialTypeIndFinal=cell2mat(trialTypeIndCellBalenced);%1st col-cor; 2nd col-err
                            ind_trial=reshape(trialTypeIndFinal,[],1);
                            label_AUC=[ones(trialNumEachType*2,1);ones(trialNumEachType*2,1)*2];
                        case 'raw'
                            [trialType,~,~] = fGetTrialType( obj.Data_extract,[],trialTypeVar(i_selectivity),'matrix','left',combineCorErr);
                            label_choice = fTrialType2Label(trialType,2);
                            if contains(trialTypeStr,'cor') % 'cor'| 'cor and err'
                                ind_trial=logical(reshape(sum(trialType(1,:,:),2),[],1));%only correct trials
                            end
                            if strcmp(trialTypeStr,'err')
                                ind_trial=logical(reshape(sum(trialType(2,:,:),2),[],1));%only error trials
                            end
                            ind_trial=logical(ind_chosen.*ind_trial);%include trials by manual define
                            label_AUC=label_choice(ind_trial);
                            %used for orthogonal subtraction
                            [trialType_orthogonal,~,~] = fGetTrialType( obj.Data_extract,[],7-trialTypeVar(i_selectivity),'matrix','left',combineCorErr);
                            label_choice_orthogonal = fTrialType2Label(trialType_orthogonal,2);
                            label_AUC_orthogonal=label_choice_orthogonal(ind_trial);
                    end
                    for roiNo = 1:size(obj.SavedCaTrials.f_raw{1},1) %SavedCaTrials.nROIs may be not true
                        disp([obj.metadata.animal,obj.metadata.date,'rioNo',num2str(roiNo)]);
                        [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,obj.dff(roiNo,:),frT,str_nFrames);
                        if strcmp(AUCCorrectedMethod,'SensoryChoiceOrthogonalSubtraction')
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
                        [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( obj.dff(roiNo,:), behEventFrameIndex,  frameNum );
                        dff_late_delay=nanmean(dff_aligned(:,frameNum(1)+round(1*1000/frT):frameNum(1)+round(1.5*1000/frT)),2);
                        dff_mid_delay=nanmean(dff_aligned(:,frameNum(1)+round(0.3*1000/frT):frameNum(1)+round(1*1000/frT)),2);
                        if strcmp(AUCCorrectedMethod,'SensoryChoiceOrthogonalSubtraction')
                            dff_aligned_ortho_corrected=fOrthogonalSubtraction(dff_aligned,ind_trial,label_AUC,label_AUC_orthogonal);
                            dff_mid_delay4auc=fOrthogonalSubtraction(dff_mid_delay,ind_trial,label_AUC,label_AUC_orthogonal);
                            dff_late_delay4auc=fOrthogonalSubtraction(dff_late_delay,ind_trial,label_AUC,label_AUC_orthogonal);
                        else
                            dff_aligned=dff_aligned(ind_trial,:);
                            dff_mid_delay4auc=dff_mid_delay(ind_trial,:);
                            dff_late_delay4auc=dff_late_delay(ind_trial,:);
                        end
                        [mid_delay(roiNo),pmid_delay(roiNo)]=fAUC(label_AUC,dff_mid_delay4auc,poslabel,nshuffle);
                        [late_delay(roiNo),plate_delay(roiNo)]=fAUC(label_AUC,dff_late_delay4auc,poslabel,nshuffle);
                        
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
                                    if strcmp(AUCCorrectedMethod,'SensoryChoiceOrthogonalSubtraction')
                                        [delayMovingAUC{roiNo}.do,pdelayMovingAUC{roiNo}.do] = fMovingAUC(label_AUC,dff_aligned_ortho_corrected,2,nshuffle,binsize,binstep);
                                    else
                                        [delayMovingAUC{roiNo}.do,pdelayMovingAUC{roiNo}.do] = fMovingAUC(label_AUC,dff_aligned,2,nshuffle,binsize,binstep);
                                    end
                                elseif strcmp(trialTypeStr,'cor')
                                    [delayMovingAUC{roiNo}.cor,pdelayMovingAUC{roiNo}.cor] = fMovingAUC(label_AUC,dff_aligned,2,nshuffle,binsize,binstep);
                                end
                            elseif nResult==2 && strcmp(trialTypeStr,'err')
                                [delayMovingAUC{roiNo}.err,pdelayMovingAUC{roiNo}.err] = fMovingAUC(label_AUC,dff_aligned,2,nshuffle,binsize,binstep);
                            end
                        end
                    end
                    
                elseif strcmp(AUCtype,'stimuli') %here, compare auc of cor/err for each stimuli
                    [trialType,~,~] = fGetTrialType( obj.Data_extract,[],trialTypeVar(i_selectivity),'matrix','left','divideCorErr');
                    [ITI, sound,delay, response, lick, pITI, psound, pdelay, presponse, plick]= deal(nan(nROI,size(trialType,2)));
                    for roiNo = 1:size(obj.SavedCaTrials.f_raw{1},1)
                        disp([animal,date,'rioNo',num2str(roiNo)]);
                        [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,obj.dff(roiNo,:),frT,str_nFrames);
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
                        [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( obj.dff(roiNo,:), behEventFrameIndex,  frameNum );
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
                %add task performance
                nROI=size(obj.SavedCaTrials.f_raw{1},1);
                
                Tperformance=obj.mSummaryBehavior(3);
                [ipsi_performance{1:nROI}]=deal(Tperformance{'ipsi','correct'});
                ipsi_performance=reshape(ipsi_performance,[],1);
                [contra_performance{1:nROI}]=deal(Tperformance{'contra','correct'});
                contra_performance=reshape(contra_performance,[],1);
                [overall_performance{1:nROI}]=deal(Tperformance{'total','correct'});
                overall_performance=reshape(overall_performance,[],1);
                varROI=(1:size(obj.SavedCaTrials.f_raw{1},1));
                varROI=reshape(varROI,[],1);
                TAUC=table(varanimal,vardate,varfield,varcelltype,varROItype,varROI,...
                    ipsi_performance, contra_performance,overall_performance,...
                    ITI, sound,delay, response, lick, mid_delay,late_delay,...
                    pITI, psound, pdelay, presponse, plick,pmid_delay,plate_delay,delayMovingAUC,pdelayMovingAUC,...
                    'VariableNames',{'animal','date','field','celltype','ROItype','nROI',...
                    'ipsi_performance', 'contra_performance','overall_performance',...
                    'ITI','sound','delay','response','lick','mid_delay','late_delay',...
                    'pITI','psound','pdelay','presponse','plick','pmid_delay','plate_delay','delayMovingAUC','pdelayMovingAUC'});
            end
            %calculate mean activities for each time epoch and whether they
            %different from ITI as baseline
            if ~exist('Tmean','var')
                nROI=size(obj.SavedCaTrials.f_raw{1},1);
                %     [ITI_m, sound_m,delay_m, response_m, lick_m, pITI_m, psound_m, pdelay_m, presponse_m, plick_m]= deal(nan(nROI,1));
                [trialType,~,~] = fGetTrialType( obj.Data_extract,[],trialTypeVar(i_selectivity),'matrix','left',combineCorErr);
                [trialType_byChoice,~,~] = fGetTrialType( obj.Data_extract,[],3,'matrix','left','combineCorErr');
                if contains(trialTypeStr,'cor') % 'cor'| 'cor and err'
                    ind_trial=logical(reshape(sum(trialType(1,:,:),2),[],1));%only correct trials
                end
                ind_ipsi=logical(reshape(sum(trialType_byChoice(1,1,:),2),[],1));
                ind_contra=logical(reshape(sum(trialType_byChoice(1,2,:),2),[],1));
                ind_trial=logical(ind_chosen.*ind_trial);
                ind_ipsi=logical(ind_chosen.*ind_ipsi);
                ind_contra=logical(ind_chosen.*ind_contra);
                ind_tr_1=1;
                ntr=length(obj.SavedCaTrials.f_raw);
                frT = obj.SavedCaTrials.FrameTime;
                % align to behavior event
                nFrameEachTrial=cellfun(@(x) size(x,2),obj.SavedCaTrials.f_raw);
                ind_1stFrame=zeros(1,length(nFrameEachTrial));
                ind_1stFrame(1)=1;
                ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
                ind_1stFrame=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials
                frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
                frameNum=double(round(frameNumTime*1000/frT));
                [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( obj.Data_extract, ind_1stFrame, obj.SavedCaTrials.FrameTime );%get behavior event time
                [sound_pselectivity,delay_pselectivity,response_pselectivity,lick_pselectivity]=deal(ones(size(obj.SavedCaTrials.f_raw{1},1),1));
                for roiNo = 1:size(obj.SavedCaTrials.f_raw{1},1) %SavedCaTrials.nROIs may be not true
                    [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,obj.dff(roiNo,:),frT,str_nFrames);
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
                    sound_pselectivity(roiNo)=ranksum(T_SigbyEpoch_contra.sound,T_SigbyEpoch_ipsi.sound);
                    delay_pselectivity(roiNo)=ranksum(T_SigbyEpoch_contra.delay,T_SigbyEpoch_ipsi.delay);
                    response_pselectivity(roiNo)=ranksum(T_SigbyEpoch_contra.response,T_SigbyEpoch_ipsi.response);
                    lick_pselectivity(roiNo)=ranksum(T_SigbyEpoch_contra.lick,T_SigbyEpoch_ipsi.lick);
                end
                [varanimal{1:nROI}]=deal(obj.metadata.animal);
                varanimal=reshape(varanimal,[],1);
                [vardate{1:nROI}]=deal(obj.metadata.date);
                vardate=reshape(vardate,[],1);
                [varfield{1:nROI}]=deal(obj.metadata.field);
                varfield=reshape(varfield,[],1);
                [varcelltype{1:nROI}]=deal(celltype);
                varcelltype=reshape(varcelltype,[],1);
                [varROItype{1:nROI}]=deal(ROItype);
                varROItype=reshape(varROItype,[],1);
                varROI=(1:size(obj.SavedCaTrials.f_raw{1},1));
                varROI=reshape(varROI,[],1);
                
                Tmean=table(varanimal,vardate,varfield,varcelltype,varROItype,varROI,...
                    ITI_m', sound_m',delay_m', response_m', lick_m', ...
                    pITI_m', psound_m', pdelay_m', presponse_m', plick_m',...
                    sound_pselectivity,delay_pselectivity,response_pselectivity,...
                    lick_pselectivity,'VariableNames',...
                    {'animal','date','field','celltype','ROItype','nROI',...
                    'ITI','sound','delay','response','lick',...
                    'pITI','psound','pdelay','presponse','plick',...
                    'sound_pselectivity','delay_pselectivity','response_pselectivity',...
                    'lick_pselectivity'});
            end
            switch str_trialtype
                case 'ctrl'
                    obj.TAUC_ctrl=TAUC;
                    obj.Tmean_ctrl=Tmean;
                    obj.TAUC=TAUC;
                    obj.Tmean=Tmean;
                case 'opto'
                    obj.TAUC_opto=TAUC;
                    obj.Tmean_opto=Tmean;
            end
            save(fileNameT,'TAUC','Tmean');
            save(fileNameT_inSummary,'TAUC','Tmean')
        end
        
        function [figRasterMean,figBehavior] = mPlotRasterPSTHoptoVSctrl(obj,ind_ROI,activity_type,behEventAlignPool, masklickPool, behEventSortPool,i_selectivity)
            %   using other existed function
            obj.metadata.behEventAlign=behEventAlignPool;
            obj.metadata.masklick=masklickPool;
            obj.metadata.i_selectivity=i_selectivity;
            if isempty(ind_ROI) || min(ind_ROI)>obj.metadata.nROI || max(ind_ROI)<1
                iStart=1;
                iEnd=obj.metadata.nROI;
            else
                iStart=max(1,min(ind_ROI));
                iEnd=min(obj.metadata.nROI,max(ind_ROI));
            end
            for i_align=1:length(behEventAlignPool)
                subplot(1,length(behEventAlignPool),i_align);
                hold on;
                %decide some global variable
                behEventAlign=behEventAlignPool{i_align};%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
                behEventSort=behEventSortPool{i_align};% string can be in{'first lick','reward','go cue'};
                masklick=masklickPool{i_align};
                
                selectivitystr={'stimuli','sensory difficulty','sensory','choice'};%sensory means grouping difficulties;
                trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
                trialTypeStr=selectivitystr{i_selectivity};
                %plotting settings
                if strcmp(behEventAlign,'stim onset')
                    frameNumTime=[0.5,4];%from 5s before align point to 5s after align point
                elseif strcmp(behEventAlign,'delay onset')
                    frameNumTime=[1,3.5];
                else
                    frameNumTime=[1,3];%from 5s before align point to 5s after align point
                end
                mkdir([obj.path.root,filesep,obj.name]);
                if strcmp(activity_type,'dff')
                    savename_fig_activity=[obj.path.root,filesep,obj.name,filesep,obj.name,'-cat_f_inclusion'];
                    [figDff] = fLabelROIsDffIncludedTrialRange(obj.dff,savename_fig_activity,obj.SavedCaTrials,obj.metadata.trial2include);
                    for iROI=iStart:iEnd
                        savename_fig=[obj.path.root,filesep,obj.name,filesep,obj.name,'ROI-',num2str(iROI),'-alignTo',behEventAlign,'-sort',behEventSort,'-dff-cmpOptoCtrl-rasterPSTH'];
                        disp(savename_fig);
                        IDstr=[obj.name,'ROI-',num2str(iROI)];
                        IDstr=strrep(IDstr,'_','\_');%下标变为转义字符的下划线
                        [figRasterMean,figBehavior] = fPlotRasterPSTHoptoVSctrl(obj.Data_extract,behEventAlign, masklick, behEventSort,...
                            obj.dff(iROI,:),IDstr,frameNumTime,obj.metadata.ind_1stFrame,1000/obj.metadata.frT,obj.metadata.ind_tr_1,i_selectivity);
                        saveas(figRasterMean,[savename_fig,'.png'],'png');
                    end
                elseif strcmp(activity_type,'spkr')
                    savename_fig_activity=[obj.path.root,filesep,obj.name,filesep,obj.name,'-spkr_inclusion'];
                    [figDff] = fLabelROIsDffIncludedTrialRange(obj.spkr,savename_fig_activity,obj.SavedCaTrials,obj.metadata.trial2include);
                    for iROI=iStart:iEnd
                        savename_fig=[obj.path.root,filesep,obj.name,filesep,obj.name,'ROI-',num2str(iROI),'-alignTo',behEventAlign,'-sort',behEventSort,'-spkr-cmpOptoCtrl-rasterPSTH'];
                        disp(savename_fig);
                        IDstr=[obj.name,'ROI-',num2str(iROI)];
                        IDstr=strrep(IDstr,'_','\_');%下标变为转义字符的下划线
                        [figRasterMean,figBehavior] = fPlotRasterPSTHoptoVSctrl(obj.Data_extract,behEventAlign, masklick, behEventSort,...
                            obj.spkr(iROI,:),IDstr,frameNumTime,obj.metadata.ind_1stFrame,1000/obj.metadata.frT,obj.metadata.ind_tr_1,i_selectivity);
                        saveas(figRasterMean,[savename_fig,'.png'],'png');
                    end
                end

            end
        end
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
