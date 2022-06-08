function [TSVM,figMeanSigEpoch,figMeanSig,figMeanSigGo] = fGetEpochSVMscoreASession(filepath,animal,date,field,celltype,dataForm,indROI,strROI,indTrial2use,trialTypeStr,SVMtype,nRepeat, pTraining,CorrectedMethod,smooth_binsize,CTflag)
%FGETEPOCHSVMSCOREASESSION calculate epoch SVM from mulitiple sessions,
%with variables- ITI,sound,delay,response,lick, etc. only for control
%session, not opto session
%Input-filepath, path of one imaging session where beh, dff exist
%   CorrectedMethod-method to correct AUC value if correct and error 
%   trials are combined to calculate sensory/choice
%   AUC£¬{'balencedCorErrTrialNum'£¬'SensoryChoiceOrthogonalSubtraction','raw'},
%   indROI- choose subpopulation of ROI to perform calculation, for
%   futher compare different cell's function
%   strROI- a label to name indROI
%   CTflag- {'CT','TT','CE','EE',1,2,3,4,5}, 'CT' all cross time SVM, 'TT' only
%       moving SVM, 'CE' cross epoch, 'EE' moving epoch SVM, 
%       1-5, only ITI, sound, delay, respone, lick epoch; other data are
%       nan, but the variable exist to make table consistent in size.
%   ref fGetEpochAUCtableASession
%Output- TSVM with varibles{animal, date, field, ROI, ITI, sound,
%delay, response, lick, pITI, psound, pdelay, presponse, plick, celltype};
%   V22.5.26
cd(filepath);
dirmat=strcat(filepath,filesep,'*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
file_imaging=cellfun(@(x) contains(x,'CaTrialsSIM_'), filenames);
load(filenames{file_imaging});%load imaging data
file_beh=cellfun(@(x) contains(x,'Virables'), filenames);
load(filenames{file_beh});

session=[animal,'_',date,'_',field];
trial2include='all';
trial2exclude=[];
objsession=Session2P(session,filepath,trial2include,trial2exclude);
activity_data_raw=objsession.mSmoothFR(dataForm,smooth_binsize);%{objsession.zscored_spkr,objsession.dff,objsession.spkr}%choose among several acitivities form
if ~isempty(indROI)
    activity_data=activity_data_raw(indROI,:);%use part of the ROI to perform SVM
else
    activity_data=activity_data_raw;
end
indTrial2use=reshape(indTrial2use,[],1);
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
str_nFrames='1s';%'500ms';%'1s'
frameNumTime_delay=[1,1.2];%from 1s before delay onset to 1s after that, since 1-1.5
frameNumTime_go=[1,1.5];%from 5s before align point to 5s after align point
%calculate moving SVM
binsize=5;
binstep=5;
    
%naming the temp file by method of AUC calculation correction
filepath_inSummary='E:\2P\summary\SVM\cases';%backup another copy in the summary folder
fileName_datapars=[animal,'-',datestr,'-',strROI,'_',dataForm,'-trialType',trialTypeStr,...
    '-epochBin',str_nFrames,'-Tdelay',num2str(frameNumTime_delay(1)),'-',num2str(frameNumTime_delay(2)),...
    '_go',num2str(frameNumTime_go(1)),'-',num2str(frameNumTime_go(2))];
if strcmp(trialTypeStr,'cor') || strcmp(trialTypeStr,'err')
    combineCorErr='divideCorErr';%and only use correct trial
    corErrTrialNumber='raw';
    fileName_processpars=[SVMtype,'-movingTimeBin',binsize,'-',CTflag,'-pTraining',num2str(pTraining)];
elseif strcmp(trialTypeStr,'cor and err')%used to test whether it is sensory or choice AUC,in order to be more fair, keep trial numbers balenced for correct and error trials
    combineCorErr='combineCorErr';%{'combineCorErr','divideCorErr'}
    fileName_processpars=[SVMtype,'-',CorrectedMethod,'-movingTimeBin',binsize,'-',CTflag,'-pTraining',num2str(pTraining)];
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
fileNameT=[filepath,filesep,fileName_datapars,'-',fileName_processpars,'-EpochSVM.mat'];
fileNameT_inSummary=[filepath_inSummary,filesep,fileName_datapars,'-',fileName_processpars,'-EpochSVM.mat'];

if exist(fileNameT,'file')
    load(fileNameT);
    disp(['Table exist, use ',animal,date]);
elseif exist(fileNameT_inSummary,'file')
    load(fileNameT_inSummary);
    disp(['Table exist, use ',animal,date]);
else
    %
    disp(['analyzing',fileNameT]);
    nROI=1;%each session was viewed as one ROI
    varNROI=size(activity_data,1);
    [varindROI{1:nROI}]=deal(indROI);
    varindROI=reshape(varindROI,[],1);
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
    %time range for delay and go algnment
    frameNum_delay=double(round(frameNumTime_delay*1000/frT));
    ts_delay=-frameNumTime_delay(1):frT/1000:frameNumTime_delay(2);
    frameNum_go=double(round(frameNumTime_go*1000/frT));
    ts_go=-frameNumTime_go(1):frT/1000:frameNumTime_go(2);
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
                ind_trial=logical(ind_trial.*indTrial2use);
                label_SVM=[ones(trialNumEachType*2,1);ones(trialNumEachType*2,1)*2];
            case 'raw'
                [trialType,~,~] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix','left',combineCorErr);
                label_choice = fTrialType2Label(trialType,2);
                if contains(trialTypeStr,'cor') % 'cor'| 'cor and err'
                    ind_trial=logical(reshape(sum(trialType(1,:,:),2),[],1));%only correct trials
                elseif strcmp(trialTypeStr,'err')
                    ind_trial=logical(reshape(sum(trialType(2,:,:),2),[],1));%only correct trials
                end
                ind_trial=logical(ind_trial.*indTrial2use);
                label_SVM=label_choice(ind_trial);
                %used for orthogonal subtraction
                [trialType_orthogonal,~,~] = fGetTrialType( Data_extract,[],7-trialTypeVar(i_selectivity),'matrix','left',combineCorErr);
                label_choice_orthogonal = fTrialType2Label(trialType_orthogonal,2);
                label_SVM_orthogonal=label_choice_orthogonal(ind_trial);
        end

        sigbyEpoch=[];
        sigMoving=[];
        sigMovingGo=[];
        for roiNo = 1:size(activity_data,1)
            %data for epoch activity
            [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,activity_data(roiNo,:),frT,str_nFrames);
            if strcmp(CorrectedMethod,'SensoryChoiceOrthogonalSubtraction')
                T_SigbyEpoch=fOrthogonalSubtraction(T_SigbyEpoch,ind_trial,label_SVM,label_SVM_orthogonal);
            else
                T_SigbyEpoch=T_SigbyEpoch(ind_trial,:);
            end
            sigbyEpoch=cat(3,sigbyEpoch,table2array(T_SigbyEpoch));
            %data for moving activity
            [ activity_data_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( activity_data(roiNo,:), behEventFrameIndex,  frameNum_delay );
            if strcmp(CorrectedMethod,'SensoryChoiceOrthogonalSubtraction')
                activity_data_aligned_ortho_corrected=fOrthogonalSubtraction(activity_data_aligned,ind_trial,label_SVM,label_SVM_orthogonal);
            else
                activity_data_aligned_ortho_corrected=activity_data_aligned(ind_trial,:);
            end
            sigMoving=cat(3,sigMoving,activity_data_aligned_ortho_corrected);
            %data for moving activity around go cue
            behEventAlign='go cue';
            [ activity_data_aligned_go, ~, ~ ] = fAlignSigBehEvent( activity_data(roiNo,:), behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum_go );
            if strcmp(CorrectedMethod,'SensoryChoiceOrthogonalSubtraction')
                activity_data_aligned_ortho_corrected=fOrthogonalSubtraction(activity_data_aligned_go,ind_trial,label_SVM,label_SVM_orthogonal);
            else
                activity_data_aligned_ortho_corrected=activity_data_aligned_go(ind_trial,:);
            end
            sigMovingGo=cat(3,sigMovingGo,activity_data_aligned_ortho_corrected);
        end
        sigbyEpoch=permute(sigbyEpoch,[1,3,2]);%now, 1d trial number, 2d roi, 3d frames
        sigMoving=permute(sigMoving,[1,3,2]);
        sigMovingGo=permute(sigMovingGo,[1,3,2]);
        
        [CESVM.score, CESVM.score_shuffle,CESVM.p,CESVM.p_shuffled,CESVM.ts] = fCrossTimeSVM(sigbyEpoch,label_SVM, nRepeat, pTraining,1,1,1:5,CTflag);
        %calculate moving SVM
        if strcmp(CTflag,'CT') ||strcmp(CTflag,'TT')
            [delayCTSVM.score, delayCTSVM.score_shuffle, delayCTSVM.p, delayCTSVM.p_shuffled, delayCTSVM.ts]=fCrossTimeSVM(sigMoving,label_SVM, nRepeat, pTraining,binsize,binstep,ts_delay,CTflag);
            [goCTSVM.score, goCTSVM.score_shuffle, goCTSVM.p, goCTSVM.p_shuffled, goCTSVM.ts]=fCrossTimeSVM(sigMovingGo,label_SVM, nRepeat, pTraining,binsize,binstep,ts_go,CTflag);
        else
            [delayCTSVM.score, delayCTSVM.score_shuffle, delayCTSVM.p, delayCTSVM.p_shuffled, delayCTSVM.ts]=fCrossTimeSVM(sigMoving,label_SVM, nRepeat, pTraining,binsize,binstep,ts_delay,'nan');
            [goCTSVM.score, goCTSVM.score_shuffle, goCTSVM.p, goCTSVM.p_shuffled, goCTSVM.ts]=fCrossTimeSVM(sigMovingGo,label_SVM, nRepeat, pTraining,binsize,binstep,ts_go,'nan');
        end
    elseif strcmp(SVMtype,'stimuli') %here, compare auc of cor/err for each stimuli
        [ITI, sound,delay, response, lick, pITI, psound, pdelay, presponse, plick]= deal(nan(nROI,size(trialType,2)));
    end
    
    TSVM=table(varanimal,vardate,varfield,varcelltype,varNROI,varindROI,CESVM,delayCTSVM,goCTSVM,...
        'VariableNames',{'animal','date','field','celltype','nROI','indROI','CESVM','delayCTSVM','goCTSVM'});
%}
end
save(fileNameT,'TSVM');
save(fileNameT_inSummary,'TSVM');%save another copy in the summary folder

%plot cross time 2D_SVM results
pSig2show=0.05;
%plot epoch SVM
[figMeanSigEpoch] = plotCrossTimeSVM(TSVM.CESVM,pSig2show,SVMtype,{'I','S','D','R','L'});
saveas(figMeanSigEpoch,[filepath,filesep,animal,'-',datestr,'trialType',trialTypeStr,'-',SVMtype,'-',CorrectedMethod,'-timeBin',str_nFrames,'-pTraining',num2str(pTraining),'-EpochSVM.pdf'],'pdf');
saveas(figMeanSigEpoch,[filepath,filesep,animal,'-',datestr,'trialType',trialTypeStr,'-',SVMtype,'-',CorrectedMethod,'-timeBin',str_nFrames,'-pTraining',num2str(pTraining),'-EpochSVM.png'],'png');
%plot time moving SVM
[figMeanSig] = plotCrossTimeSVM(TSVM.delayCTSVM,pSig2show,SVMtype);
[figMeanSigGo] = plotCrossTimeSVM(TSVM.goCTSVM,pSig2show,SVMtype);
%save plots
saveas(figMeanSig,[filepath,filesep,animal,'-',datestr,'trialType',trialTypeStr,'-',SVMtype,'-',CorrectedMethod,'-timeBin',str_nFrames,'-pTraining',num2str(pTraining),'-DelaySVM.pdf'],'pdf');
saveas(figMeanSig,[filepath,filesep,animal,'-',datestr,'trialType',trialTypeStr,'-',SVMtype,'-',CorrectedMethod,'-timeBin',str_nFrames,'-pTraining',num2str(pTraining),'-DelaySVM.png'],'png');
saveas(figMeanSigGo,[filepath,filesep,animal,'-',datestr,'trialType',trialTypeStr,'-',SVMtype,'-',CorrectedMethod,'-timeBin',str_nFrames,'-pTraining',num2str(pTraining),'-GoSVM.pdf'],'pdf');
saveas(figMeanSigGo,[filepath,filesep,animal,'-',datestr,'trialType',trialTypeStr,'-',SVMtype,'-',CorrectedMethod,'-timeBin',str_nFrames,'-pTraining',num2str(pTraining),'-GoSVM.png'],'png');
    
end

function [figMeanSig] = plotCrossTimeSVM(CTSVM,pSig,SVMtype,varargin)
%plot accuracy of each session
score=CTSVM.score;
score_shuffle=CTSVM.score_shuffle;
pCT=CTSVM.p;
pCT_shuffle=CTSVM.p_shuffled;
[score_tt,score_shuffle_tt]=deal(zeros(size(score,1),size(score,2)));
for i=1:size(score,1)
    score_tt(i,:)=diag(squeeze(score(i,:,:)));
    score_shuffle_tt(i,:)=diag(squeeze(score_shuffle(i,:,:)));
end
pTT=diag(pCT);
pTT_shuffle=diag(pCT_shuffle);

figMeanSig=figure;
set(gcf,'position',[200,200,700,200]);
colorScore={[1,0,0],[0,0,1],[0.5,0.5,0.5]};

subplot(1,3,1);
if strcmp(SVMtype,'choice')
    fPlotMean_CI(CTSVM.ts,score_tt,colorScore{1},pSig);
elseif strcmp(SVMtype,'sensory')
    fPlotMean_CI(CTSVM.ts,score_tt,colorScore{2},pSig);
end
fPlotMean_CI(CTSVM.ts,score_shuffle_tt,colorScore{3},pSig);
%label significance
indSig=pTT<0.05;
markersize=fSigMarkerSize(pTT);
y_lim=get(gca,'Ylim');
ySig=ones(size(pTT))*y_lim(2)*0.9;
xSig=CTSVM.ts(1:length(pTT));
xSig(~indSig)=[];
ySig(~indSig)=[];
markersize(~indSig)=[];
if strcmp(SVMtype,'choice')
    scatter(xSig,ySig,'Marker','*','SizeData', markersize,'MarkerEdgeColor','r');
elseif strcmp(SVMtype,'sensory')
    scatter(xSig,ySig,'Marker','*','SizeData', markersize,'MarkerEdgeColor','b');
end
hold on;
%color plot of the cross temporal decoding
ax2=subplot(1,3,2);
fPlotMeanSig2D_SVM(ax2, score,pCT,CTSVM.ts,pSig,'no bar');
ylabel('Training time');
xlabel('Testing time');
title('Real data');
ax3=subplot(1,3,3);
fPlotMeanSig2D_SVM(ax3, score_shuffle,pCT_shuffle, CTSVM.ts,pSig,'colorbar');
ylabel('Training time');
xlabel('Testing time');
title('Shuffled data');
%for epoch SVM, label epochs
if ~isempty(varargin)
    ts_label=varargin{1};
    subplot(1,3,1);
    set(gca,'XTick',1:length(ts_label),'XTickLabel',ts_label);
    for i=2:3
        subplot(1,3,i);
        set(gca,'XTick',1:length(ts_label),'XTickLabel',ts_label);
        set(gca,'YTick',1:length(ts_label),'YTickLabel',ts_label);
    end
end
end

function markersize=fSigMarkerSize(p)
markersize(p>0.05)=0;
markersize(logical((p<=0.05).*(p>0.01)))=4;
markersize(logical((p<=0.01).*(p>0.001)))=8;
markersize(logical(p<=0.001))=12;
end

