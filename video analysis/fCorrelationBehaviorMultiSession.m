function [ output_args ] = fCorrelationBehaviorMultiSession(filepath, bodyparts,coordinates,pSigTtest,treshold4likelihood,behEventAlign,saveNameID,yrange,titlestr,nonOverlappingBin )
%FCORRELATIONBEHAVIORMULTISESSION plot coordinates shift or likelihood for
%multiple session multiple body parts, all video files are stored in one
%folder
%   filepath='F:\video tracking\M2 imaging video';%or for SC data
%   bodyparts- n-by-m cells, indicate which body parts to be plot
%   coordinates- n-by-m cells(same as bodyparts), indicate which
%       coordinates/likelihood to be plot
%   pSigTtest- 0.05(default), threshold for moving ttest p value to be
%       significant
%   treshold4likelihood- 0.5(default), threshold of likelihood to discard
%       data
%   behEventAlign-cell array, align to which event(string can be in {'Stim onset',
%       'go cue','first lick','first left lick','first right lick',
%       'answer','reward','start'},
%   saveNameID-save name for figures and .mat files for the plot
%   yrange- same size as bodyparts, flexible to set ylim
%   titlestr- same size as bodyparts, flexible to set title



savepath=[filepath,filesep,'summary figures'];
summaryFile=[filepath,filesep,'imaging_video_data_summary.xlsx'];
[~,~,temp]=xlsread(summaryFile,1);
dataSummaryT=cell2table(temp(2:end,:));
dataSummaryT.Properties.VariableNames =temp(1,:);
nSession=size(dataSummaryT,1);
%choosing right data or exluded trash data
[~,datasource,~]=xlsread(summaryFile,2);%%%%%%%%%%%%%%%choosing which DLC model 

%% global settings
%plot part from the pool
nbodyparts=size(bodyparts,1)*size(bodyparts,2);
% indBodyPartsChosen=1:nbodyparts;
% bodyparts=bodyparts(indBodyPartsChosen);
% coordinates=coordinates(indBodyPartsChosen);

% treshold4likelihood=0.5;

i_selectivity=4;%*********variable**************
selectivitystr={'stimuli','difficulty','sensory','choice'};%sensory means grouping difficulties;
trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
outcomeType='combineCorErr';%'combineCorErr','divideCorErr',

binsize=1;%bin size for moving ttest
binstep=1;%bin step for moving ttest
% pSigTtest=0.01;%significance level for moving ttest

meanEachSession=cell(nSession,nbodyparts);%each cell represent one session mean
meanDiffEachSession=cell(nSession,nbodyparts);%each cell represent one session mean
fr=24/nonOverlappingBin;
for iSession=1:nSession%can be a loop
    led_file=[filepath,filesep,dataSummaryT.OLEDFileName{iSession},'.csv'];
    behdata=[filepath,filesep,dataSummaryT.session{iSession},'_Virables.mat'];
    
    %% from .beh file calculated a vector of trial start, for trial start from video use
    load(behdata);
    timeTrialStartBeh=cellfun(@fTrialStartStr2double,Data_extract.Time_trialStart);% time(s) relative to start of session day
    name_file_led=strsplit(led_file,filesep);
    name_file_led2=strsplit(name_file_led{end},'.');
    if ~exist([filepath,filesep,name_file_led2{1},'.mat'],'file')
        T=readtable(led_file);%different sessions may need mannually set parameters and save result separately
        frameTrialStartVideo=fTrialStartByOLEDrefArduino( T, 24 ,timeTrialStartBeh);
        save([filepath,filesep,name_file_led2{1},'.mat'],'frameTrialStartVideo');
    else
        load([filepath,filesep,name_file_led2{1},'.mat'])
    end
    frameTrialStartVideo=reshape(frameTrialStartVideo,1,[]);
    frameTrialStartVideo=floor(frameTrialStartVideo/nonOverlappingBin);
    % figure;
    % plot(frameTrialStartVideo);
    
    
    %% get behavior event from .beh
    ind_tr_1=1;%using data from trial 1
    [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( Data_extract, frameTrialStartVideo, 1000/fr ,ind_tr_1);%get behavior event time
    [trialType,rule,trialTypeStr] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix','left',outcomeType);%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
    
    for iBodyPart=1:nbodyparts    
        %align to behavior event
        behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
        if strcmp(behEventAlign{iBodyPart},'stim onset') || strcmp(behEventAlign{iBodyPart},'delay onset')
            masklick='yes';
        else
            masklick='no';
        end
        if strcmp(behEventAlign{iBodyPart},'stim onset')
            frameNumTime=[0.5,2];%plotting settings, from 0.5s before align point to 2s after align point
        elseif strcmp(behEventAlign{iBodyPart},'delay onset')
            frameNumTime=[1,1.4];
        else
            frameNumTime=[1,2];
        end
        frameNum=double(floor(frameNumTime*fr));
        % import DLC result
        col_datasource=cellfun(@(x) strcmp(bodyparts{iBodyPart},x),datasource(1,:));
        row_datasource=cellfun(@(x) strcmp(dataSummaryT.session{iSession},x),datasource(:,1));
        DLCiteration=datasource{row_datasource,col_datasource};
        switch DLCiteration
            case 'iteration-1'   
                file_trace=[filepath,filesep,'iteration1',filesep,dataSummaryT.DLCFileName1{iSession},'.csv'];
                if ~exist(file_trace,'file')%when no multiple files, just place the file at root path
                    file_trace=[filepath,filesep,dataSummaryT.DLCFileName1{iSession},'.csv'];
                end
                flagSetNAN=0;
            case 'iteration-2'
                file_trace=[filepath,filesep,'iteration2',filesep,dataSummaryT.DLCFileName2{iSession},'.csv'];
                flagSetNAN=0;
            case 'iteration-3'
                file_trace=[filepath,filesep,'iteration3',filesep,dataSummaryT.DLCFileName3{iSession},'.csv'];
                flagSetNAN=0;
            case 'iteration-4'
                file_trace=[filepath,filesep,'iteration4',filesep,dataSummaryT.DLCFileName4{iSession},'.csv'];
                flagSetNAN=0;
            case 'iteration-5'
                file_trace=[filepath,filesep,'iteration5',filesep,dataSummaryT.DLCFileName5{iSession},'.csv'];
                flagSetNAN=0;
            otherwise %no data, skip this loop
                file_trace=[filepath,filesep,dataSummaryT.DLCFileName1{iSession},'.csv'];%for  iteration-1~3
%                 file_trace=[filepath,filesep,'iteration4',filesep,dataSummaryT.DLCFileName4{iSession},'.csv'];%for  iteration-4+
     
                flagSetNAN=1;
        end
        name_file_trace=strsplit(file_trace,'.');
        if ~exist([name_file_trace{1},'.mat'],'file')
            [dcnum,dctxt,~] = xlsread(file_trace,1);%this process is time consuming, so not do very often
            save([name_file_trace{1},'.mat'],'dcnum','dctxt');
        else
            load([name_file_trace{1},'.mat'])
        end
        %% get coordinates of body parts, analogy to dff
        indcol1=cellfun(@(x) strcmp(bodyparts{iBodyPart},x),dctxt(2,:));
        indcol2=cellfun(@(x) strcmp(coordinates{iBodyPart},x),dctxt(3,:));
        indcol_likelihood=cellfun(@(x) strcmp('likelihood',x),dctxt(3,:));
        indcol_co=find(indcol1.*indcol2);
        indcol_like=find(indcol1.*indcol_likelihood);
        bodyco=dcnum(:,indcol_co);
        bodycoli=dcnum(:,indcol_like);
        bodyco(bodycoli<treshold4likelihood)=nan;%rule out those low likelihood data
        if ~strcmp(coordinates{iBodyPart},'likelihood')
            if contains(bodyparts{iBodyPart},'LickPort')
                bodyco=fBaselineCorrection(bodyco,5*fr);%40s as span
            end
            bodyco=bodyco-nanmean(bodyco);%calculate pixel shift, if it is likelihood, no need for normalization
%             % rule out those body parts with low likelihood
%             if mean(bodycoli)<treshold4likelihood
%                 bodyco(:)=nan;
%             end

        end
        %bin tongue data with 2 frames, non-overlapping
        bodyco=bodyco(1:floor(length(bodyco)/nonOverlappingBin)*nonOverlappingBin,1);
        bodyco=reshape(bodyco,nonOverlappingBin,[]);
        bodyco=nanmean(bodyco,1);
%         %replace nan as zeros
%         bodyco(isnan(bodyco))=0;
        %remove whole session if low likelihood
        if flagSetNAN==1
            bodyco(:)=nan;
        end
        %% aligned the coordinates change of body parts to behavior events
        if strcmp(behEventAlign{iBodyPart},'delay onset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
            disp(bodyparts{iBodyPart});
            [ bodyco_aligned, ~,~ ] = fAlignDelaySigal( bodyco, behEventFrameIndex,  frameNum ,'raw');
        else
            [ bodyco_aligned, ~, ~ ] = fAlignSigBehEvent( bodyco, behEventFrameIndex,lickingFrameIndex,behEventAlign{iBodyPart},frameNum );%decide which behavior event to align to
        end
        nfigcol=size(trialType,1)-2;%4 column(correct/error/miss/violation),companied with 4 lick raster
        meanEachSession{iSession,iBodyPart}=cell(size(trialType,2),nfigcol);
        meanDiffEachSession{iSession,iBodyPart}=cell(1,nfigcol);%diff between two groups mean
        for nStim=1:size(trialType,2) %for each stimulus
            for  nResult=1:nfigcol
                selectedTrialInd=trialType(nResult,nStim,:);
                selectedTrialInd=logical(squeeze(selectedTrialInd))';
                neuralActivity=bodyco_aligned(selectedTrialInd,:);
                [neuralActivityMean, ~]=fMean_CI(neuralActivity,0.05);
                meanEachSession{iSession,iBodyPart}{nStim,nResult}=neuralActivityMean;
            end
        end
        for  nResult=1:nfigcol
            tempmeanEachSession1=cell2mat(meanEachSession{iSession,iBodyPart}(1:size(trialType,2)/2,nResult));
            tempmeanEachSession2=cell2mat(meanEachSession{iSession,iBodyPart}(size(trialType,2)/2+1:end,nResult));
            meanDiffEachSession{iSession,iBodyPart}{1,nResult}=nanmean(tempmeanEachSession1,1)-nanmean(tempmeanEachSession2,1);
        end
    end
end
save([savepath,filesep,saveNameID],'meanEachSession','meanDiffEachSession');
%% settings for plot
%plot color plot,4 column(correct/error/miss/violation)
nfigcol=size(trialType,1)-2;%4 column(correct/error/miss/violation),companied with 4 lick raster
figMeanTrace=figure;%plot mean trace
% set(gcf, 'position', [0 0 250*nfigcol 200*nbodyparts]);
set(gcf, 'PaperPosition', [0 0 2*size(bodyparts,1) 2*size(bodyparts,2)]);
figMeanDiff=figure;%plot difference of mean trace for different choice etc. in one session
set(gcf, 'PaperPosition', [0 0 2*size(bodyparts,1) 2*size(bodyparts,2)]);
load([savepath,filesep,saveNameID]);
for iBodyPart=1:nbodyparts
    %% plot coordinates change, ref plotDffPSTH
%     if strcmp(outcomeType,'divideCorErr')
%         titlestr={'Correct','Error','Miss','Violation'};
%     else
%         titlestr={'Do','Miss','Violation'};
%     end
%     titlestr=strcat(bodyparts{iBodyPart},'-',coordinates{iBodyPart},'-',titlestr);
    if size(trialType,2)==6
        color_mean_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};
    elseif size(trialType,2)==2
        color_mean_trace={[0 0 1],[1 0 0]};
        color_cases={[0.5 0.5 1],[1 0.5 0.5]};
    end
    if strcmp(rule,'low click rate-right')
        color_mean_trace=fliplr(color_mean_trace);
    end
    for nResult=1:nfigcol
        matCasesStimResult=cell(size(meanEachSession,1),size(trialType,2));
        for  nStim=1:size(trialType,2) %for each stimulus
            for iSession=1:nSession %plot individual traces
                meanEachSessionCase=meanEachSession{iSession,iBodyPart}{nStim,nResult};
                matCasesStimResult{iSession,nStim}=meanEachSessionCase;
            end
        end
        %using t-test/ AUC to test its significance
        %pTtest = fMovingTtest(cell2mat(matCasesStimResult(:,1)),cell2mat(matCasesStimResult(:,2)),binsize,binstep,'ttest');
        n_shuffle=1000;    %%%%%%%%%%%%%%%%%%%%%%%%%%%change for shuffle times  
        category1=cell2mat(matCasesStimResult(:,1));%note here that only useful for 2 stimuli, if multiple stimuli, then this cause error
        category2=cell2mat(matCasesStimResult(:,2));
        activity=cat(1,category1,category2);
        label = cat(1,ones(size(category1,1),1),2*ones(size(category2,1),1));
        poslabel=2;
        fileAUC=strcat(saveNameID,bodyparts{iBodyPart},coordinates{iBodyPart},behEventAlign{iBodyPart},'binsize',num2str(binsize),'nshuffle',num2str(n_shuffle),'pAUC.mat');
        if exist([savepath,filesep,fileAUC],'file')
            load([savepath,filesep,fileAUC]);
        else
            [auc,pAUC,activitySmoothed] = fMovingAUC(label,activity,poslabel,n_shuffle,binsize,binstep);
        end
        save([savepath,filesep,fileAUC],'auc','pAUC','activitySmoothed');
        %plot mean, ci, individual stimuli etc after AUC analysis.
        for  nStim=1:size(trialType,2) %for each stimulus; in fact only support 2 category
            figure(figMeanTrace);%save mean trace
            subplot(size(bodyparts,2),size(bodyparts,1),iBodyPart);
            for iSession=1:nSession %plot individual traces
%                 curve_case(nStim)=plot(1:size(activitySmoothed{nStim},2),activitySmoothed{nStim}(iSession,:),'Color',color_cases{nStim},'linewidth',1);
                curve_case(nStim)=plot(1:size(matCasesStimResult{iSession,nStim},2),matCasesStimResult{iSession,nStim},'Color',color_cases{nStim},'linewidth',1);
                hold on;
            end
            %plot mean trace
            figure(figMeanTrace);%save mean trace
            subplot(size(bodyparts,2),size(bodyparts,1),iBodyPart);
%             [neuralActivityMean,neuralActivityCI]=fMean_CI(activitySmoothed{nStim},0.05);
            [neuralActivityMean,neuralActivityCI]=fMean_SE(activitySmoothed{nStim});
            activity_Std=nanstd(activitySmoothed{nStim},1,1);
            neuralActivityCI(1,:)=nanmean(activitySmoothed{nStim},1)+activity_Std;
            neuralActivityCI(2,:)=nanmean(activitySmoothed{nStim},1)-activity_Std;
%             %plot CI
%             if size(trialType,2)==2%if only two stimuli, plot CI, otherwise, do not plot for the reason of simplicity
%                 xpatch=[1:size(neuralActivityCI,2), fliplr(1:size(neuralActivityCI,2))];
%                 ypatch=[neuralActivityCI(1,:),fliplr(neuralActivityCI(2,:))];
%                 p=patch(xpatch,ypatch,color_mean_trace{nStim});%plot confidence interval
%                 p.FaceAlpha=0.1;
%                 p.EdgeColor=color_mean_trace{nStim};%'none';
%                 hold on;
%             end
            %plot mean
            curve_meanTrace(nStim)=plot(1:size(neuralActivityMean,2),neuralActivityMean,'Color',color_mean_trace{nStim},'linewidth',2);
            hold on;
            title(titlestr{iBodyPart});
        end
        if ~exist('yrange','var') || isempty(yrange)
            y_lim=get(gca,'Ylim');
        else
            y_lim=yrange{iBodyPart};
        end
%       %plot significant bar (indicating time point of significant)
        indSig=logical((pAUC<pSigTtest/2)+(pAUC>1-pSigTtest/2));
%         indSig=(pTtest<pSigTtest);
        ySig=ones(size(pAUC))*y_lim(2)*0.9;
        xSig=(1:length(pAUC));
        xSig(~indSig)=nan;
        ySig(~indSig)=nan;
        xSig=fRuleOutOccasional(xSig,1);
        ySig=fRuleOutOccasional(ySig,1);
        plot(xSig,ySig,'k-','LineWidth',1);
        
        %plot meandiff
        matMeanDiffEachSession=cell(size(meanDiffEachSession));
        for iSession=1:nSession
            matMeanDiffEachSession{iSession,iBodyPart}=meanDiffEachSession{iSession,iBodyPart}{1,nResult};
        end
        figure(figMeanDiff);
        subplot(size(bodyparts,2),size(bodyparts,1),iBodyPart);
        matCasesStimResult=cell2mat(matMeanDiffEachSession(:,iBodyPart));
        %calculate mean of the distribution comparing to zero by bootstrp
        nboot=n_shuffle;
        activitydiff=matCasesStimResult;
        fileBoot=strcat(saveNameID,bodyparts{iBodyPart},coordinates{iBodyPart},behEventAlign{iBodyPart},'binsize',num2str(binsize),'nboot',num2str(nboot),'pBoot.mat');
        if exist([savepath,filesep,fileBoot],'file')
            load([savepath,filesep,fileBoot]);
        else
            [ meanBoot, ciBoot, prctileBoot,activityBootBin ] = fMovingBootstrpMeanCIPrctile( activitydiff,nboot,binsize,binstep );
        end
        save([savepath,filesep,fileBoot],'meanBoot', 'ciBoot', 'prctileBoot','activityBootBin');
        activitydiff_Std=nanstd(activitydiff,1,1);
        ciBoot(1,:)=nanmean(activitydiff,1)+activitydiff_Std;
        ciBoot(2,:)=nanmean(activitydiff,1)-activitydiff_Std;
%         %plot individual traces (smoothed)
%         for iSession=1:nSession
%             figure(figMeanDiff);
%             subplot(size(bodyparts,2),size(bodyparts,1),iBodyPart);
%             curve_Diffcase(nStim)=plot(1:size(activityBootBin,2),activityBootBin(iSession,:),'Color',[0.5,0.5,0.5],'linewidth',1);
%             hold on;
%         end
        %plot CI/STD
        if size(trialType,2)==2%if only two stimuli, plot CI, otherwise, do not plot for the reason of simplicity
            xpatch=[1:size(ciBoot,2), fliplr(1:size(ciBoot,2))];
            ypatch=[ciBoot(1,:),fliplr(ciBoot(2,:))];
            p=patch(xpatch,ypatch,[0,0,0]);%plot confidence interval
            p.FaceAlpha=0.1;
            p.EdgeColor=[0,0,0];%'none';
            hold on;
        end
        %plot mean
        curve_meanDiff(nStim)=plot(1:size(nanmean(activitydiff,1),2),nanmean(activitydiff,1),'k-','linewidth',2);
        hold on;
        plot([1,size(activitydiff,2)],[0,0],'k-','linewidth',1);
        title(titlestr{iBodyPart});
        indSig=logical((prctileBoot<pSigTtest/2)+(prctileBoot>1-pSigTtest/2));
        ySig=ones(size(prctileBoot))*y_lim(2)*0.9;
        xSig=(1:length(prctileBoot));
        xSig(~indSig)=nan;
        ySig(~indSig)=nan;
        xSig=fRuleOutOccasional(xSig,1);
        ySig=fRuleOutOccasional(ySig,1);
        plot(xSig,ySig,'b-','LineWidth',1);
    end
    
end
figures=cell(1,2);
figures{1}=figMeanTrace;
figures{2}=figMeanDiff;
for ifig=1:2
    for iBodyPart=1:nbodyparts
        if strcmp(behEventAlign{iBodyPart},'stim onset')
            frameNumTime=[0.5,2];%plotting settings, from 0.5s before align point to 2s after align point
            tempXticklabel=[0,1,2];
        elseif strcmp(behEventAlign{iBodyPart},'delay onset')
            frameNumTime=[1,1.4];
            tempXticklabel=[-.5,0,.5,1];
        else
            frameNumTime=[1,2];
            tempXticklabel=[-1,0,1,2];
        end
        frameNum=double(floor(frameNumTime*fr));
        figure(figures{ifig});
        subplot(size(bodyparts,2),size(bodyparts,1),iBodyPart);
        if strcmp(coordinates{iBodyPart},'likelihood')
            yrange={[ 0 , 1.1 ]};
        end
        if ~exist('yrange','var')|| isempty(yrange)
            y_lim=get(gca,'Ylim');
        else
            y_lim=yrange{iBodyPart};
        end
        xlim([1,sum(frameNum)+1]);
        plot([frameNum(1)+1,frameNum(1)+1],[y_lim(1),y_lim(2)],'k-');%align to a behavior event
        hold on;
        if strcmp(behEventAlign{iBodyPart},'stim onset')
            plot((frameNum(1)+1+0.5*fr)*ones(1,2),[y_lim(1),y_lim(2)],'k-');%plot delay onset
        elseif strcmp(behEventAlign{iBodyPart},'delay onset')
            plot((frameNum(1)+1-0.5*fr)*ones(1,2),[y_lim(1),y_lim(2)],'k-');%plot stim onset
        end
        if ceil(iBodyPart/size(bodyparts,1))==size(bodyparts,2)
            xlabel(['Time (s) from ',behEventAlign{iBodyPart}]);
        end
        if mod(iBodyPart,size(bodyparts,1))==1 || size(bodyparts,1)==1
            if strcmp(coordinates{iBodyPart},'likelihood')
                ylabel('likelihood');
            else
                if ifig==1
                    ylabel('{\it\Delta}pixel');
                elseif ifig==2
                    ylabel('{\it\Delta}pixel (ipsi-contra)');
                end
            end
        end
        set(gca,'xtick',frameNum(1)+fr*tempXticklabel+1,'xticklabel',tempXticklabel);
%         set(gca,'xtick',[floor(fr*(frameNumTime(1)-floor(frameNumTime(1)))):floor(fr):size(neuralActivity,2)],'xticklabel',tempXticklabel);
        set(gca,'FontName','Arial','FontSize',12);
        set(gca,'Ylim',y_lim);
        box off;
    end
end

figure(figMeanTrace);
subplot(size(bodyparts,2),size(bodyparts,1),1);
if contains(trialTypeStr,'stimuli')
    h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','','contra hard','','contra easy'},'Location','best');
    %legend(curve_meanTrace(:),num2str(Data_extract.Stimuli(:)),'Location','best');
elseif contains(trialTypeStr,'difficulty')
    h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','contra hard','contra easy'},'Location','best');
elseif contains(trialTypeStr,'first lick')
    h=legend(curve_meanTrace(:),{'ipsi lick first','contra lick first'},'Location','best');
elseif contains(trialTypeStr,'sensory')
    h=legend(curve_meanTrace(:),{'low click rate','high click rate'},'Location','best');
end
set(h,'box','off');

saveas(figMeanTrace,[savepath,filesep,saveNameID,'-coordinates-likelihood-algin to ',behEventAlign{iBodyPart},'-sort ',behEventSort,'mean trace.pdf'],'pdf');
saveas(figMeanDiff,[savepath,filesep,saveNameID,'-Diff-coordinates-likelihood-p','01.pdf'],'pdf');

end

%% auxiliary fuction
function [timeTrialStart]=fTrialStartStr2double(strTrialStart)
%transform the trial start time from string form to double form
t=strsplit(strTrialStart,'_');
t=str2double(t);
timeTrialStart=t(4)*3600+t(5)*60+t(6)+t(7)/1000; %time(s) relative to trial1
end
function [out]=fBaselineCorrection(in,span)
%span = round(20*1000/frT/2); %every 40s, get a mean and substrate to compensate for baseline change
x = [];
for i = 1:length(in)
    %     ind_x = ind_x + 1;
    ind1 = max(i- span,1);
    ind2 = min(i+ span,length(in));
    x(i) = prctile(in(ind1:ind2),8);
end
x=reshape(x,[],1);
out = in - x + nanmean(x);
end
