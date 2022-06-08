close all;
[num,txt,raw] =xlsread('C:\Users\PYX\Documents\DataSummary\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
savepath='E:\2P\summary';%'H:\2P\summary\summary_DLCfiltered';%
clear TSVM_combine;
% trialTypeStr={'cor'};%'cor and err';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
% SVMtype={'choice'};%'stimuli' means comparing cor and err for each stimuli
% AUCCorrectedMethod='None';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
trialTypeStr={'cor and err'};%'cor and err';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
SVMtype={'sensory','choice'};%'stimuli' means comparing cor and err for each stimuli
AUCCorrectedMethod='SensoryChoiceOrthogonalSubtraction';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
nRepeat=100; %shuffled times
pTraining=0.9;
%total loop
activity_form={'dff','spkr'};
activity_smooth_binsize={30,30};%ms, 200for dff, 400 for spkr
% % 1-part
% activity_form={'dff'};
% activity_smooth_binsize={200};%ms, 200for dff, 400 for spkr
% %2-part
% activity_form={'spkr'};
% activity_smooth_binsize={400};%ms, 200for dff, 400 for spkr
strROI='allROI';
CTflag='CT';%{'CT','TT',1,2,3,4,5};
%% calculate SVM accuracy
%
for i_form=1:length(activity_form)
    dataForm = activity_form{i_form};
    binsize=activity_smooth_binsize{i_form};
    T=cell2table(raw(2:end,1:15));
    T.Properties.VariableNames=strrep(raw(1,1:15),' ','_');%table variable name can't have ' ',so replace them
    ind_session=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*contains(T.ROI_type,'soma').*(~contains(T.behavior_performance,'probe')).*(~contains(T.field,'soma'));
    ind_session=find(ind_session);
    n_session=length(ind_session);
    animal_unique=unique(T.animal(ind_session));
    n_animal=length(animal_unique);
    score=cell(n_session,1);
    session_name=cell(n_session*length(trialTypeStr)*length(SVMtype),1);
    celltype=cell(n_session*length(trialTypeStr)*length(SVMtype),1);
    i_fig=1;
    for i_trialtype=1:length(trialTypeStr)
        for i_SVMtype=1:length(SVMtype)
            clear TSVM_combine;
            for i_session=1:n_session
                indrow=ind_session(i_session);
                session_name{i_fig} =strcat(T.animal{indrow},'_',T.date{indrow});
                celltype{i_fig} =T.cell_type{indrow};
                trial2include='all';
                trial2exclude=[];
                objSession=Session2P(session_name{i_fig},T.file_path{indrow},trial2include,trial2exclude);
                nROI_total=objSession.metadata.nROI;
                indTrial2use=objSession.metadata.indTrial2use;
                indROI=[];
%                 [TSVM_currentSession,figMeanSigEpoch(i_fig),figMeanSig(i_fig),figMeanSigGo(i_fig)]= ...
%                     fGetEpochSVMscoreASession(T.file_path{indrow},T.animal{indrow},T.date{indrow},T.field{indrow},T.cell_type{indrow},dataForm,...
%                     indROI,strROI,indTrial2use,trialTypeStr{i_trialtype},SVMtype{i_SVMtype},nRepeat, pTraining,AUCCorrectedMethod,binsize,CTflag);
                %using class member function
                [TSVM_currentSession,objSession]=objSession.mSVMscore(T.cell_type{indrow},dataForm,indROI,strROI,trialTypeStr{i_trialtype},...
                    SVMtype{i_SVMtype},nRepeat, pTraining,AUCCorrectedMethod,binsize,CTflag)
                %plot cross time 2D_SVM results
                pSig2show=0.05;
                figMeanSigEpoch(i_fig) = objSession.mPlotCrossTimeSVM('CESVM',pSig2show,SVMtype{i_SVMtype},{'I','S','D','R','L'});
                figMeanSig(i_fig) = objSession.mPlotCrossTimeSVM('delayCTSVM',pSig2show,SVMtype{i_SVMtype});
                figMeanSigGo(i_fig) = objSession.mPlotCrossTimeSVM('goCTSVM',pSig2show,SVMtype{i_SVMtype});
                
                if exist('TSVM_combine','var')
                    TSVM_combine=vertcat(TSVM_combine,TSVM_currentSession);
                else
                    TSVM_combine=TSVM_currentSession;
                end
                i_fig=i_fig+1;
            end
            save([savepath,filesep,'trialType',trialTypeStr{i_trialtype},'-',SVMtype{i_SVMtype},'pTraining',num2str(pTraining),'-TepochSVM.mat'],'TSVM_combine');
        end
    end
    %}
    %%
    %save individual cases in a PPT file
    pptFileName=['E:\2P\summary\SVM\SVM_cases_',activity_form{i_form},'_smooth_binsize',num2str(activity_smooth_binsize{i_form}),'_trial-',trialTypeStr{i_trialtype},'_decode-',SVMtype{i_SVMtype},'_pTraining',num2str(pTraining),'.pptx'];
    isOpen  = exportToPPTX();
    if ~isempty(isOpen)
        % If PowerPoint already started, then close first and then open a new one
        exportToPPTX('close');
    end
    slidesize=[16 9];
    if exist(pptFileName,'file')
        exportToPPTX('open',pptFileName);
    else
        exportToPPTX('new','Dimensions',slidesize, ...
            'Title','Population_activities_cases', ...
            'Author','PYX', ...
            'Comments','This file has been automatically generated by exportToPPTX');
    end
    
    for isession=1:length(session_name)
        slideId = exportToPPTX('addslide');
        exportToPPTX('addpicture',figMeanSigEpoch(isession),'Position',[1 1 7 2]);
        exportToPPTX('addpicture',figMeanSig(isession),'Position',[8 1 7 2]);
        exportToPPTX('addpicture',figMeanSigGo(isession),'Position',[8 4 7 2]);
        exportToPPTX('addtext',session_name{isession},'Position',[1 0 3 0.5],'Vert','top');
        exportToPPTX('addtext',celltype{isession},'Position',[slidesize(1)-1 0 3 0.5],'Vert','top');
    end
    
    % Save and close (in one command)
    newFile = exportToPPTX('saveandclose',pptFileName);
    fprintf('New file has been saved: <a href="matlab:open(''%s'')">%s</a>\n',newFile,newFile);
    close all;
    
end
%}

%% plot SVM accuracy cases, first load data
celltype='vgat';
trialTypeStr='cor';
SVMtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
load([savepath,filesep,'trialType',trialTypeStr,'-',SVMtype,'pTraining',num2str(pTraining),'-TepochSVM.mat']);
ind_session=logical((strcmp(TSVM_combine.celltype,celltype)+strcmp(TSVM_combine.celltype,[celltype,'-flpo'])).*(TSVM_combine.nROI>10)...
    .*(~contains(TSVM_combine.field,'spine')).*(~contains(TSVM_combine.field,'field')).*(~contains(TSVM_combine.field,'dendrite')));
SVMsensory=TSVM_combine(ind_session,:);

SVMtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
load([savepath,filesep,'trialType',trialTypeStr,'-',SVMtype,'pTraining',num2str(pTraining),'-TepochSVM.mat']);
ind_session=logical((strcmp(TSVM_combine.celltype,celltype)+strcmp(TSVM_combine.celltype,[celltype,'-flpo'])).*(TSVM_combine.nROI>10)...
    .*(~contains(TSVM_combine.field,'spine')).*(~contains(TSVM_combine.field,'field')).*(~contains(TSVM_combine.field,'dendrite')));
SVMchoice=TSVM_combine(ind_session,:);

%% plot SVM accuracy of each session and SVM accuracy overall
%
colorScore={[1,0,0],[0,0,1],[0.5,0.5,0.5]};%choice,sensory,shuffle
pSig=0.05;
nroistr=num2cell(SVMchoice.nROI);
textstr=cellfun(@(x) [num2str(x),'cells'],nroistr,'UniformOutput',false);
y_lim=[0.4,1];
% plot by epochs
[figScoreCase,figMean]=fScoreCaseMean(SVMchoice,SVMsensory,'delay',pSig,textstr);
% plot moving accuracy during delay
[figScoreCase,figMean]=fScoreCaseMean(SVMchoice,SVMsensory,'delay',pSig,textstr);
% plot moving accuracy around go
[figScoreCase,figMean]=fScoreCaseMean(SVMchoice,SVMsensory,'go',pSig,textstr);
%}
%% plot cross temporal decoding accuracy of SVM around go cue
%{
figMeanSig=figure;
set(gcf,'position',[200,200,600,250]);
pSig=0.05;
%when calculate moving SVM
binsize=3;
binstep=3;
frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
frT=30;%SavedCaTrials.FrameTime;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%need further tune
ts_raw=-frameNumTime(1):frT/1000:frameNumTime(2);
tempgoSVM=SVMchoice.goMovingSVM;
n_session=length(tempgoSVM);
clear score score_shuffle;
for i_session=2%1:n_session 
    if strcmp(trialTypeStr,'cor and err')
        if exist('score','var')
            score=vertcat(score,tempgoSVM(i_session).doCT);
            score_shuffle=vertcat(score_shuffle,tempgoSVM(i_session).shuffle_doCT);
        else
            score=tempgoSVM(i_session).doCT;
            score_shuffle=tempgoSVM(i_session).shuffle_doCT;
        end
    elseif strcmp(trialTypeStr,'cor')
        if exist('score','var')
            score=vertcat(score,tempgoSVM(i_session).corCT);
            score_shuffle=vertcat(score_shuffle,tempgoSVM(i_session).shuffle_corCT);
        else
            score=tempgoSVM(i_session).corCT;
            score_shuffle=tempgoSVM(i_session).shuffle_corCT;
        end
    end

end
nframe=size(score,3);
ts=ts_raw(ceil(binsize/2):binsize:end);%ts also need to be binned
ts=ts(1:nframe);

ax1=subplot(1,2,1);%real data
[ meandata, sigmat,pmat] = fPlotMeanSig2D_SVM(ax1, score,ts,pSig, 0.5,'no bar');
xlabel('Training time');
ylabel('Testing time');
title('Real data');
ax2=subplot(1,2,2);%shuffled data
[ meandata, sigmat,pmat] = fPlotMeanSig2D_SVM(ax2, score_shuffle,ts,pSig, 0.5,'colorbar');
xlabel('Training time');
ylabel('Testing time');
title('Shuffled data');
%}
%% helper function

function markersize=fSigMarkerSize(p)
markersize(p>0.05)=0;
markersize(logical((p<=0.05).*(p>0.01)))=4;
markersize(logical((p<=0.01).*(p>0.001)))=8;
markersize(logical(p<=0.001))=12;
end

function [figScoreCase,figMean]=fScoreCaseMean(SVMchoice,SVMsensory,epochstr,pSig,textstr)
%plot SVM accuracy cases
%Input- scoredataC struct array with fields:score,score_shuffle,p,p_shuffled,ts
%
%
figScoreCase=figure;
n_col=4;
n_row=ceil(size(SVMchoice,1)/n_col);
set(gcf,'Position',[100,100,200*n_col,200*n_row]);
colorScore={[1,0,0],[0,0,1],[0.5,0.5,0.5]};%choice,sensory,shuffle
if strcmp(epochstr,'epoch')
    %x_tick_label={'ITI','sound','delay','response','lick'};
    x_tick_label={'I','S','D','R','L'};
    x_tick=1:5;
    [choice_score_tt,choice_score_shuffle_tt,sensory_score_tt,sensory_score_shuffle_tt]=deal(zeros(size(SVMchoice,1),5));
    scoreC=SVMchoice.CESVM;
    scoreS=SVMsensory.CESVM;
    ts=scoreC(1).ts;
    xlabelstr='Behavioral epoch';
else
    if strcmp(epochstr,'delay')
        scoreC=SVMchoice.delayCTSVM;
        scoreS=SVMsensory.delayCTSVM;
        xlabelstr='Time (s) from delay onset';
    elseif strcmp(epochstr,'go')
        scoreC=SVMchoice.goCTSVM;
        scoreS=SVMsensory.goCTSVM;
        xlabelstr='Time (s) from go cue';
    end
    ts=scoreC(1).ts;
    [choice_score_tt,choice_score_shuffle_tt,sensory_score_tt,sensory_score_shuffle_tt]=deal(zeros(size(SVMchoice,1),length(ts)));
    x_tick_label={num2str(ceil(ts(1))),'0',num2str(floor(ts(end)))};
    x_tick=[ceil(ts(1)),0,floor(ts(end))];
end

for i=1:size(SVMchoice,1)
    ax=subplot(n_row,n_col,i);
    [choice_score_tt(i,:),choice_score_shuffle_tt(i,:)]=plotMovingSVM(ax,scoreC,pSig,'choice');
    [sensory_score_tt(i,:),sensory_score_shuffle_tt(i,:)]=plotMovingSVM(ax,scoreS,pSig,'sensory');
    plot([ts(1),ts(end)],[0.5,0.5],'k--');
    title([SVMchoice.animal{i},'-',SVMchoice.date{i}]);
    text(0.1,0.9,[num2str(SVMchoice.nROI(i)),'cells']);
    set(gca,'Ylim',[0.3,1],'Xlim',[0.9,5.1],'XTick',x_tick,'XTickLabel',x_tick_label);
    ylabel('Decoding accuracy');
    xlabel(xlabelstr);
    set(gca,'FontSize',12);
%     h=gca;
%     h.XTickLabelRotation=45;
    box off;
    text(0.1,0.9,textstr{i},'Units','normalized');
end
set(gcf,'PaperPosition',[1,1,2,2]);
saveas(figScoreCase,[savepath,filesep,celltype,'-trialType',trialTypeStr,'-',SVMtype,'-',epochstr,'-pTraining',num2str(pTraining),'-TepochSVM_score_cases.pdf'],'pdf');
% plot SVM accuracy overall
figMean=figure;
set(gcf,'Position',[100,100,200,200]);
fPlotMean_SE(1:size(sensory_score_shuffle_tt,2),sensory_score_shuffle_tt,colorScore{3},'cases');
fPlotMean_SE(1:size(choice_score_shuffle_tt,2),choice_score_shuffle_tt,colorScore{3},'cases');
fPlotMean_SE(1:size(sensory_score_tt,2),sensory_score_tt,colorScore{2},'only cases');
fPlotMean_SE(1:size(choice_score_tt,2),choice_score_tt,colorScore{1},'only cases');
fPlotMean_SE(1:size(sensory_score_tt,2),sensory_score_tt,colorScore{2});
fPlotMean_SE(1:size(choice_score_tt,2),choice_score_tt,colorScore{1});
plot([ts(1),ts(end)],[0.5,0.5],'k--');
set(gca,'Ylim',[0.4,1],'Xlim',[0.9,5.1],'XTick',x_tick,'XTickLabel',x_tick_label);
ylabel('Decoding accuracy');
xlabel(xlabelstr);
set(gca,'FontSize',12);
% h=gca;
% h.XTickLabelRotation=45;
box off;
nboot=1000;
binsize=1;
binstep=1;
pSigTtest=0.01;
baseline=0.5;
ts=1:5;
[ meanBoot, ciBoot, prctileBoot,activityBootBin ] = fMovingBootstrpMeanCIPrctile( choice_score_tt,nboot,binsize,binstep,baseline);
pMovingChoice=(0.5-abs(0.5-prctileBoot))*2;
indSig=pMovingChoice<0.05;
markersize=fSigMarkerSize(pMovingChoice);
y_lim=get(gca,'Ylim');
ySig=ones(size(prctileBoot))*y_lim(2)*0.9;
xSig=ts(1:length(prctileBoot));
xSig(~indSig)=[];
ySig(~indSig)=[];
markersize(~indSig)=[];
scatter(xSig,ySig,'Marker','*','SizeData', markersize,'MarkerEdgeColor','r');
hold on;
[ meanBoot, ciBoot, prctileBoot,activityBootBin ] = fMovingBootstrpMeanCIPrctile( sensory_score_tt,nboot,binsize,binstep,baseline);
pMovingSensory=(0.5-abs(0.5-prctileBoot))*2;
indSig=pMovingSensory<0.05;
markersize=fSigMarkerSize(pMovingSensory);
y_lim=get(gca,'Ylim');
ySig=ones(size(prctileBoot))*y_lim(2)*0.95;
xSig=ts(1:length(prctileBoot));
xSig(~indSig)=[];
ySig(~indSig)=[];
markersize(~indSig)=[];
scatter(xSig,ySig,'Marker','*','SizeData', markersize,'MarkerEdgeColor','b');

set(gcf,'PaperPosition',[1,1,2,2]);
saveas(figMean,[savepath,filesep,celltype,'trialType',trialTypeStr,'-',SVMtype,'-',epochstr,'-pTraining',num2str(pTraining),'-TepochSVM_scoreMean.pdf'])

end



function [mean_score_tt,mean_score_shuffle_tt]=plotMovingSVM(ax,CTSVM,pSig,SVMtype,varargin)
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
mean_score_tt=mean(score_tt);
mean_score_shuffle_tt=mean(score_shuffle_tt);
axes(ax);
colorScore={[1,0,0],[0,1,0],[0.5,0.5,0.5]};
if strcmp(SVMtype,'choice')
    fPlotMean_CI(1:size(score_tt,2),score_tt,colorScore{1},pSig);
elseif strcmp(SVMtype,'sensory')
    fPlotMean_CI(1:size(score_tt,2),score_tt,colorScore{2},pSig);
end
fPlotMean_CI(1:size(score_tt,2),score_shuffle_tt,colorScore{3},pSig);
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
%for epoch SVM, label epochs
if ~isempty(varargin)
    ts_label=varargin{1};
    set(gca,'XTick',1:length(ts_label),'XTickLabel',ts_label);
end
end