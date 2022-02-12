close all;
[num,txt,raw] =xlsread('C:\Users\PYX\Documents\DataSummary\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
savepath='E:\2P\summary';%'H:\2P\summary\summary_DLCfiltered';%
clear TSVM_combine;
% trialTypeStr={'cor'};%'cor and err';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
% SVMtype={'choice'};%'stimuli' means comparing cor and err for each stimuli
% AUCCorrectedMethod='None';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
trialTypeStr={'cor and err'};%'cor and err';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
SVMtype={'sensory'};%'stimuli' means comparing cor and err for each stimuli
AUCCorrectedMethod='SensoryChoiceOrthogonalSubtraction';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
nRepeat=100;
pTraining=0.9;
dataForm = 'spkr';
%% calculate SVM accuracy
%
T=cell2table(raw(2:end,1:15));
T.Properties.VariableNames=strrep(raw(1,1:15),' ','_');%table variable name can't have ' ',so replace them
ind_session=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*contains(T.ROI_type,'soma');
ind_session=find(ind_session);
n_session=length(ind_session);
animal_unique=unique(T.animal(ind_session));
n_animal=length(animal_unique);
score=cell(n_session,1);
for i_trialtype=1:length(trialTypeStr)
    for i_SVMtype=1:length(SVMtype)
        clear TSVM_combine;
        for i_session=1:n_session
            indrow=ind_session(i_session);
            TSVM_currentSession= fGetEpochSVMscoreASession(T.file_path{indrow},T.animal{indrow},T.date{indrow},T.field{indrow},T.cell_type{indrow},trialTypeStr{i_trialtype},SVMtype{i_SVMtype},nRepeat, pTraining,AUCCorrectedMethod,dataForm);
            if exist('TSVM_combine','var')
                TSVM_combine=vertcat(TSVM_combine,TSVM_currentSession);
            else
                TSVM_combine=TSVM_currentSession;
            end
        end
        save([savepath,filesep,'trialType',trialTypeStr{i_trialtype},'-',SVMtype{i_SVMtype},'pTraining',num2str(pTraining),'-TepochSVM.mat'],'TSVM_combine');
    end
end
%}
%% plot SVM accuracy cases, first load data
celltype='vgat-flpo';
trialTypeStr='cor and err';
SVMtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
load([savepath,filesep,'trialType',trialTypeStr,'-',SVMtype,'pTraining',num2str(pTraining),'-TepochSVM.mat']);
ind_session=logical(strcmp(TSVM_combine.celltype,celltype).*(TSVM_combine.nROI>10)...
    .*(~contains(TSVM_combine.field,'spine')).*(~contains(TSVM_combine.field,'field')).*(~contains(TSVM_combine.field,'dendrite')));
SVMsensory=TSVM_combine(ind_session,:);

SVMtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
load([savepath,filesep,'trialType',trialTypeStr,'-',SVMtype,'pTraining',num2str(pTraining),'-TepochSVM.mat']);
ind_session=logical(strcmp(TSVM_combine.celltype,celltype).*(TSVM_combine.nROI>10)...
    .*(~contains(TSVM_combine.field,'spine')).*(~contains(TSVM_combine.field,'field')).*(~contains(TSVM_combine.field,'dendrite')));
SVMchoice=TSVM_combine(ind_session,:);

%% plot accuracy of each session
%{
colorScore={[1,0,0],[0,0,1],[0.5,0.5,0.5]};%choice,sensory,shuffle
pSig=0.05;
figScoreCase=figure;
n_col=4;
n_row=ceil(size(SVMchoice,1)/n_col);
set(gcf,'Position',[100,100,200*n_col,200*n_row]);
for i=1:size(SVMchoice,1)
    subplot(n_row,n_col,i);
    fPlotMean_SE(1:size(SVMchoice.shuffleEpochSVMscore{i},2),SVMchoice.shuffleEpochSVMscore{i},colorScore{3});
    fPlotMean_SE(1:size(SVMsensory.shuffleEpochSVMscore{i},2),SVMsensory.shuffleEpochSVMscore{i},colorScore{3});
    fPlotMean_SE(1:size(SVMchoice(i,6:10),2),cell2mat(table2array(SVMchoice(i,6:10))),colorScore{1});
    fPlotMean_SE(1:size(SVMsensory(i,6:10),2),cell2mat(table2array(SVMsensory(i,6:10))),colorScore{2}); 
    plot([1,5],[0.5,0.5],'k--');
    title([SVMchoice.animal{i},'-',SVMchoice.date{i}]);
    text(0.1,0.9,[num2str(SVMchoice.nROI(i)),'cells']);
    set(gca,'Ylim',[0.3,1],'Xlim',[0.9,5.1],'XTick',1:5,'XTickLabel',{'ITI','sound','delay','response','lick'});
    set(gca,'FontSize',12);
    h=gca;
    h.XTickLabelRotation=45;
    box off;
end
set(gcf,'PaperPosition',[1,1,2,2]);
saveas(figScoreCase,[savepath,filesep,celltype,'-trialType',trialTypeStr,'-',SVMtype,'pTraining',num2str(pTraining),'-TepochSVM_score',SVMchoice.animal{i},'-',strrep(SVMchoice.date{i},'/',''),'.pdf']);
%}
%% plot SVM accuracy overall
%{
figMean=figure;
set(gcf,'Position',[100,100,200,200]);
temp=table2array(SVMsensory(:,6:10));
SVMsensoryMean=cellfun(@nanmean,temp);
fPlotMean_SE(1:size(SVMsensoryMean,2),SVMsensoryMean,colorScore{2},'only cases');
temp=table2array(SVMchoice(:,6:10));
SVMchoiceMean=cellfun(@nanmean,temp);
fPlotMean_SE(1:size(SVMchoiceMean,2),SVMchoiceMean,colorScore{1},'only cases');
fPlotMean_SE(1:size(SVMsensoryMean,2),SVMsensoryMean,colorScore{2});
fPlotMean_SE(1:size(SVMchoiceMean,2),SVMchoiceMean,colorScore{1});
plot([1,5],[0.5,0.5],'k--');
set(gca,'Ylim',[0.4,1],'Xlim',[0.9,5.1],'XTick',1:5,'XTickLabel',{'ITI','sound','delay','response','lick'});
ylabel('Decoding accuracy');
set(gca,'FontSize',12);
h=gca;
h.XTickLabelRotation=45;
box off;
nboot=1000;
binsize=1;
binstep=1;
pSigTtest=0.01;
baseline=0.5;
ts=1:5;
[ meanBoot, ciBoot, prctileBoot,activityBootBin ] = fMovingBootstrpMeanCIPrctile( SVMchoiceMean,nboot,binsize,binstep,baseline);
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
[ meanBoot, ciBoot, prctileBoot,activityBootBin ] = fMovingBootstrpMeanCIPrctile( SVMsensoryMean,nboot,binsize,binstep,baseline);
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

% set(gcf,'PaperPosition',[1,1,2,2]);
% saveas(figMean,[savepath,filesep,celltype,'trialType',trialTypeStr,'-',SVMtype,'pTraining',num2str(pTraining),'-TepochSVM_scoreMean.pdf'])
%}
%% plot SVM accuracy cases during delay
%plot accuracy of each session
%{
frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
xlabelstr='Time(s) from delay onset';
y_lim=[0.4,1];
titlestr=cellfun(@(x,y) [x,'-',y],SVMchoice.animal,SVMchoice.date,'UniformOutput',false);
nroistr=num2cell(SVMchoice.nROI);
textstr=cellfun(@(x) [num2str(x),'cells'],nroistr,'UniformOutput',false);
figScoreCase=fScoreCase(SVMchoice.delayMovingSVM,SVMsensory.delayMovingSVM,frameNumTime,xlabelstr,y_lim,titlestr,textstr);
% plot SVM accuracy overall 
figScoreMean=fScoreMean(SVMchoice.delayMovingSVM,SVMsensory.delayMovingSVM,frameNumTime,xlabelstr,y_lim);
%}
%% plot SVM accuracy cases around go cue
%plot accuracy of each session
%{
frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
xlabelstr='Time(s) from go cue';
y_lim=[0.4,1];
figScoreCase=fScoreCase(SVMchoice.goMovingSVM,SVMsensory.goMovingSVM,frameNumTime,xlabelstr,y_lim,titlestr,textstr);
% plot SVM accuracy overall 
figScoreMean=fScoreMean(SVMchoice.goMovingSVM,SVMsensory.goMovingSVM,frameNumTime,xlabelstr,y_lim);
%}
%% plot cross temporal decoding accuracy of SVM around go cue
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
%% helper function

function markersize=fSigMarkerSize(p)
markersize(p>0.05)=0;
markersize(logical((p<=0.05).*(p>0.01)))=4;
markersize(logical((p<=0.01).*(p>0.001)))=8;
markersize(logical(p<=0.001))=12;
end

function figScoreCase=fScoreCase(scoredataC,scoredataS,frameNumTime,xlabelstr,y_lim,titlestr,textstr)
%plot SVM accuracy cases
colorScore={[1,0,0],[0,0,1],[0.5,0.5,0.5]};%choice,sensory,shuffle
frameNum=size(scoredataC(1).shuffle_do,2);
ts=-frameNumTime(1):sum(frameNumTime)/frameNum:frameNumTime(2);
figScoreCase=figure;
n_col=4;
n_row=ceil(size(scoredataC,1)/n_col);
set(gcf,'Position',[100,100,200*n_col,200*n_row]);
for i=1:size(scoredataC,1)
    subplot(n_row,n_col,i);
    fPlotMean_SE(ts,scoredataC(i).shuffle_do,colorScore{3});
    fPlotMean_SE(ts,scoredataS(i).shuffle_do,colorScore{3});
    fPlotMean_SE(ts,scoredataC(i).do,colorScore{1});
    fPlotMean_SE(ts,scoredataS(i).do,colorScore{2}); 
    plot([-1,ts(end)],[0.5,0.5],'k--');
    plot([0,0],y_lim,'k--');
    if contains(xlabelstr,'delay')
        plot([-0.5,-0.5],y_lim,'k--');
    end
    title(titlestr{i});
    text(0.1,0.9,textstr{i});
    set(gca,'Ylim',y_lim,'XTick',[-frameNumTime(1),0,floor(frameNumTime(2))],'XTickLabel',{num2str(-frameNumTime(1)),'0',num2str(floor(frameNumTime(2)))});
    box off;
    if mod(i,n_col)==1
        ylabel('Decoding accuracy');
        xlabel(xlabelstr);
    end
    set(gca,'FontSize',12);
end
end

function figScoreMean=fScoreMean(scoredataC,scoredataS,frameNumTime,xlabelstr,y_lim)
frameNum=size(scoredataC(1).shuffle_do,2);
ts=-frameNumTime(1):sum(frameNumTime)/frameNum:frameNumTime(2);
figScoreMean=figure;
colorScore={[1,0,0],[0,0,1],[0.5,0.5,0.5]};%choice,sensory,shuffle
set(gcf,'Position',[100,100,200,200]);
SVMsensoryMean_temp=arrayfun(@(X) nanmean(X.do),scoredataS,'UniformOutput',false);
n_2d=cellfun(@(x) size(x,2),SVMsensoryMean_temp);
SVMsensoryMean=cellfun(@(x) x(1:min(n_2d)),SVMsensoryMean_temp,'UniformOutput',false);
SVMsensoryMean=cell2mat(SVMsensoryMean);
fPlotMean_SE(ts,SVMsensoryMean,colorScore{2},'only cases');
hold on;
SVMchoiceyMean_temp=arrayfun(@(X) nanmean(X.do),scoredataC,'UniformOutput',false);
n_2d=cellfun(@(x) size(x,2),SVMchoiceyMean_temp);
SVMchoiceMean=cellfun(@(x) x(1:min(n_2d)),SVMchoiceyMean_temp,'UniformOutput',false);
SVMchoiceMean=cell2mat(SVMchoiceMean);
fPlotMean_SE(ts,SVMchoiceMean,colorScore{1},'only cases');
fPlotMean_SE(ts,SVMsensoryMean,colorScore{2});
fPlotMean_SE(ts,SVMchoiceMean,colorScore{1});
plot([-1,ts(end)],[0.5,0.5],'k--');
plot([0,0],y_lim,'k--');
if contains(xlabelstr,'delay')
    plot([-0.5,-0.5],y_lim,'k--');
end
set(gca,'Ylim',y_lim,'Xlim',[-1.1,1.6],'XTick',[-0.5,0,1],'XTickLabel',{'-0.5','0','1'});
box off;
ylabel('Decoding accuracy');
xlabel(xlabelstr);
set(gca,'FontSize',12);
box off;
nboot=1000;
binsize=1;
binstep=1;
pSigTtest=0.01;
baseline=0.5;
[ meanBoot, ciBoot, prctileBoot,activityBootBin ] = fMovingBootstrpMeanCIPrctile( SVMchoiceMean,nboot,binsize,binstep,baseline);
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
[ meanBoot, ciBoot, prctileBoot,activityBootBin ] = fMovingBootstrpMeanCIPrctile( SVMsensoryMean,nboot,binsize,binstep,baseline);
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
% 
% indSig=logical((prctileBoot<pSigTtest/2)+(prctileBoot>1-pSigTtest/2));
% 
% xSig(~indSig)=nan;
% ySig(~indSig)=nan;
% xSig=fRuleOutOccasional(xSig,1);
% ySig=fRuleOutOccasional(ySig,1);
% plot(xSig,ySig,'b-','LineWidth',1);
end

function [ meandata, sigmat,pmat] = fPlotMeanSig2D_SVM(ax,data,ts,pSig, baseline,strbar)
%FPLOTMEANSIG2D plot mean in color and draw contour lines ot represent
%significant range; Currently mainly used for showing SVM decoding accuracy
%Input-
%   data- r-by-m-by-n matrix, m rows, n columns and r repeats, m=n
%   pSig- p value that decide whether data was significant from baseline
%   baseline- double, compared with data
%Output-
%   ax- axis to plot
%   meandata- m-by-n matrix store the mean of each point
%   ts- time stamps 
%   sigmat- m-by-n matrix made of 0-1, which indicates which position is
%   significant
%   pmat- m-by-n matrix with each point indicates the p value of that
%   posision

meandata = nanmean(data,1);
pmat=ones(size(data,2),size(data,3));
% %this method seems to slow
% for i=1:size(data,2)
%     for j=1:size(data,3)
%         test_mean=bootstrp(1000,@mean,data(:,i,j));
%         temp_p=nansum(test_mean>baseline)/length(test_mean);
%         pmat(i,j)=1-2*abs(temp_p-0.5);
%     end
% end
for i=1:size(data,2)
    for j=1:size(data,3)
        [h,temp_p]=ttest(data(:,i,j)-baseline);
        pmat(i,j)=temp_p;
    end
end
sigmat= (pmat<pSig);
%plot mean
axes(ax);
clims=[0,1];
meandata=squeeze(meandata);
imagesc(ts,ts,meandata,clims);
hold on;
if strcmp(strbar,'colorbar')
    colormap(jet);
    chandle=colorbar;
    set(chandle,'Position',[0.91,0.1,0.02,0.8]);
    chandle.Label.String= 'Decoding accuracy';
end
%plot significance contour
contour(ts,ts,sigmat,1,'--','LineColor','w');
plot([0,0],[ts(1),ts(end)],'k--');
plot([ts(1),ts(end)],[0,0],'k--');
set(gca,'FontSize',12);

end


