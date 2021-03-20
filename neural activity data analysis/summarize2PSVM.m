close all;
[num,txt,raw] =xlsread('D:\xulab\project\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
savepath='H:\2P\summary\summary_DLCfiltered';%'H:\2P\summary';
clear TSVM_combine;
trialTypeStr='cor and err';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
SVMtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
AUCCorrectedMethod='balencedCorErrTrialNum';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
nRepeat=100;
pTraining=0.9;
%% calculate SVM accuracy
%
T=cell2table(raw(2:end,1:14));
T.Properties.VariableNames=strrep(raw(1,1:14),' ','_');%table variable name can't have ' ',so replace them
ind_session=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control');
ind_session=find(ind_session);
n_session=length(ind_session);
animal_unique=unique(T.animal(ind_session));
n_animal=length(animal_unique);

score=cell(n_session,1);
for i_session=1:n_session
    indrow=ind_session(i_session);
    TSVM_currentSession= fGetEpochSVMscoreASession(T.file_path{indrow},T.animal{indrow},T.date{indrow},T.field{indrow},T.cell_type{indrow},trialTypeStr,SVMtype,nRepeat, pTraining,AUCCorrectedMethod);
    if exist('TSVM_combine','var') 
        TSVM_combine=vertcat(TSVM_combine,TSVM_currentSession);
    else
        TSVM_combine=TSVM_currentSession;
    end
end
save([savepath,filesep,'trialType',trialTypeStr,'-',SVMtype,'pTraining',num2str(pTraining),'-TepochSVM.mat'],'TSVM_combine');
%}
%% plot SVM accuracy cases
celltype='vglut2';
SVMtype='sensory';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
load([savepath,filesep,'trialType',trialTypeStr,'-',SVMtype,'pTraining',num2str(pTraining),'-TepochSVM.mat']);
SVMsensory=TSVM_combine(strcmp(TSVM_combine.celltype,celltype),:);

SVMtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
load([savepath,filesep,'trialType',trialTypeStr,'-',SVMtype,'pTraining',num2str(pTraining),'-TepochSVM.mat']);
SVMchoice=TSVM_combine(strcmp(TSVM_combine.celltype,celltype),:);

%plot accuracy of each session
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
    fPlotMean_SE(1:size(SVMchoice(i,5:9),2),cell2mat(table2array(SVMchoice(i,5:9))),colorScore{1});
    fPlotMean_SE(1:size(SVMsensory(i,5:9),2),cell2mat(table2array(SVMsensory(i,5:9))),colorScore{2}); 
    plot([1,5],[0.5,0.5],'k--');
    title([SVMchoice.animal{i},'-',SVMchoice.date{i}]);
    set(gca,'Ylim',[0.3,1],'Xlim',[0.9,5.1],'XTick',1:5,'XTickLabel',{'ITI','sound','delay','response','lick'});
    set(gca,'FontSize',12);
    h=gca;
    h.XTickLabelRotation=45;
    box off;
end
set(gcf,'PaperPosition',[1,1,2,2]);
saveas(figScoreCase,[savepath,filesep,celltype,'-trialType',trialTypeStr,'-',SVMtype,'pTraining',num2str(pTraining),'-TepochSVM_score',SVMchoice.animal{i},'-',strrep(SVMchoice.date{i},'/',''),'.pdf']);

%% plot SVM accuracy overall
figMean=figure;
set(gcf,'Position',[100,100,200,200]);
temp=table2array(SVMsensory(:,5:9));
SVMsensoryMean=cellfun(@nanmean,temp);
fPlotMean_SE(1:size(SVMsensoryMean,2),SVMsensoryMean,colorScore{2},'only cases');
temp=table2array(SVMchoice(:,5:9));
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

% set(gcf,'PaperPosition',[1,1,2,2]);
% saveas(figMean,[savepath,filesep,celltype,'trialType',trialTypeStr,'-',SVMtype,'pTraining',num2str(pTraining),'-TepochSVM_scoreMean.pdf'])

%% plot SVM accuracy cases during delay
%plot accuracy of each session
colorScore={[1,0,0],[0,0,1],[0.5,0.5,0.5]};%choice,sensory,shuffle
frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
frameNum=size(SVMchoice.delayMovingSVM(i).shuffle_do,2);
ts=-1:2.5/frameNum:1.5;
pSig=0.05;
figScoreCasedelay=figure;
n_col=4;
n_row=ceil(size(SVMchoice,1)/n_col);
set(gcf,'Position',[100,100,200*n_col,200*n_row]);
for i=1:size(SVMchoice,1)
    subplot(n_row,n_col,i);
    fPlotMean_SE(ts,SVMchoice.delayMovingSVM(i).shuffle_do,colorScore{3});
    fPlotMean_SE(ts,SVMsensory.delayMovingSVM(i).shuffle_do,colorScore{3});
    fPlotMean_SE(ts,SVMchoice.delayMovingSVM(i).do,colorScore{1});
    fPlotMean_SE(ts,SVMsensory.delayMovingSVM(i).do,colorScore{2}); 
    plot([-1,ts(end)],[0.5,0.5],'k--');
    plot([0,0],[0.3,0.8],'k--');
    plot([-0.5,-0.5],[0.3,0.8],'k--');
    set(gca,'Ylim',[0.3,0.8],'XTick',[-0.5,0,1],'XTickLabel',{'-0.5','0','1'});
    box off;
    if mod(i,n_col)==1
        ylabel('Decoding accuracy');
        xlabel('Time(s) from delay onset');
    end
    set(gca,'FontSize',12);
end

%% plot SVM accuracy overall during delay
figMeanDelay=figure;
set(gcf,'Position',[100,100,250,250]);
SVMsensoryMean_temp=arrayfun(@(X) nanmean(X.do),SVMsensory.delayMovingSVM,'UniformOutput',false);
n_2d=cellfun(@(x) size(x,2),SVMsensoryMean_temp);
SVMsensoryMean=cellfun(@(x) x(1:min(n_2d)),SVMsensoryMean_temp,'UniformOutput',false);
SVMsensoryMean=cell2mat(SVMsensoryMean);
fPlotMean_SE(ts,SVMsensoryMean,colorScore{2},'only cases');
hold on;
SVMchoiceyMean_temp=arrayfun(@(X) nanmean(X.do),SVMchoice.delayMovingSVM,'UniformOutput',false);
n_2d=cellfun(@(x) size(x,2),SVMchoiceyMean_temp);
SVMchoiceMean=cellfun(@(x) x(1:min(n_2d)),SVMchoiceyMean_temp,'UniformOutput',false);
SVMchoiceMean=cell2mat(SVMchoiceMean);
fPlotMean_SE(ts,SVMchoiceMean,colorScore{1},'only cases');
fPlotMean_SE(ts,SVMsensoryMean,colorScore{2});
fPlotMean_SE(ts,SVMchoiceMean,colorScore{1});
plot([-1,ts(end)],[0.5,0.5],'k--');
plot([0,0],[0.3,0.8],'k--');
plot([-0.5,-0.5],[0.3,0.8],'k--');
set(gca,'Ylim',[0.3,0.8],'Xlim',[-1.1,1.6],'XTick',[-0.5,0,1],'XTickLabel',{'-0.5','0','1'});
box off;
ylabel('Decoding accuracy');
xlabel('Time(s) from delay onset');
set(gca,'FontSize',12);
box off;
nboot=1000;
binsize=1;
binstep=1;
pSigTtest=0.01;
baseline=0.5;
[ meanBoot, ciBoot, prctileBoot,activityBootBin ] = fMovingBootstrpMeanCIPrctile( SVMchoiceMean,nboot,binsize,binstep,baseline);

indSig=logical((prctileBoot<pSigTtest/2)+(prctileBoot>1-pSigTtest/2));
y_lim=get(gca,'Ylim');
ySig=ones(size(prctileBoot))*y_lim(2)*0.9;
xSig=ts(1:length(prctileBoot));
xSig(~indSig)=nan;
ySig(~indSig)=nan;
xSig=fRuleOutOccasional(xSig,1);
ySig=fRuleOutOccasional(ySig,1);
plot(xSig,ySig,'b-','LineWidth',1);

