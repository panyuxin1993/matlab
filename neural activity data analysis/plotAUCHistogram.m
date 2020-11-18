cd 'D:\data\pyx224_20190729\result_save';
CurrFolder = pwd;
rootname='pyx224_20190729';
savefolder=rootname;%'1-200trials';%'segNP';
load([CurrFolder,'\',rootname, '-imaging_Virables.mat']);%load behavior data
dirmat=strcat(CurrFolder,'\*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
file_imaging=cellfun(@(x) contains(x,'Ca'), filenames);
i_file_imaging=find(file_imaging);
load(filenames{1,i_file_imaging});%load imaging data
load('dff.mat');

if ~exist(savefolder)
    mkdir(savefolder);
end
varTypes = {'double','double','double','double','double'};
nROI=size(SavedCaTrials.f_raw{1},1);
T_AUC=table('Size',[nROI,5],'VariableTypes',varTypes,...
    'VariableNames',{'ITI','sound','delay','response','lick'});
T_pvalues=table('Size',[nROI,5],'VariableTypes',varTypes,...
    'VariableNames',{'ITI','sound','delay','response','lick'});
ind_tr_1=1;
ntr=length(SavedCaTrials.f_raw);
frT = SavedCaTrials.FrameTime;
for roiNo = 1:size(SavedCaTrials.f_raw{1},1) %SavedCaTrials.nROIs may be not true
%%
    % align to behavior event
    nFrameEachTrial=cellfun(@(x) size(x,2),SavedCaTrials.f_raw);
    ind_1stFrame=zeros(1,length(nFrameEachTrial));
    ind_1stFrame(1)=1;
    ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
    ind_1stFrame=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials
    frameNumTime=[2,2];%from 5s before align point to 5s after align point
    frameNum=double(round(frameNumTime*1000/frT));
    [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( Data_extract, ind_1stFrame, SavedCaTrials.FrameTime );%get behavior event time
    T_SigbyEpoch = fGetSigBehEpoch(behEventFrameIndex,dff(roiNo,:),frT);
    [trialType,rule] = fGetTrialType( Data_extract);%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
    label_choice=double(Data_extract.Action_choice)';%0,1,2,3
    ind_trial=(label_choice<2);%correct/error trials
%     ind_trial=logical(reshape(sum(trialType(1,:,:),2),[],1));%only correct trials
    poslabel=1;
    [T_AUC.ITI(roiNo),T_pvalues.ITI(roiNo)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.ITI(ind_trial),poslabel,1000);
    [T_AUC.sound(roiNo),T_pvalues.sound(roiNo)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.sound(ind_trial),poslabel,1000);
    [T_AUC.delay(roiNo),T_pvalues.delay(roiNo)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.delay(ind_trial),poslabel,1000);
    [T_AUC.response(roiNo),T_pvalues.response(roiNo)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.response(ind_trial),poslabel,1000);
    [T_AUC.lick(roiNo),T_pvalues.lick(roiNo)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.lick(ind_trial),poslabel,1000);
end
save('table_AUC_Correct','T_AUC','T_pvalues');
%%
load('table_AUC.mat');
psig=[0.005,0.995];
figAUC=figure;
set(gcf,'Position',[100,100,800,200]);
subplot(1,5,1);
histogram(T_AUC.ITI,'BinWidth',0.05);
hold on;
histogram(T_AUC.ITI(T_pvalues.ITI<psig(1)),'FaceColor','r','BinWidth',0.05);
histogram(T_AUC.ITI(T_pvalues.ITI>psig(2)),'FaceColor','m','BinWidth',0.05);
title('AUC during ITI');
set(gca,'Xlim',[0,1]);
subplot(1,5,2);
histogram(T_AUC.sound,'BinWidth',0.05);
hold on;
histogram(T_AUC.sound(T_pvalues.sound<psig(1)),'FaceColor','r','BinWidth',0.05);
histogram(T_AUC.sound(T_pvalues.sound>psig(2)),'FaceColor','m','BinWidth',0.05);
title('AUC during sound');
set(gca,'Xlim',[0,1]);
subplot(1,5,3);
histogram(T_AUC.delay,'BinWidth',0.05);
hold on;
histogram(T_AUC.delay(T_pvalues.delay<psig(1)),'FaceColor','r','BinWidth',0.05);
histogram(T_AUC.delay(T_pvalues.delay>psig(2)),'FaceColor','m','BinWidth',0.05);
title('AUC during delay');
set(gca,'Xlim',[0,1]);
subplot(1,5,4);
histogram(T_AUC.response,'BinWidth',0.05);
hold on;
histogram(T_AUC.response(T_pvalues.response<psig(1)),'FaceColor','r','BinWidth',0.05);
histogram(T_AUC.response(T_pvalues.response>psig(2)),'FaceColor','m','BinWidth',0.05);
title('AUC during response');
set(gca,'Xlim',[0,1]);
subplot(1,5,5);
histogram(T_AUC.lick,'BinWidth',0.05);
hold on;
histogram(T_AUC.lick(T_pvalues.lick<psig(1)),'FaceColor','r','BinWidth',0.05);
histogram(T_AUC.lick(T_pvalues.lick>psig(2)),'FaceColor','m','BinWidth',0.05);
title('AUC during lick');
set(gca,'Xlim',[0,1]);
saveas(figAUC,['AUC of at significance level of ',num2str(2*psig(1)),'.png'],'png');