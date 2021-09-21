%compare opto inhibition effect, both amplitude and selectivity
%only analyze correct trials
%amplitude- mean of both side, opto vs. control
%selectivity- difference of contra and ipsi (absolute value)/ AUC
close all;
%decide some global variable
behEventAlign='delay onset';%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
masklick='no';
frameNumTime=[1,3];%from 5s before align point to 5s after align point
yrange=[ -0.5 , 1 ];
i_selectivity=4;%*********variable**************
selectivitystr={'stimuli','difficulty','sensory','choice'};%sensory means grouping difficulties;

rootname={'pyx349_20210504'};
celltypestr={'vgat'};%{'vgat','vglut2'};
figCmpCtrlOpto=figure;
set(gcf,'PaperPosition',[1,1,5,2]);
for i=1:length(celltypestr)
    rootpath=['H:\2P\',rootname{i},'\im_data_reg\result_save'];
    cd(rootpath);
    
    CurrFolder = pwd;
    
    savefolder=rootname{i};%'1-200trials';%'segNP';
    % load([CurrFolder,'\',rootname, '_imaging_Virables.mat']);%load behavior data
    dirmat=strcat(CurrFolder,'\*.mat');
    dirs=dir(dirmat);
    dircell=struct2cell(dirs);
    filenames=dircell(1,:);
    file_imaging=cellfun(@(x) contains(x,'Ca'), filenames);
    i_file_imaging=find(file_imaging);
    load(filenames{1,i_file_imaging});%load imaging data
    file_beh=cellfun(@(x) contains(x,'imaging_Virables'), filenames);
    i_file_beh=find(file_beh);
    load(filenames{1,i_file_beh});%load behavior data
    
    if ~exist(savefolder,'dir')
        mkdir(savefolder);
    end
    if exist('dff.mat','file')
        load([CurrFolder,'\','dff.mat']);%load dff
    else
        dff=fnx_getDff(CurrFolder,savefolder,'save figure');
    end
    
    ind_tr_1=1;%using data from trial 1
    ntr = length(SavedCaTrials.f_raw); %use all data, or  just set the number of trials to use
    frT = SavedCaTrials.FrameTime;
    frameRate=1000/frT;
    % align to behavior event
    nFrameEachTrial=cellfun(@(x) size(x,2),SavedCaTrials.f_raw);
    ind_1stFrame=zeros(1,length(nFrameEachTrial));
    ind_1stFrame(1)=1;
    ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
    ind_1stFrame=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials
    frameNum=double(round(frameNumTime*frameRate));
    frameNumGo=double(round(1*frameRate));%delay length is 1s
    [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( Data_extract, ind_1stFrame, 1000/frameRate ,ind_tr_1);%get behavior event time
    
    
    [trialType,rule] = fGetTrialType( Data_extract,[],i_selectivity,'matrix','left','divideCorErr','divideOpto');%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
    
    meanResponse=cell(size(SavedCaTrials.f_raw{1},1),size(trialType,2),size(trialType,4));
    
    for roiNo = 1:size(SavedCaTrials.f_raw{1},1)%for each cell
        if strcmp(behEventAlign,'delay onset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
            [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff(roiNo,:), behEventFrameIndex,  frameNum );
        else
            [ dff_aligned, behEvent_aligned, licking_aligned ] = fAlignSigBehEvent( dff(roiNo,:), behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
        end
        for iChoice=1:size(trialType,2)
            for iOpto=1:size(trialType,4)
                selectedTrialInd=trialType(1,iChoice,:,iOpto);
                selectedTrialInd=logical(squeeze(selectedTrialInd))';
                meanResponse{roiNo,iChoice,iOpto}=dff_aligned(selectedTrialInd,frameNum(1)+1:frameNum(1)+1+frameNumGo);
            end
        end
    end
    
    meanActivity=cellfun(@(x) nanmean(nanmean(x)), meanResponse);
    meanActivityCtrlOpto=squeeze(nanmean(meanActivity,2));
    selectivityCtrlOpto=squeeze(abs(meanActivity(:,2,:)-meanActivity(:,1,:)));
    figure(figCmpCtrlOpto);
    subplot(2,2,i*2-1);
    curveMean=fScatterPlot(meanActivityCtrlOpto(:,1),meanActivityCtrlOpto(:,2));
    xlabel('\DeltaF/F');
    ylabel('{\DeltaF/F} with M2-SC opto-inactivation');
    title(celltypestr{i});
    subplot(2,2,i*2);
    curveSelectivity=fScatterPlot(selectivityCtrlOpto(:,1),selectivityCtrlOpto(:,2));
    xlabel('selectivity (\DeltaF/F)');
    ylabel('selectivity {\DeltaF/F} with M2-SC opto-inactivation');
    title(celltypestr{i});
end
% saveas(figCmpCtrlOpto,[savefolder,filesep,'opto effect.pdf'],'pdf');

function [curve]=fScatterPlot(x,y)
curve=scatter(x,y,10,'k','filled');hold on;
temp1=min(min(x),min(y));
temp2=max(max(x),max(y));
plot([temp1,temp2],[temp1,temp2]);

set(gca,'Xlim',[temp1,temp2],'Ylim',[temp1,temp2]);
if adtest(x) && adtest(y)
    [h,p]=ttest(x,y);
    text(0,1,['t test p=',num2str(p)],'Unit','Normalized');
else
    p=signrank(x,y);
    text(0,1,['signed rank test p=',num2str(p)],'Unit','Normalized');
end
set(gca,'FontSize',12,'FontName','Arial');
end