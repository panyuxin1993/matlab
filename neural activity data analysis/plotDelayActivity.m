%only plot mean trace of delay activity

%decide some global variable
behEventAlign='delay onset';%align to which event(string can be in {'delay onset', 'go cue'},

cd 'F:\pyx095_20180810\im_data_reg_cpu\result_save';
CurrFolder = pwd;
rootname='pyx095_20180810';
load([CurrFolder,'\',rootname, '-imaging_Virables.mat']);%load behavior data
dirmat=strcat(CurrFolder,'\*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
file_imaging=cellfun(@(x) contains(x,'Ca'), filenames);
i_file_imaging=find(file_imaging);
load(filenames{1,i_file_imaging});%load imaging data

if ~exist(rootname)
    mkdir(rootname);
end

f_cat  = [];
ind_1stFr(1) = 1;
ntr = length(SavedCaTrials.f_raw);
for i = 1:ntr
    f_cat = [f_cat SavedCaTrials.f_raw{i}];
    ind_1stFr(i+1) = size(f_cat,2) + 1;
end
ind_1stFr(i+1) = [];

for roiNo = 1:size(SavedCaTrials.f_raw{1},1) %SavedCaTrials.nROIs may be not true
%%
%substract changing baseline
    totalFr = size(f_cat,2);
    frT = SavedCaTrials.FrameTime;
    span = round(20*1000/frT/2); %every 40s, get a mean and substrate to compensate for baseline change
    ind_x = 0;
    x = [];
    for i = 1:totalFr
        %     ind_x = ind_x + 1;
        ind1 = i- span;
        ind2 = i+span;
        if ind1 <= 0
            ind1 = 1;
        end
        if ind2 > totalFr
            ind2 = totalFr;
        end
        x(i) = prctile(f_cat(roiNo,ind1:ind2),5);
    end
    f = f_cat(roiNo,:) - x + mean(x);
    % figure; histogram(f_cat_sub(roiNo,:),100)   
    %%
    %get f mode as f0
    [N,edges,bin] = histcounts(f,100);
    f0 = edges(N == max(N));
    if length(f0)==2
        f0=mean(f0);
    end
    %get dff
    dff = (f - f0)/f0*100;
%     B=figure; 
%     set(gcf, 'position', [0 0 1500 500]);
%     plot(f_cat(roiNo,:),'k-');
%     hold on;
%     plot(x,'c-');
%     plot(f,'g-');
%     plot(dff,'b-');
%     legend('f cat','moving f mean','f baseline correction','dff');
%     set(gca,'FontName','Arial','FontSize',14);
%     saveas(B,[CurrFolder,'\',rootname,'\',rootname,'ROI-',num2str(roiNo),'-cat_f.jpg'],'jpg');
    % align to behavior event
    nFrameEachTrial=cellfun(@(x) size(x,2),SavedCaTrials.f_raw);
    ind_1stFrame=zeros(1,length(nFrameEachTrial));
    ind_1stFrame(1)=1;
    for i=2:length(nFrameEachTrial)
        ind_1stFrame(i)=ind_1stFrame(i-1)+nFrameEachTrial(i-1);
    end
    frameNumTime=[1.5,1.5];%from 5s before align point to 5s after align point
    frameNum=double(round(frameNumTime*1000/frT));
    [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( Data_extract, ind_1stFrame, SavedCaTrials.FrameTime );%get behavior event time
   %   align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
    [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff, behEventFrameIndex,frameNum );%extract delay activity
    trialType = fGetTrialType( Data_extract);%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
    %%
    %plot  mean trace
    figMeanTrace=figure;%plot mean trace
    set(gcf, 'position', [0 0 1400 300]);
    ColLimit = prctile(dff,90);
    titlestr={'Correct','Error','Miss','Violation'};
    titlestr=strcat('ROI-',num2str(roiNo),titlestr);
    color_mean_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};
    for nStim=1:size(trialType,2) %for each stimulus
        for  nResult=1:size(trialType,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
            selectedTrialInd=trialType(nResult,nStim,:);
            selectedTrialInd=logical(squeeze(selectedTrialInd))';
%             sot=behEvent_aligned.stimOnset(selectedTrialInd);% stim onset, white
%             flt=behEvent_aligned.lickFirst(selectedTrialInd);% first lick, black
%             flt_l=behEvent_aligned.lickFirst_left(selectedTrialInd);%left lick, magenta dots
%             flt_r=behEvent_aligned.lickFirst_right(selectedTrialInd);%right lick, red dots
%             rwt=behEvent_aligned.rewTime(selectedTrialInd);%reward, cyan dots
%             go=behEvent_aligned.go(selectedTrialInd);%go cue,white
            neuralActivity=dff_aligned(selectedTrialInd,:);
            %plot mean trace
            figure(figMeanTrace);%save mean trace
            %subplot(size(trialType,2)+3,2*size(trialType,1),2*nResult-1+2*size(trialType,1)*(size(trialType,2)+2));
            subplot(1,size(trialType,1),nResult);
            [neuralActivityMean, neuralActivityCI]=fMean_CI(neuralActivity,0.05);
%             xpatch=[1:size(neuralActivity,2), fliplr(1:size(neuralActivity,2))];
%             ypatch=[neuralActivityCI(1,:),fliplr(neuralActivityCI(2,:))];
%             p=patch(xpatch,ypatch,color_mean_trace{nStim});%plot confidence interval
%             p.FaceAlpha=0.1;
%             p.EdgeColor=color_mean_trace{nStim};%'none';
%             hold on;
            curve_meanTrace(nStim)=plot(1:size(neuralActivity,2),neuralActivityMean,'Color',color_mean_trace{nStim},'linewidth',2);
            hold on;
            title(titlestr{nResult});
            y_lim=get(gca,'Ylim');
%             xlim([1,size(neuralActivity,2)]);
            plot([frameNumTime(1)*round(1000/frT),frameNumTime(1)*round(1000/frT)],[y_lim(1),y_lim(2)],'k-');%align to a behavior event
            xlabel(['time(s) from ',behEventAlign]);
            ylabel('df/f');
            set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);    
            set(gca,'FontName','Arial','FontSize',14);
        end
    end
    figure(figMeanTrace);
    subplot(1,size(trialType,1),1);
    legend(curve_meanTrace(:),num2str(Data_extract.Stimuli(:)),'Location','best');
    saveas(figMeanTrace,[CurrFolder,'\',rootname,'\',rootname,'ROI-',num2str(roiNo),'delay activity-algin to ',behEventAlign,'-mean trace.jpg'],'jpg');
    close all;
end
%%
% figure(gcf);
% plot((1:size(dff_aligned,2))*frT/1000,dff_aligned(55,:))
