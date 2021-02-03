function [ fig ] = fPlotImagingExample( fig,position )
%FPLOTIMAGINGEXAMPLE Summary of this function goes here

%decide some global variable
behEventAlign='first lick';%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
masklick='no';
frameNumTime=[2,3];%from 5s before align point to 5s after align point
roiNo=1;

CurrFolder = 'H:\2P\pyx095_20180810\im_data_reg\result_save';
rootname='pyx095_20180810';
load([CurrFolder,'\',rootname, '-imaging_Virables.mat']);%load behavior data
dirmat=strcat(CurrFolder,'\*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
file_imaging=cellfun(@(x) contains(x,'Ca'), filenames);
i_file_imaging=find(file_imaging);
load([CurrFolder,'\',filenames{1,i_file_imaging}]);%load imaging data

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
dff = (f - f0)/f0;
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

frameNum=double(round(frameNumTime*1000/frT));
[behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( Data_extract, ind_1stFrame, SavedCaTrials.FrameTime );%get behavior event time
%   align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
if strcmp(behEventAlign,'stimOnset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
    [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff, behEventFrameIndex,  frameNum );
else
    [ dff_aligned, behEvent_aligned, licking_aligned ] = fAlignSigBehEvent( dff, behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
end

trialType = fGetTrialType( Data_extract);%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
trialType(3:4,:,:)=[];%do not plot vio and miss
%%
%plot color plot,4 column(correct/error/miss/violation), 7 row(for each
%stimlus)+ mean trace
figure(fig);
set(fig, 'PaperPosition', position);
ColLimit = prctile(dff,90);
titlestr={'Correct','Error','Miss','Violation'};
titlestr=strcat('ROI-',num2str(roiNo),titlestr);
color_mean_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};
for nStim=1:size(trialType,2) %for each stimulus
    for  nResult=1:size(trialType,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
        selectedTrialInd=trialType(nResult,nStim,:);
        selectedTrialInd=logical(squeeze(selectedTrialInd))';
        sot=behEvent_aligned.stimOnset(selectedTrialInd);% stim onset, white
        flt=behEvent_aligned.lickFirst(selectedTrialInd);% first lick, black
        flt_l=behEvent_aligned.lickFirst_left(selectedTrialInd);%left lick, magenta dots
        flt_r=behEvent_aligned.lickFirst_right(selectedTrialInd);%right lick, red dots
        rwt=behEvent_aligned.rewTime(selectedTrialInd);%reward, cyan dots
        go=behEvent_aligned.go(selectedTrialInd);%go cue,white
        neuralActivity=dff_aligned(selectedTrialInd,:);
        temp=find(selectedTrialInd);
        leftLick=cell(1,length(temp));
        rightLick=cell(1,length(temp));
        for i=1:length(temp)
            leftLick{i}= licking_aligned.leftLick{temp(i)};
            rightLick{i}=licking_aligned.rightLick{temp(i)};
        end
        %decide which behavior event to be sort
        if strcmp(behEventSort,'first lick')
            [B,I]=sort(flt);
        elseif strcmp(behEventSort,'reward')
            [B,I]=sort(rwt);
        elseif strcmp(behEventSort,'go cue')
            [B,I]=sort(go);
        end
        go=go(I);
        sot=sot(I);
        flt=flt(I);
        flt_l=flt_l(I);
        flt_r=flt_r(I);
        rwt=rwt(I);
        subplot(size(trialType,2)+1,2*size(trialType,1),2*nResult-1+2*size(trialType,1)*(nStim-1));
        hold on;
        imagesc(neuralActivity(I,:));
        x_lim=[0,size(neuralActivity,2)];%get(gca,'Xlim');
        plot(sot,[1:sum(selectedTrialInd)],'w.','markersize',8);
        hold on;
        plot(go,[1:sum(selectedTrialInd)],'w.','markersize',8);
        hold on;
        plot(flt,[1:sum(selectedTrialInd)],'k.','markersize',8);
        hold on;
        plot(flt_l,[1:sum(selectedTrialInd)],'m.','markersize',8);
        hold on;
        plot(flt_r,[1:sum(selectedTrialInd)],'r.','markersize',8);
        hold on;
        plot(rwt,[1:sum(selectedTrialInd)],'c.','markersize',8);
        hold on;
        set(gca,'clim',[0 ColLimit]);
        set(gca,'ytick',sum(selectedTrialInd),'yticklabel',sum(selectedTrialInd),'ydir','normal');
        set(gca,'xtick',[]);
        ColBar = colorbar;
        set(ColBar,'position',[0.91 0.2400 0.02 0.1445]);
        if sum(selectedTrialInd)>0
            ylim([0.5 sum(selectedTrialInd)+0.5]);
        end
        if nStim == 1
            title(titlestr{nResult});
        end
        if nResult==1
            ylabel([num2str(Data_extract.Stimuli(nStim))]); %' clicks/s'
        end
        xlim(x_lim);

        subplot(size(trialType,2)+1,2*size(trialType,1),2*nResult+2*size(trialType,1)*(nStim-1));
        for i=1:length(leftLick)%each trial as a row
            for jl=1:length(leftLick{I(i)})
                line([leftLick{I(i)}(jl),leftLick{I(i)}(jl)],[i-0.5,i+0.5],'color','b','linewidth',1);%left lick
                hold on;
            end
            for jr=1:length(rightLick{I(i)})
                line([rightLick{I(i)}(jr),rightLick{I(i)}(jr)],[i-0.5,i+0.5],'color','r','linewidth',1);%right lick
                hold on;
            end
        end
        set(gca,'ytick',[]);%licking raster have corresponding color plot label, so not need
        set(gca,'xtick',[]);
        xlim(x_lim);
        if sum(selectedTrialInd)>0
            ylim([0.5 sum(selectedTrialInd)+0.5]);
        end
        
        %plot mean trace
        subplot(size(trialType,2)+1,2*size(trialType,1),2*nResult-1+2*size(trialType,1)*(size(trialType,2)));
        [neuralActivityMean, neuralActivityCI]=fMean_CI(neuralActivity,0.05);
        %plot CI
        %{
            xpatch=[1:size(neuralActivity,2), fliplr(1:size(neuralActivity,2))];
            ypatch=[neuralActivityCI(1,:),fliplr(neuralActivityCI(2,:))];
            p=patch(xpatch,ypatch,color_mean_trace{nStim});%plot confidence interval
            p.FaceAlpha=0.1;
            p.EdgeColor=color_mean_trace{nStim};%'none';
            hold on;
        %}
        %plot mean
        curve_meanTrace(nStim)=plot(1:size(neuralActivity,2),neuralActivityMean,'Color',color_mean_trace{nStim},'linewidth',2);
        hold on;
        y_lim=get(gca,'Ylim');
        xlim(x_lim);
        plot([frameNumTime(1)*round(1000/frT),frameNumTime(1)*round(1000/frT)],[y_lim(1),y_lim(2)],'k-');%align to a behavior event
        xlabel(['time(s) from ',behEventAlign]);
        if nResult==1
            ylabel('\it\DeltaF/F');
        end
        set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
        set(gca,'ylim',[0,3]);
%         set(gca,'FontName','Arial','FontSize',14);
    end
end
subplot(size(trialType,2)+1,2*size(trialType,1),2*size(trialType,1)*(size(trialType,2)+1)-1);
leg=legend(curve_meanTrace(:),num2str(Data_extract.Stimuli(:)));
set(leg,'position',[0.91 0.05 0.02 0.14]);
legend boxoff 
end

