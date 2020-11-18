%continuing sampling, so link each trials one by one to get a 'long trial'
%averaging all dff of each cell 
% load('F:\pyx083-20180522\im_data_reg_cpu\result_save\CaTrialsSIM_pyx083_20180522_920nm_power50_zoom4x_dftReg_.mat');
% load('D:\xulab\behavior\pyx083\Data_Virables\2018_05_22_pyx083-imaging_Virables.mat');

%decide some global variable
behEventAlign='start';%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
activitySort='peak';% string can be in{'peak'};
masklick='no';
frameNumTime=[2,5];%from 5s before align point to 5s after align point
    
cd 'F:\2P\pyx255_20191206\im_data_reg\result_save';
CurrFolder = pwd;
rootname='pyx255_20191206';
savefolder=rootname;%'1-200trials';%'segNP';
load([CurrFolder,'\',rootname, '_imaging_Virables.mat']);%load behavior data
dirmat=strcat(CurrFolder,'\*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
file_imaging=cellfun(@(x) contains(x,'Ca'), filenames);
i_file_imaging=find(file_imaging);
load(filenames{1,i_file_imaging});%load imaging data

if ~exist(savefolder)
    mkdir(savefolder);
end
load([CurrFolder,'\','dff.mat']);%load dff

ind_tr_1=1;%using data from trial 1
ntr = length(SavedCaTrials.f_raw); %use all data, or  just set the number of trials to use
frT = SavedCaTrials.FrameTime;
% align to behavior event
nFrameEachTrial=cellfun(@(x) size(x,2),SavedCaTrials.f_raw);
ind_1stFrame=zeros(1,length(nFrameEachTrial));
ind_1stFrame(1)=1;
ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
ind_1stFrame=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials

frameNum=double(round(frameNumTime*1000/frT));
[behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( Data_extract, ind_1stFrame, SavedCaTrials.FrameTime ,ind_tr_1);%get behavior event time
[trialType,rule] = fGetTrialType( Data_extract);%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials

dff_all=zeros(size(SavedCaTrials.f_raw{1},1),sum(frameNum)+1,size(trialType,2),size(trialType,1));
%%
%plot color plot,4 column(correct/error/miss/violation), 7 row(for each
%stimlus)+ mean trace
A=figure;
set(gcf, 'position', [0 0 1400 600]);

if size(trialType,2)==6
    color_mean_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};
elseif size(trialType,2)==2
    color_mean_trace={[0 0 1],[1 0 0]};
end
if strcmp(rule,'low click rate-right')
    color_mean_trace=fliplr(color_mean_trace);
end
for roiNo = 1:size(SavedCaTrials.f_raw{1},1) %SavedCaTrials.nROIs may be not true
   %   align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
   if strcmp(behEventAlign,'stimOnset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
       [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff(roiNo,:), behEventFrameIndex,  frameNum );
   else
       [ dff_aligned, behEvent_aligned, licking_aligned ] = fAlignSigBehEvent( dff(roiNo,:), behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
   end
   

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
            dff_all(roiNo,:,nStim,nResult)=mean(neuralActivity);
            temp=find(selectedTrialInd);
            leftLick=cell(1,length(temp));
            rightLick=cell(1,length(temp));
            for i=1:length(temp)
                leftLick{i}= licking_aligned.leftLick{temp(i)};
                rightLick{i}=licking_aligned.rightLick{temp(i)};
            end            
        end
    end
end
for nStim=1:size(trialType,2) %for each stimulus
    for  nResult=1:size(trialType,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
        figure(A);
        subplot(size(trialType,2),2*size(trialType,1),2*nResult-1+2*size(trialType,1)*(nStim-1));
        hold on;
        for roiNo=1:size(dff_all,1)
            t_peak(roiNo)=find(dff_all(roiNo,:,nStim,nResult)==max(dff_all(roiNo,:,nStim,nResult)));
        end
        %decide which behavior event to be sort
        if strcmp(activitySort,'peak')
            [B,I]=sort(t_peak);
        end
        dff_stim_result=dff_all(:,:,nStim,nResult);
        imagesc(dff_stim_result);
        x_lim=[0,size(dff_all,2)];%get(gca,'Xlim');
        ColLimit = prctile(dff(roiNo,:),95);
        titlestr={'Correct','Error','Miss','Violation'};
        set(gca,'clim',[0 ColLimit]);
        set(gca,'ytick',size(dff_all,1),'yticklabel',size(dff_all,1),'ydir','normal');
        if nStim==size(trialType,2)
            set(gca,'xtick',[1:round(1000/frT):size(dff_all,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
            xlabel(['time(s) from ',behEventAlign],'FontName','Arial','FontSize',14);
        else
            set(gca,'xtick',[]);
        end
        ColBar = colorbar;
        set(ColBar,'position',[0.91 0.1100 0.02 0.1445]);
        
        if nStim == 1
            title(titlestr{nResult});
        end
        if nResult==1
            ylabel([num2str(Data_extract.Stimuli(nStim))]); %' clicks/s'
        end
        xlim(x_lim);
        ylim([0,size(dff_all,1)]);
    end
end
%             if nStim == 1 %plot neural activity; here end trials have more trials so wider
%                 subplot(size(trialType,2)+3,2*size(trialType,1),[2*nResult-1,2*nResult-1+2*size(trialType,1)]);
%             elseif nStim  == size(trialType,2)
%                 subplot(size(trialType,2)+3,2*size(trialType,1),[2*nResult-1+2*size(trialType,1)*size(trialType,2),2*nResult-1+2*size(trialType,1)*size(trialType,2)+2*size(trialType,1)]);
%             else
%                 subplot(size(trialType,2)+3,2*size(trialType,1),2*nResult-1+2*size(trialType,1)*nStim);
%             end

            
%             if nStim == 1 %plot licking raster
%                 subplot(size(trialType,2)+3,2*size(trialType,1),[2*nResult,2*nResult+2*size(trialType,1)]);
%             elseif nStim  == size(trialType,2)
%                 subplot(size(trialType,2)+3,2*size(trialType,1),[2*nResult+2*size(trialType,1)*size(trialType,2),2*nResult+2*size(trialType,1)*size(trialType,2)+2*size(trialType,1)]);
%             else
%                 subplot(size(trialType,2)+3,2*size(trialType,1),2*nResult+2*size(trialType,1)*nStim);
%             end
%             subplot(size(trialType,2),2*size(trialType,1),2*nResult+2*size(trialType,1)*(nStim-1));
%             for i=1:length(leftLick)%each trial as a row
%                 for jl=1:length(leftLick{I(i)})
%                     line([leftLick{I(i)}(jl),leftLick{I(i)}(jl)],[i-0.5,i+0.5],'color','b','linewidth',1);%left lick
%                     hold on;
%                 end
%                 for jr=1:length(rightLick{I(i)})
%                     line([rightLick{I(i)}(jr),rightLick{I(i)}(jr)],[i-0.5,i+0.5],'color','r','linewidth',1);%right lick
%                     hold on;
%                 end
%             end
%             set(gca,'ytick',[]);%licking raster have corresponding color plot label, so not need
%             if nStim==size(trialType,2)
%                 set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
%             else
%                 set(gca,'xtick',[]);
%             end
%             xlim(x_lim);
%             if sum(selectedTrialInd)>0
%                 ylim([0.5 sum(selectedTrialInd)+0.5]);
%             end

%             %plot mean trace
%             figure(figMeanTrace);%save mean trace
%             %subplot(size(trialType,2)+3,2*size(trialType,1),2*nResult-1+2*size(trialType,1)*(size(trialType,2)+2));
%             subplot(1,size(trialType,1),nResult);
%             [neuralActivityMean, neuralActivityCI]=fMean_CI(neuralActivity,0.05);
%             %plot CI
%             %{
%             if size(trialType,2)==2%if only two stimuli, plot CI, otherwise, do not plot for the reason of simplicity
%                 xpatch=[1:size(neuralActivity,2), fliplr(1:size(neuralActivity,2))];
%                 ypatch=[neuralActivityCI(1,:),fliplr(neuralActivityCI(2,:))];
%                 p=patch(xpatch,ypatch,color_mean_trace{nStim});%plot confidence interval
%                 p.FaceAlpha=0.1;
%                 p.EdgeColor=color_mean_trace{nStim};%'none';
%                 hold on;
%             end
%             %}
%             plot(1:size(neuralActivity,2),neuralActivityCI(1,:),'Color',color_mean_trace{nStim},'linewidth',0.5);
%             hold on;
%             plot(1:size(neuralActivity,2),neuralActivityCI(2,:),'Color',color_mean_trace{nStim},'linewidth',0.5);
%             hold on;
%             %plot mean
%             curve_meanTrace(nStim)=plot(1:size(neuralActivity,2),neuralActivityMean,'Color',color_mean_trace{nStim},'linewidth',2);
%             hold on;
%             title(titlestr{nResult});
%             y_lim=get(gca,'Ylim');
%             xlim(x_lim);
%             plot([frameNumTime(1)*round(1000/frT),frameNumTime(1)*round(1000/frT)],[y_lim(1),y_lim(2)],'k-');%align to a behavior event
%             xlabel(['time(s) from ',behEventAlign]);
%             if nResult==1
%                 ylabel('\it\DeltaF/F');
%             end
%             set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
%             set(gca,'FontName','Arial','FontSize',14);
%         end
%     end
%     figure(figMeanTrace);
%     subplot(1,size(trialType,1),1);
%     legend(curve_meanTrace(:),num2str(Data_extract.Stimuli(:)),'Location','best');
    saveas(A,[CurrFolder,'\',savefolder,'\',rootname,'ROI-',num2str(roiNo),'-algin to ',behEventAlign,'-sort ',activitySort,'colorplot.jpg'],'jpg');
%     saveas(figMeanTrace,[CurrFolder,'\',savefolder,'\',rootname,'ROI-',num2str(roiNo),'-algin to ',behEventAlign,'-sort ',activitySort,'mean trace.jpg'],'jpg');
%     close all;


%%
% figure(gcf);
% plot((1:size(dff_aligned,2))*frT/1000,dff_aligned(55,:))
