%continuing sampling, so link each trials one by one to get a 'long trial'
% load('F:\pyx083-20180522\im_data_reg_cpu\result_save\CaTrialsSIM_pyx083_20180522_920nm_power50_zoom4x_dftReg_.mat');
% load('D:\xulab\behavior\pyx083\Data_Virables\2018_05_22_pyx083-imaging_Virables.mat');
clear;
close all;
% sessionNamePool={'pyx290_20200528'};%roiNo = 15
% yrange=[-0.5,4];
sessionNamePool={'pyx224_20190729'};%roiNo = 23,4
% yrange=[-2,6];%roiNo =4
behEventAlignPool={'delay onset','go cue','first lick'};
masklickPool={'yes','no','no'};

%plotting settings
frameNumTime=[1,1.5];%from 1s before align point to 1.5s after align point

%plot color plot+ mean trace
A=figure;
set(gcf, 'position', [0 0 300*length(behEventAlignPool)+200 300]);
figMeanTrace=figure;%plot mean trace
set(gcf, 'position', [0 300 300*length(behEventAlignPool)+200 300]);
for i_align=1:length(behEventAlignPool)
    %decide some global variable
    behEventAlign=behEventAlignPool{i_align};%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
    behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
    masklick=masklickPool{i_align};
    
    i_selectivity=4;%*********variable**************
    selectivitystr={'stimuli','sensory difficulty','sensory','choice'};%sensory means grouping difficulties;
    trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
    trialTypeStr=selectivitystr{i_selectivity};
    
    for i_session=1:length(sessionNamePool)
        sessionName=sessionNamePool{i_session};
        datapath=['H:\2P\',sessionName];
        clear CurrFolder;
        cd([datapath,'\im_data_reg\result_save']);
        CurrFolder=pwd;
        savefolder=sessionName;%'1-200trials';%'segNP';
        % savefolder='F:\2P\example\';
        
        % load([CurrFolder,'\',sessionName, '_imaging_Virables.mat']);%load behavior data
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
        
        if ~exist(savefolder)
            mkdir(savefolder);
        end
        
        if exist('dff.mat','file')
            load([CurrFolder,'\','dff.mat']);%load dff
        else
            dff=fnx_getDff(CurrFolder,sessionName,'save figure');
        end
        
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
        [trialType,rule] = fGetTrialType( Data_extract,[],i_selectivity,'matrix','left','divideCorErr');%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
        
        for roiNo = 33%1:size(SavedCaTrials.f_raw{1},1) %SavedCaTrials.nROIs may be not true
            %   align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
            if strcmp(behEventAlign,'delay onset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
                [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff(roiNo,:), behEventFrameIndex,  frameNum );
            else
                [ dff_aligned, behEvent_aligned, licking_aligned ] = fAlignSigBehEvent( dff(roiNo,:), behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
            end
            
            ColLimit = prctile(dff(roiNo,:),90);
            titlestr={'Correct','Error','Miss','Violation'};
            titlestr=strcat('ROI-',num2str(roiNo),titlestr);
            if size(trialType,2)==6
                color_mean_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};
                color_mean_trace_err={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};
            elseif size(trialType,2)==2
                color_mean_trace={[0 0 1],[1 0 0]};
                if strcmp(trialTypeStr,'choice')
                    color_mean_trace_err={[1 0.5 0.5],[0.5 0.5 1]};
                elseif strcmp(trialTypeStr,'sensory')
                    color_mean_trace_err={[0.5 0.5 1],[1 0.5 0.5]};
                end
            end
            %     if strcmp(rule,'low click rate-right')
            %         color_mean_trace=fliplr(color_mean_trace);
            %     end
            for nStim=1:size(trialType,2) %for each stimulus
                nResultPlot=size(trialType,1)-2; %2 column(correct/error/miss/violation),companied with 4 lick raster,only plot cor and err
                for  nResult=1:nResultPlot
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
                    figure(A);
                    %             if nStim == 1 %plot neural activity; here end trials have more trials so wider
                    %                 subplot(size(trialType,2)+3,2*size(trialType,1),[2*nResult-1,2*nResult-1+2*size(trialType,1)]);
                    %             elseif nStim  == size(trialType,2)
                    %                 subplot(size(trialType,2)+3,2*size(trialType,1),[2*nResult-1+2*size(trialType,1)*size(trialType,2),2*nResult-1+2*size(trialType,1)*size(trialType,2)+2*size(trialType,1)]);
                    %             else
                    %                 subplot(size(trialType,2)+3,2*size(trialType,1),2*nResult-1+2*size(trialType,1)*nStim);
                    %             end
                    subplot(size(trialType,2),2*nResultPlot*length(behEventAlignPool),2*nResult+4*i_align-5+(nStim-1)*2*nResultPlot*length(behEventAlignPool));
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
                    if nStim==size(trialType,2)
                        set(gca,'xtick',[round(1000/frT*(frameNumTime(1)-floor(frameNumTime(1))))+1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
                        xlabel(['time(s) from ',behEventAlign],'FontName','Arial','FontSize',14);
                    else
                        set(gca,'xtick',[]);
                    end
                    ColBar = colorbar;
                    set(ColBar,'position',[0.91 0.1100 0.02 0.1445]);
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
                    
                    %             if nStim == 1 %plot licking raster
                    %                 subplot(size(trialType,2)+3,2*size(trialType,1),[2*nResult,2*nResult+2*size(trialType,1)]);
                    %             elseif nStim  == size(trialType,2)
                    %                 subplot(size(trialType,2)+3,2*size(trialType,1),[2*nResult+2*size(trialType,1)*size(trialType,2),2*nResult+2*size(trialType,1)*size(trialType,2)+2*size(trialType,1)]);
                    %             else
                    %                 subplot(size(trialType,2)+3,2*size(trialType,1),2*nResult+2*size(trialType,1)*nStim);
                    %             end
                    subplot(size(trialType,2),2*nResultPlot*length(behEventAlignPool),2*nResult+4*i_align-4+(nStim-1)*2*nResultPlot*length(behEventAlignPool));
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
                    if nStim==size(trialType,2)
                        %set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
                        %set(gca,'xtick',[-round(1000/frT*frameNumTime(1))+1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
                        set(gca,'xtick',[round(1000/frT*(frameNumTime(1)-floor(frameNumTime(1))))+1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
                    else
                        set(gca,'xtick',[]);
                    end
                    xlim(x_lim);
                    if sum(selectedTrialInd)>0
                        ylim([0.5 sum(selectedTrialInd)+0.5]);
                    end
                    
                    %plot mean trace
                    figure(figMeanTrace);%save mean trace
                    %subplot(size(trialType,2)+3,2*size(trialType,1),2*nResult-1+2*size(trialType,1)*(size(trialType,2)+2));
                    subplot(1,length(behEventAlignPool),i_align);
                    %             [neuralActivityMean, neuralActivityCI]=fMean_CI(neuralActivity,0.05);
                    ts=double((-frameNum(1):1:frameNum(2))*frT/1000);
                    hold on;
                    if nResult==1
                        curve_meanTrace(nStim)=fPlotMean_SE(ts,neuralActivity,color_mean_trace{nStim});
                    elseif nResult==2
                        curve_meanTrace(nStim+size(trialType,2))=fPlotMean_SE(ts,neuralActivity,color_mean_trace_err{nStim});
                    end
                end
            end
        end  
    end
    % label and general plots
    subplot(1,length(behEventAlignPool),i_align);
    flagclearyrange=0;
    if ~exist('yrange','var')
        flagclearyrange=1;
        if ~exist('y_lim','var')
            y_lim=get(gca,'Ylim');
        else
            y_lim2=get(gca,'Ylim');
            y_lim=[min(y_lim(1),y_lim2(1)),max(y_lim(2),y_lim2(2))];
        end
    else
        y_lim=yrange;
    end
    xlim([ts(1),ts(end)]);
    plot([0,0],[y_lim(1),y_lim(2)],'k-');%align to a behavior event
    if strcmp(behEventAlign,'stim onset')
        plot(0.5*ones(1,2),[y_lim(1),y_lim(2)],'k-');%align to a behavior event
    elseif strcmp(behEventAlign,'delay onset')
        plot(-0.5*ones(1,2),[y_lim(1),y_lim(2)],'k-');%align to a behavior event
    end
    xlabel(['Time (s) from ',behEventAlign]);
    if nResult==1
        ylabel('\it\DeltaF/F');
    end
    %set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
    set(gca,'Ylim',y_lim,'xtick',[-floor(frameNumTime(1)):1:frameNumTime(2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
    set(gca,'FontName','Arial','FontSize',14);
end


% figure(figMeanTrace);
% subplot(1,length(behEventAlignPool),1);
% if contains(trialTypeStr,'stimuli')
%     %         h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','','contra hard','','contra easy'},'Location','best');
%     h=legend(curve_meanTrace(:),[strcat(num2str(Data_extract.Stimuli(:)),'cor'),strcat(num2str(Data_extract.Stimuli(:)),'err')],'Location','best');
% elseif contains(trialTypeStr,'difficulty')
%     h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','contra hard','contra easy'},'Location','best');
% elseif contains(trialTypeStr,'choice')
%     h=legend(curve_meanTrace(:),{'ipsi choice correct','contra choice correct','contra choice error','ipsi choice error'},'Location','best');
% elseif contains(trialTypeStr,'sensory')
%     h=legend(curve_meanTrace(:),{'low click rate correct','high click rate correct','low click rate error','high click rate error'},'Location','best');
% end
% set(h,'box','off');

% saveas(A,[CurrFolder,'\',savefolder,'\',sessionName,'ROI-',num2str(roiNo),'-algin to ',behEventAlign,'-sort ',behEventSort,'colorplot.jpg'],'jpg');
% saveas(figMeanTrace,[CurrFolder,'\',savefolder,'\',sessionName,'ROI-',num2str(roiNo),'-algin to ',behEventAlign,'-sort ',behEventSort,'mean trace.jpg'],'jpg');

% set(A,'PaperPosition',[0,0,4,2]);
% saveas(A,['F:\2P\example\',sessionName,'ROI-',num2str(roiNo),'-algin to ',behEventAlign,'-sort ',behEventSort,'colorplot.pdf'],'pdf');
set(figMeanTrace,'PaperPosition',[0,0,6,2]);
saveas(figMeanTrace,['H:\2P\example\',sessionName,'ROI-',num2str(roiNo),'mean trace.pdf'],'pdf');
close all;
if flagclearyrange==1
    clear y_lim;
end
