%%
clear
close all
clc
rootpath='E:\FP\pyx262_20200106';
SavingFolder=rootpath;
DatePath=[];
DatePath{1}=rootpath;
files = dir(strcat(rootpath,'\*Virables.mat'));
if length(files)==1
    behaviorFile=files(1).name;
else
    warning('Multiple beh files');
end
%%
% the session to be plotted
for n_path=1:length(DatePath)
    clearvars -except n_path SavingFolder DatePath behaviorFile
    cd(DatePath{n_path});
    %
    FrameInfo = dir('*.log');
%         [TrialCount,TrialStart_FrameCount]=textread(FrameInfo.name,'%d %d','headerlines',12);
    fileID = fopen(FrameInfo.name);
    C=textscan(fileID,'%d %d','HeaderLines',16);
    fclose(fileID);
    TrialCount=C{1,1};
    TrialStart_FrameCount=C{1,2};
    nTrial = length(TrialCount);
%     ImagingSetup
    ImagingSetup=FrameInfo.name(1:end-4);
    disp(ImagingSetup);
    ImagingSetup(ImagingSetup=='_')='-';
    ImagingSetup_L=sprintf('%s_Left', ImagingSetup);
    ImagingSetup_R=sprintf('%s_Right', ImagingSetup);
    ImagingSetup_Ctrl=sprintf('%s_410Ctrl', ImagingSetup);
    
    Session470LEDstartFrame=1;% usually 205
    Session410LEDstartFrame=2;

    fiberstr={'Left','Right'};
    FrameRate=40;
    FrameTime=1000/FrameRate;
    %% compute dff or just load
    if exist('dff_temp.mat','file')
        load('dff_temp');
    else
        dff=fGetFPdff(SavingFolder);
    end
    %% load Beh mat data and extract behavior events
    load(behaviorFile);%load behavior data
    [trialType,behrule] = fGetTrialType( Data_extract);%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
    frT=FrameTime*2;%2 channel, so framerate and frametime should be half
    frameNumTime=[1,2];%from 2s before align point to 5s after align point
    frameNum=double(round(frameNumTime*1000/frT));
    [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime(Data_extract,round(double(TrialStart_FrameCount')/2),frT);
    %% for each dff, plot mean trace, color plot, etc.
    dffstr={'dff baseline correction first','dff motion correction first','dff470','dff410'};
    behEventAlign='first lick';%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
    behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
    masklick='no';%'yes' if mask lickings or 'no' when want to see all activity
    for nfiber=1:2
%         A=figure;
%         set(gcf, 'position', [0 0 1400 600]);
%         figMeanTrace=figure;%plot mean trace
%         set(gcf, 'position', [0 0 1400 300]);

        for ndff=[1] %plot for each dff, see which is better and test whether 410 signals have similar trend
            if strcmp(behEventAlign,'stimOnset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
                [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff{1,ndff}(nfiber,:), behEventFrameIndex,  frameNum );          
            else
                [ dff_aligned, behEvent_aligned, licking_aligned ] = fAlignSigBehEvent( dff{1,ndff}(nfiber,:), behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
            end
            A=figure;
            set(gcf, 'position', [0 0 1400 600]);
            figMeanTrace=figure;%plot mean trace
            set(gcf, 'position', [0 0 1400 300]);
            ColLimit = prctile(dff{1,ndff}(nfiber,:)',98);
            titlestr={'Correct','Error','Miss','Violation'};
            titlestr=strcat(ImagingSetup,'-',fiberstr{nfiber},'-',titlestr);
            color_mean_trace=fMeanTraceColor(behrule,size(trialType,2));
            for nStim=1:size(trialType,2) %for each stimulus%[1,6]%just 2 end trials 
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
                    figure(A);
%                     suptitle(dffstr{ndff});
%                     if nStim == 1 %plot neural activity
%                         subplot(size(trialType,2)+3,2*size(trialType,1),[2*nResult-1,2*nResult-1+2*size(trialType,1)]);
%                     elseif nStim  == size(trialType,2)
%                         subplot(size(trialType,2)+3,2*size(trialType,1),[2*nResult-1+2*size(trialType,1)*size(trialType,2),2*nResult-1+2*size(trialType,1)*size(trialType,2)+2*size(trialType,1)]);
%                     else
%                         subplot(size(trialType,2)+3,2*size(trialType,1),2*nResult-1+2*size(trialType,1)*nStim);
%                     end
                    subplot(size(trialType,2),2*size(trialType,1),2*nResult-1+2*size(trialType,1)*(nStim-1));
                    hold on;
                    imagesc(neuralActivity(I,:));
                    %x_lim=get(gca,'Xlim');
                    x_lim=[0,sum(frameNum)];
                    %plot behavior event
                    plot(sot,[1:sum(selectedTrialInd)],'w.','markersize',8);
                    hold on;
                    plot(go,[1:sum(selectedTrialInd)],'w.','markersize',8);
                    hold on;
%                     plot(flt,[1:sum(selectedTrialInd)],'k.','markersize',8);
%                     hold on;
%                     plot(flt_l,[1:sum(selectedTrialInd)],'m.','markersize',8);
%                     hold on;
%                     plot(flt_r,[1:sum(selectedTrialInd)],'r.','markersize',8);
%                     hold on;
                    plot(rwt,[1:sum(selectedTrialInd)],'c.','markersize',8);
                    hold on;
                    set(gca,'clim',[0 ColLimit]);
                    set(gca,'ytick',sum(selectedTrialInd),'yticklabel',sum(selectedTrialInd),'ydir','normal');
                    if nStim==size(trialType,2)
                        set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
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
                        ylabel([num2str(Data_extract.Stimuli(nStim))]);% ' clicks/s'
                    end
                    xlim(x_lim);
%                     if nStim == 1 %plot licking raster
%                         subplot(size(trialType,2)+3,2*size(trialType,1),[2*nResult,2*nResult+2*size(trialType,1)]);
%                     elseif nStim  == size(trialType,2)
%                         subplot(size(trialType,2)+3,2*size(trialType,1),[2*nResult+2*size(trialType,1)*size(trialType,2),2*nResult+2*size(trialType,1)*size(trialType,2)+2*size(trialType,1)]);
%                     else
%                         subplot(size(trialType,2)+3,2*size(trialType,1),2*nResult+2*size(trialType,1)*nStim);
%                     end
                    subplot(size(trialType,2),2*size(trialType,1),2*nResult+2*size(trialType,1)*(nStim-1));
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
                        set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
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
                    subplot(1,size(trialType,1),nResult);
                    [neuralActivityMean, neuralActivityCI]=fMean_CI(neuralActivity,0.05);
                    %shadow as ci
                    %
                    xpatch=[1:size(neuralActivity,2), fliplr(1:size(neuralActivity,2))];
                    ypatch=[neuralActivityCI(1,:),fliplr(neuralActivityCI(2,:))];
                    p=patch(xpatch,ypatch,color_mean_trace{nStim});%plot confidence interval
                    p.FaceAlpha=0.1;
                    p.EdgeColor=color_mean_trace{nStim};%'none';
                    hold on;
                    %}
%                     for temp=1:size(neuralActivity,1)
%                         plot(1:size(neuralActivity,2),neuralActivity(temp,:),'Color',color_mean_trace{nStim},'linewidth',0.5);
%                     end
%                     plot(1:size(neuralActivity,2),neuralActivityCI(1,:),'Color',color_mean_trace{nStim},'linewidth',0.5);
%                     hold on;
%                     plot(1:size(neuralActivity,2),neuralActivityCI(2,:),'Color',color_mean_trace{nStim},'linewidth',0.5);
%                     hold on;
                    curve_meanTrace(nStim)=plot(1:size(neuralActivity,2),neuralActivityMean,'Color',color_mean_trace{nStim},'linewidth',2);
                    hold on;
                    title(titlestr{nResult});
                    y_lim=get(gca,'Ylim');
                    xlim(x_lim);
                    plot([frameNumTime(1)*round(1000/frT),frameNumTime(1)*round(1000/frT)],[y_lim(1),y_lim(2)],'k-');%align to a behavior event
                    xlabel(['time(s) from ',behEventAlign]);
                    if nResult==1
                        ylabel('\it\DeltaF/F');
                    end
                    set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
                    set(gca,'FontName','Arial','FontSize',14);                    
                end
            end
            figure(figMeanTrace);
            subplot(1,size(trialType,1),1);
            legend(curve_meanTrace(:),num2str(Data_extract.Stimuli(:)),'Location','best');
            if strcmp(masklick,'yes')
                suptitle([dffstr{ndff},'-masklick']);
            else
                suptitle(dffstr{ndff});
            end
            figure(A);
            suptitle(dffstr{ndff});
            saveas(A,[ImagingSetup,'-',fiberstr{nfiber},dffstr{ndff},'-algin to ',behEventAlign,'-sort ',behEventSort,'-masklick',masklick,'colorplot.fig'],'fig');
            saveas(figMeanTrace,[ImagingSetup,'-',fiberstr{nfiber},dffstr{ndff},'-algin to ',behEventAlign,'-sort ',behEventSort,'-masklick',masklick,'mean trace.fig'],'fig');
            close all;
        end
    end
    
end
