function [ fig ] = fPlotFPExample( fig,position )
%FPLOTFPEXAMPLE plot example site from FP recording data
% color plot and meantrace are in same fig
SavingFolder='D:\FP\pyx191_20190602';
rootpath='D:\FP\pyx191_20190602\';
behaviorFile='D:\FP\pyx191_20190602\2019_06_02_pyx191_FP_Virables.mat';
FrameInfo = dir([rootpath,'*.log']);
fileID = fopen([rootpath,FrameInfo.name]);
C=textscan(fileID,'%d %d','HeaderLines',16);
%     ImagingSetup
ImagingSetup=FrameInfo.name(1:end-4);
disp(ImagingSetup);
ImagingSetup(ImagingSetup=='_')='-';
fclose(fileID);
TrialCount=C{1,1};
TrialStart_FrameCount=C{1,2};
fiberstr={'Left','Right'};
%% load Beh mat data and extract behavior events
load([rootpath,'dff_temp']);
load(behaviorFile);%load behavior data
trialType = fGetTrialType( Data_extract);%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
trialType(3:4,:,:)=[];%do not plot vio and miss
FrameRate=40;
FrameTime=1000/FrameRate;
frT=FrameTime*2;%2 channel, so framerate and frametime should be half
frameNumTime=[1,2];%from 2s before align point to 5s after align point
frameNum=double(round(frameNumTime*1000/frT));
[behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime(Data_extract,round(double(TrialStart_FrameCount')/2),frT);
%% compare different dff('dff baseline correction first','dff motion correction first','dff470','dff410')
%for each dff, plot mean trace, color plot, etc.
dffstr={'dff baseline correction first','dff motion correction first','dff470','dff410'};
behEventAlign='first lick';%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
masklick='no';%'yes' if mask lickings or 'no' when want to see all activity
for nfiber=1:1
    %         A=figure;
    %         set(gcf, 'position', [0 0 1400 600]);
    %         figMeanTrace=figure;%plot mean trace
    %         set(gcf, 'position', [0 0 1400 300]);
    
    for ndff=1:1 %plot for each dff, see which is better and test whether 410 signals have similar trend
        if strcmp(behEventAlign,'stimOnset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
            [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff{1,ndff}(nfiber,:), behEventFrameIndex,  frameNum );
        else
            [ dff_aligned, behEvent_aligned, licking_aligned ] = fAlignSigBehEvent( dff{1,ndff}(nfiber,:), behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
        end
        figure(fig);
        set(fig, 'PaperPosition', position);
        ColLimit = prctile(dff{1,ndff}(nfiber,:)',98);
        titlestr={'Correct','Error','Miss','Violation'};
%         titlestr=strcat(ImagingSetup,'-',fiberstr{nfiber},'-',titlestr);
        if contains(SavingFolder,'pyx159')
            color_mean_trace={[0 0 1],[1 0 0]};
        elseif contains(SavingFolder,'pyx214')
            color_mean_trace={[1 0 0],[0 0 1]};
        elseif contains(SavingFolder,'pyx172')
            color_mean_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};
        elseif contains(SavingFolder,'pyx161')||contains(SavingFolder,'pyx19')
            %                 color_mean_trace={[1 0 0],[0 0 1]};
            color_mean_trace={[1 0 0],[0.8 0 0.2],[0.6 0 0.4],[0.4 0 0.6],[0.2 0 0.8],[0 0 1]};
        end
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
                figure(fig);
                
                subplot(size(trialType,2)+1,2*size(trialType,1),2*nResult-1+2*size(trialType,1)*(nStim-1));
                hold on;
                imagesc(neuralActivity(I,:));
                %x_lim=get(gca,'Xlim');
                x_lim=[0,sum(frameNum)];
                %plot behavior event
                plot(sot,[1:sum(selectedTrialInd)],'w.','markersize',8);
                hold on;
                plot(go,[1:sum(selectedTrialInd)],'w.','markersize',8);
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
                    ylabel([num2str(Data_extract.Stimuli(nStim))]);% ' clicks/s'
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
                subplot(size(trialType,2)+1,2*size(trialType,1),2*nResult-1+2*size(trialType,1)*size(trialType,2));
                %subplot(1,size(trialType,1),nResult);
                [neuralActivityMean, neuralActivityCI]=fMean_CI(neuralActivity,0.05);
%                 xpatch=[1:size(neuralActivity,2), fliplr(1:size(neuralActivity,2))];
%                 ypatch=[neuralActivityCI(1,:),fliplr(neuralActivityCI(2,:))];
%                 p=patch(xpatch,ypatch,color_mean_trace{nStim});%plot confidence interval
%                 p.FaceAlpha=0.1;
%                 p.EdgeColor=color_mean_trace{nStim};%'none';
%                 hold on;
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
                set(gca,'ylim',[0,0.6]);
%                 set(gca,'FontName','Arial','FontSize',14);
                box off;
            end
        end
        subplot(size(trialType,2)+1,2*size(trialType,1),2*size(trialType,1)*(size(trialType,2)+1)-1);
        leg=legend(curve_meanTrace(:),num2str(Data_extract.Stimuli(:)));
        set(leg,'position',[0.91 0.1 0.02 0.08]);
        legend boxoff;

%         if strcmp(masklick,'yes')
%             suptitle([dffstr{ndff},'-masklick']);
%         else
%             suptitle(dffstr{ndff});
%         end
%         suptitle(dffstr{ndff});

    end
end




end

