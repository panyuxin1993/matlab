function [figRaster,figMeanTrace,figBehavior] = fPlotRasterPSTHoptoVSctrl(Data_extract,behEventAlign, masklick, behEventSort,dff,IDstr,frameNumTime,ind_1stFrame,frameRate,ind_tr_1,f,varargin)
%FPLOTPSTHOPTOVSCTRL plot raster plot of activities and mean traces,
%usually 2 stimuli, opto/ctrl condition, and fixed delay length
%Data_extract- behavior data extracted as a mat file, including behavior
%   event of each trial from a session
%behEventAlign- the behavior event to align, {'stimOnset','go cue',...
%   'first lick','first left lick','first right lick', 'answer','reward'};
%behEventSort- which event to be sort, {'first lick','reward','go cue'};
%masklick- whether mask activities after licks, {'yes','no'};
%IDstr- the ID of an ROI(2P) or fiber(FP)
%dff- neural activities,in matrix form
%frameNumTime- [time before aligned time point, time after that time point]
%ind_1stFrame- the frame index of each trial start
%frameRate- sample rate of 2p and FP (Hz)
%ind_tr_1-from which trial to start(usually the 1st trial)
%f- method for trial type 
%varargin- can be yrange
%Output- figure handle of 2 figures.

if ~isempty(varargin)
    yrange=varargin{1};
end
frameNum=double(round(frameNumTime*frameRate));
[behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( Data_extract, ind_1stFrame, 1000/frameRate ,ind_tr_1);%get behavior event time

if strcmp(behEventAlign,'delay onset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
    [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff, behEventFrameIndex,  frameNum );
else
    [ dff_aligned, behEvent_aligned, licking_aligned ] = fAlignSigBehEvent( dff, behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
end
[trialType,rule] = fGetTrialType( Data_extract,[],f,'matrix','left','divideCorErr','divideOpto');%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials

%%
%plot color plot,4 column(correct/error/miss/violation), 7 row(for each
%stimlus)+ mean trace
figRaster=figure;
set(gcf, 'position', [0 0 1400 600]);
figMeanTrace=figure;%plot mean trace
set(gcf, 'position', [0 600 1400 300]);
figBehavior=figure;%plot behavior
set(gcf, 'position', [1400 0 300 300]);
ColLimit = prctile(dff,90);
titlestr={'Correct','Error','Miss','Violation'};
titlestr=strcat(IDstr,titlestr);
titlestropto={'Ctrl','Opto'};
nTrial=zeros(size(trialType,1),size(trialType,2),size(trialType,4));% 3d-opto/non-opto,2d-stimuli(usually 2 or 4)
pRightChoice=zeros(size(trialType,2),size(trialType,4));% 2d-opto/non-opto,1d-stimuli(usually 2 or 4)
if size(trialType,2)==6
    color_mean_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};
elseif size(trialType,2)==2
    color_mean_trace={[0 0 1],[0,0.5,1];[1,0,0],[1,0.5,0]};
end
if strcmp(rule,'low click rate-right')
    color_mean_trace=flipud(color_mean_trace);
end
%color_adjust=[1,0.5;1,0.5;1,0.5];
for nResult=1:size(trialType,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
    n_curveM=1;
    for nStim=1:size(trialType,2) %for each stimulus
        for nOpto=1:size(trialType,4)
            selectedTrialInd=trialType(nResult,nStim,:,nOpto);
            selectedTrialInd=logical(squeeze(selectedTrialInd))';
            sot=behEvent_aligned.stimOnset(selectedTrialInd);% stim onset, white
            flt=behEvent_aligned.lickFirst(selectedTrialInd);% first lick, black
            flt_l=behEvent_aligned.lickFirst_left(selectedTrialInd);%left lick, magenta dots
            flt_r=behEvent_aligned.lickFirst_right(selectedTrialInd);%right lick, red dots
            rwt=behEvent_aligned.rewTime(selectedTrialInd);%reward, cyan dots
            go=behEvent_aligned.go(selectedTrialInd);%go cue,white
            oon=behEvent_aligned.optoOnset(selectedTrialInd);%opto onset, a bar out of the block
            ooff=behEvent_aligned.optoOffset(selectedTrialInd);
            neuralActivity=dff_aligned(selectedTrialInd,:);
            temp=find(selectedTrialInd);
            nTrial(nResult,nStim,nOpto)=length(temp);
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
            figure(figRaster);
            if nOpto==1%ctrl
                subplot(size(trialType,2),2*size(trialType,1),2*nResult-1+2*size(trialType,1)*(nStim-1));
            elseif nOpto==2%opto
                subplot(size(trialType,2),2*size(trialType,1),2*nResult+2*size(trialType,1)*(nStim-1));
            end
            %% try to remove opto artifacts
%             if nOpto>1 %opto condition
%                 bin2calculateArtifact=0.5*frameRate;% s
%                 artifactAmp=mean(mean(neuralActivity(:,oon+1:oon+bin2calculateArtifact)))-mean(mean(neuralActivity(:,oon-bin2calculateArtifact:oon-1)));
%                 neuralActivity(:,oon:ooff)=neuralActivity(:,oon:ooff)-artifactAmp;
%             end
            %% continue to plot raster and mean trace
            imagesc(neuralActivity(I,:));
            hold on;
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
                set(gca,'xtick',[1:round(frameRate):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
                xlabel(['Time (s) from ',behEventAlign],'FontName','Arial','FontSize',14);
            else
                set(gca,'xtick',[]);
            end
            ColBar = colorbar;
            set(ColBar,'position',[0.91 0.1100 0.02 0.1445]);
            %                 if sum(selectedTrialInd)>0
            %                     ylim([0.5 sum(selectedTrialInd)+0.5]);
            %                 end
            if nStim == 1
                title([titlestr{nResult},titlestropto{nOpto}]);
            end
            if nResult==1
                ylabel([num2str(Data_extract.Stimuli(nStim))]); %' clicks/s'
            end
            if nOpto==2 && strcmp(behEventAlign,'delay onset')%opto and in one block they aligned well
                box off;
                set(gca,'Ylim',[0 sum(selectedTrialInd)+1]);
                y_opto=sum(selectedTrialInd)+1;
                plot([mean(oon),mean(ooff)],[y_opto,y_opto],'r-','LineWidth',1);
                hold on;
            end
            xlim(x_lim);
            
%             %plot licking raster
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
            
            %plot mean trace
            figure(figMeanTrace);%save mean trace
            %subplot(size(trialType,2)+3,2*size(trialType,1),2*nResult-1+2*size(trialType,1)*(size(trialType,2)+2));
            subplot(1,size(trialType,1),nResult);
%             ts=-frameNumTime(1):1/frameRate:frameNumTime(2);
            ts=1:size(neuralActivity,2);
            curve_meanTrace(n_curveM)=fPlotMean_SE(ts,neuralActivity,color_mean_trace{nStim,nOpto});
            hold on;
            str_lgd{n_curveM}=[num2str(Data_extract.Stimuli(nStim)),titlestropto{nOpto}];
            n_curveM=n_curveM+1;
            title(titlestr{nResult});
            %plot opto-stimuli duration
            if exist('yrange','var')
                y_lim=yrange;
            else
                if exist('y_lim','var')
                    y_lim_temp=get(gca,'Ylim');
                    y_lim=[min(y_lim(1),y_lim_temp(1)),max(y_lim(2),y_lim_temp(2))];
                else
                    y_lim=get(gca,'Ylim');
                end
            end           
        end
    end
    if nOpto==2 && strcmp(behEventAlign,'delay onset') && nStim==size(trialType,2) && nResult ~=size(trialType,1) %not plot violation
        xpatch_opto=[mean(oon),mean(ooff),mean(ooff),mean(oon)];
        ypatch_opto=[y_lim(1),y_lim(1),y_lim(end),y_lim(end)];
        p_opto=patch(xpatch_opto,ypatch_opto,'b');
        p_opto.FaceAlpha=0.1;
    end
    xlim(x_lim);
    plot((frameNumTime(1)*round(frameRate)+1)*ones(1,2),[y_lim(1),y_lim(2)],'k-');%align to a behavior event
    ylim(y_lim);
    clearvars y_lim;
end
for nResult=1:size(trialType,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
    
    xlabel(['Time (s) from ',behEventAlign]);
    if nResult==1
        ylabel('\it\DeltaF/F');
    end
    figure(figMeanTrace);%save mean trace
    subplot(1,size(trialType,1),nResult);
    set(gca,'xtick',[1:round(frameRate):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
    set(gca,'FontName','Arial','FontSize',14);
end

figure(figMeanTrace);
subplot(1,size(trialType,1),1);
legend(curve_meanTrace(:),str_lgd,'Location','best');

for nOpto=1:size(trialType,4)
    for nStim=1:size(trialType,2)/2 %for each stimulus, first half
        switch rule
            case 'low click rate-left'
                pRightChoice(nStim,nOpto)=nTrial(2,nStim,nOpto)/(nTrial(1,nStim,nOpto)+nTrial(2,nStim,nOpto));
            case 'low click rate-right'
                pRightChoice(nStim,nOpto)=nTrial(1,nStim,nOpto)/(nTrial(1,nStim,nOpto)+nTrial(2,nStim,nOpto));
        end
    end
    for nStim=size(trialType,2)/2+1:size(trialType,2) %for each stimulus, first half
        switch rule
            case 'low click rate-left'
                pRightChoice(nStim,nOpto)=nTrial(1,nStim,nOpto)/(nTrial(1,nStim,nOpto)+nTrial(2,nStim,nOpto));
            case 'low click rate-right'
                pRightChoice(nStim,nOpto)=nTrial(2,nStim,nOpto)/(nTrial(1,nStim,nOpto)+nTrial(2,nStim,nOpto));
        end
    end
end
figure(figBehavior);
curve(1)=scatter(1:size(pRightChoice,2),pRightChoice(:,1),'k');hold on;%control
curve(2)=scatter(1:size(pRightChoice,2),pRightChoice(:,2),'r');%opto
legend(curve(:),'ctrl','opto');
set(gca,'Xlim',[0,size(pRightChoice,2)+1],'Ylim',[0,1]);
ylabel('P(Right Choice)');
xlabel('stimuli');
set(gca,'FontSize',14);
%close all;
end

