function [curve_meanTrace] = fPlotDffPSTH_ASession(dff,ind_tr_1,Data_extract,SavedCaTrials,frameNumTime,behEventAlign,masklick,i_selectivity,trial2include,trial2exclude,savename_figdff,title_fig)
%FPLOTDFFPSTHASESSION plot dff PSTH of one session,only plot correct trials
%   Detailed explanation goes here

ntr = length(SavedCaTrials.f_raw); %use all data, or  just set the number of trials to use
frT = SavedCaTrials.FrameTime;
% align to behavior event
nFrameEachTrial=cellfun(@(x) size(x,2),SavedCaTrials.f_raw);
ind_1stFrame=zeros(1,length(nFrameEachTrial));
ind_1stFrame(1)=1;
ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
ind_1stFrame=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials
if isnan(trial2include)
    indTrial2include=fExcludeTrials(trial2exclude,ind_1stFrame,'logical');
    trial2include=[1,length(ind_1stFrame)];
else
    indTrial2include=fIncludeTrials(trial2include,ind_1stFrame,'logical');
end
frameNum=double(round(frameNumTime*1000/frT));
[behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( Data_extract, ind_1stFrame, SavedCaTrials.FrameTime ,ind_tr_1);%get behavior event time


% %replot f_raw picture, labeled with used trial range
% fLabelDffIncludedTrialRange(dff,savename_figdff,ind_1stFrame,trial2include)
%   align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
if strcmp(behEventAlign,'delay onset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
    [ dff_aligned, ~,~ ] = fAlignDelaySigal( dff, behEventFrameIndex,  frameNum );
else
    [ dff_aligned, ~, ~ ] = fAlignSigBehEvent( dff, behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
end

[trialType,rule] = fGetTrialType( Data_extract,[],i_selectivity,'matrix','left','divideCorErr');%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials

%%
%plot color plot,4 column(correct/error/miss/violation), 7 row(for each
%stimlus)+ mean trace

titlestr=title_fig;
if size(trialType,2)==6
    color_mean_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};
elseif size(trialType,2)==2
    color_mean_trace={[0 0 1],[1 0 0]};
end
%     if strcmp(rule,'low click rate-right')
%         color_mean_trace=fliplr(color_mean_trace);
%     end
for nStim=1:size(trialType,2) %for each stimulus
    for nResult=1%:size(trialType,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
        selectedTrialInd=trialType(nResult,nStim,:);
        selectedTrialInd=logical(squeeze(selectedTrialInd))';
        selectedTrialInd=logical(selectedTrialInd.*indTrial2include');
        neuralActivity=dff_aligned(selectedTrialInd,:);

        %plot mean trace
        %[neuralActivityMean, neuralActivityCI]=fMean_CI(neuralActivity,0.05);
        ts=double((-frameNum(1):1:frameNum(2))*frT/1000);
        curve_meanTrace(nStim)=fPlotMean_SE(ts,neuralActivity,color_mean_trace{nStim});
        title(titlestr);
    end
end

%% label and general plots
% find common ylim
for  nResult=1%:size(trialType,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
    
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
end
for  nResult=1%:size(trialType,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
    xlim([ts(1),ts(end)]);
    plot([0,0],[y_lim(1),y_lim(2)],'k-');%align to a behavior event
    if strcmp(behEventAlign,'stim onset')
        plot(0.5*ones(1,2),[y_lim(1),y_lim(2)],'k-');%align to a behavior event
    elseif strcmp(behEventAlign,'delay onset')
        plot(-0.5*ones(1,2),[y_lim(1),y_lim(2)],'k-');%align to a behavior event
    end
    

    %set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
    set(gca,'Ylim',y_lim,'xtick',[-floor(frameNumTime(1)):1:frameNumTime(2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
    %         set(gca,'FontName','Arial','FontSize',14);
end

% label and general plots
% find common ylim
clear y_lim;
for  nResult=1%:size(trialType,1) %4 column(correct/error/miss/violation),companied with 4 lick raster

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
end

if flagclearyrange==1
    clear y_lim;
end
end

%%
function []=fLabelDffIncludedTrialRange(dff,filepath,ind_1stFrame,trial2include)
% if exist(filepath,'file')
%     return;
% else
    figDff=figure;
    plot(dff,'k-');hold on;
    yrange=get(gca,'Ylim');
    for i=1:size(trial2include,1)
        tempx=[ind_1stFrame(trial2include(i,1)),ind_1stFrame(trial2include(i,2))];
        xpatch=[tempx,fliplr(tempx)];
        ypatch=[yrange(1),yrange(1),yrange(2),yrange(2)];
        patch(xpatch,ypatch,'b','FaceAlpha',0.5);
    end
    saveas(figDff,filepath,'jpg');
% end
end

