function [figDff] = fLabelROIsDffIncludedTrialRange(dff,savename_figdff,SavedCaTrials,trial2include,varargin)
%UNTITLED label all ROIs' dff with included range
%   Detailed explanation goes here
if isempty(varargin)
    ind_tr_1=1;
else
    ind_tr_1=varargin{1};
end
ntr = length(SavedCaTrials.f_raw); %use all data, or  just set the number of trials to use
% align to behavior event
nFrameEachTrial=cellfun(@(x) size(x,2),SavedCaTrials.f_raw);
ind_1stFrame=zeros(1,length(nFrameEachTrial));
ind_1stFrame(1)=1;
ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
ind_1stFrame=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials
nROI=size(SavedCaTrials.f_raw{1,1},1);
%replot f_raw picture, labeled with used trial range
figDff=figure;
n_col=8;
n_row=ceil(nROI/n_col);
set(gcf,'Position',[0,0,100*n_col,100*n_row]);
for i=1:nROI
    subplot(n_row,n_col,i);
    fLabelDffIncludedTrialRange(dff(i,:),savename_figdff,ind_1stFrame,trial2include);
    title(['ROI-',num2str(i)]);
end
saveas(figDff,savename_figdff,'jpg');
end

%%
function [figDff]=fLabelDffIncludedTrialRange(dff,filepath,ind_1stFrame,trial2include)

plot(dff,'k-');hold on;
yrange=get(gca,'Ylim');
for i=1:size(trial2include,1)
    if strcmp(trial2include,'all')
        tempx=[1,length(dff)];
    else
        tempx=[ind_1stFrame(trial2include(i,1)),ind_1stFrame(trial2include(i,2))];
    end
    xpatch=[tempx,fliplr(tempx)];
    ypatch=[yrange(1),yrange(1),yrange(2),yrange(2)];
    patch(xpatch,ypatch,'b','FaceAlpha',0.5);
end

end