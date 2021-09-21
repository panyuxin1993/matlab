function [I] = fRasterDff(neuralActivity,behEvent_aligned,selectedTrialInd,behEventSort,ColLimit)
%FRASTERDFF plot a raster plot of neural activity with behavior event markers
%   (each time only one block)
%Input-
%   neuralActivity=m-by-n matrix, m trials, each n frame number
%   behEvent_aligned= struct, behavior event frame index(including
%   stimOnset, answer time, first lick and so on),
%   selectedTrialInd = 3/4-d matrix,1d(cor/err/miss/vio or do/miss/vio),
%   2d(stimuli/sensory/choice),3d(trial numbers),4d(optional, divide opto
%   or not)
%   behEventSort=char array, 'first lick'(default)|'reward'|'go
%   cue',indicating which event to be aligned
%   ColLimit- color of scale bar
%Output-
%   I-useful when licking raster was needed, and this parameter help to
%   keep each line corresponding to lickings

sot=behEvent_aligned.stimOnset(selectedTrialInd);% stim onset, white
flt=behEvent_aligned.lickFirst(selectedTrialInd);% first lick, black
flt_l=behEvent_aligned.lickFirst_left(selectedTrialInd);%left lick, magenta dots(ipsi lick)
flt_r=behEvent_aligned.lickFirst_right(selectedTrialInd);%right lick, red dots(contra lick)
rwt=behEvent_aligned.rewTime(selectedTrialInd);%reward, cyan dots
go=behEvent_aligned.go(selectedTrialInd);%go cue,white


%plot raster
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
imagesc(neuralActivity(I,:));
hold on;

%plot behavior event
plot(sot,1:sum(selectedTrialInd),'w.','markersize',8);
hold on;
plot(go,1:sum(selectedTrialInd),'w.','markersize',8);
plot(rwt,1:sum(selectedTrialInd),'c.','markersize',8);
hold on;
set(gca,'clim',ColLimit);
set(gca,'ytick',sum(selectedTrialInd),'yticklabel',sum(selectedTrialInd),'ydir','normal');
set(gca,'xtick',[]);
ColBar = colorbar;
set(ColBar,'position',[0.91 0.2400 0.02 0.1445]);
if sum(selectedTrialInd)>0
    ylim([0.5 sum(selectedTrialInd)+0.5]);
end
end

