%compare the probability/duration of significant activitiy of different
%behavioral epoch
%% load data
behEventAlign='delay onset';
masklick='no';
behEventSort='go cue';
bodypartsPool={'Tongue'};
coordinatesPool={'x'};
titlestrPool={'Tongue (x)'};
yrange={[-40,40]};
selectivity='choice';
treshold4likelihood=0.1;
pSigTtest=0.01;
nonOverlappingBin=1;
baseline='';%'LickPort';

[num,txt,raw] =xlsread('C:\Users\PYX\Documents\DataSummary\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
T=cell2table(raw(2:end,1:15));
T.Properties.VariableNames=strrep(raw(1,1:15),' ','_');%table variable name can't have ' ',so replace them
%{
ind_session=logical(contains(T.ROI_type,'soma').*strcmp(T.manipulation,'control')...
    .*(~contains(T.behavior_performance,'probe')).*strcmp(T.used_as_data,'yes')...
    .*(strcmp(T.cell_type,'M2')));%not probe session
%}
% ind_session=logical(contains(T.ROI_type,'trunk').*strcmp(T.manipulation,'control').*(strcmp(T.cell_type,'M2')));%SC projecting M2 neurons
% ind_session=logical(strcmp(T.ROI_type,'soma').*strcmp(T.manipulation,'control').*(strcmp(T.cell_type,'M2')));%randomly sampled M2 neuorns
% ind_session=logical(strcmp(T.ROI_type,'soma').*strcmp(T.manipulation,'control').*(strcmp(T.cell_type,'syn')).*strcmp(T.used_as_data,'yes'));%SC neurons
% ind_session=logical(strcmp(T.ROI_type,'soma').*strcmp(T.manipulation,'control').*(strcmp(T.cell_type,'vglut2')).*strcmp(T.used_as_data,'yes'));%SC neurons
ind_session=logical(strcmp(T.ROI_type,'soma').*strcmp(T.manipulation,'control').*(strcmp(T.cell_type,'vgat')).*strcmp(T.used_as_data,'yes'));%SC neurons
Tchoose=T(ind_session,:);
% file_path='H:\2P\pyx349_20210418\im_data_reg\result_save';
% session='pyx349_20210418';

behEventAlignPool={'delay onset','go cue','first lick'};
masklickPool={'yes','no','no'};
i_selectivity=3;%3-choice,4-sensory
behEventSortPool={'go cue','first lick','go cue'};
trial2include='all';
trial2exclude=[];
activity_form='dff';%{'dff','spkr'};
chosen_result={'Correct','Error'};%{'Correct','Error','Miss','Violation'};
chosen_stim={'ipsi','contra'};%{'ipsi','contra'};
n_trial_show=50;
ind_ROI=1;

%method 1, directly give path and session name
%{
file_path='E:\2P\CD058_20180127\im_data_reg\result_save';
session='CD058_20180127';
%}

%plot
%method 2, define the index, useful for loop
%
[frac_event_all,frac_event_start_all]= deal(cell(size(Tchoose,1),1));
for i=1:size(Tchoose,1)
    file_path=Tchoose.file_path{i};
    session=[Tchoose.session{i},'_',Tchoose.field{i}];
    
    savepath=[file_path,filesep,session];
    objsession=Session2P(session,file_path,trial2include,trial2exclude);
    behEventFrameIndex=objsession.behaviorEvent.behEventFrameIndex;
    sample_method='all';
    criteria='2STD';
    signals='spkr';
    threshold_dur=3;%continues 3 frames above 2 times std
%     [fract_event_all{i},frac_event_start_all{i}]=objsession.mSignificant(signals,criteria,threshold_dur,sample_method,'show_case');
    [frac_event_all{i},frac_event_start_all{i}]=objsession.mSignificant(signals,criteria,threshold_dur,sample_method,'not_show_case');
end
frac_event_all=cell2mat(frac_event_all);
frac_event_start_all=cell2mat(frac_event_start_all);
%{
% calculate a significant matrix
[significant_matrix,start_sig_mat] = fSignificant(objsession.spkr,'2STD',3);%continues 3 frames above 2 times std
%% compoare duration/start probability of significant signals happen in different epochs
indITI=[];
indSound=[];
indDelay=[];
indLick=[];
%% method 1, randomly sample time points
%{
for indTrial=1:objsession.metadata.ntr-1
    frameNum=[behEventFrameIndex.stimOnset(indTrial)-behEventFrameIndex.start(indTrial),...
        behEventFrameIndex.stimOffset(indTrial)-behEventFrameIndex.stimOnset(indTrial),...
        behEventFrameIndex.go(indTrial)-behEventFrameIndex.stimOffset(indTrial)];
    nFrameEachEpoch=min(frameNum);
    indITI_current=randperm(behEventFrameIndex.stimOnset(indTrial)-behEventFrameIndex.start(indTrial),nFrameEachEpoch)+behEventFrameIndex.start(indTrial);
    indITI=[indITI,indITI_current];
    indSound_current=randperm(behEventFrameIndex.stimOffset(indTrial)-behEventFrameIndex.stimOnset(indTrial),nFrameEachEpoch)+behEventFrameIndex.stimOnset(indTrial);
    indSound=[indSound,indSound_current];
    indDelay_current=randperm(behEventFrameIndex.go(indTrial)-behEventFrameIndex.stimOffset(indTrial),nFrameEachEpoch)+behEventFrameIndex.stimOffset(indTrial);
    indDelay=[indDelay,indDelay_current];
    indLick_current=randperm(30,nFrameEachEpoch)+behEventFrameIndex.go(indTrial);
    indLick=[indLick,indLick_current];
end
%}
%% method 2, sample all the points and get a mean
%{
for indTrial=1:objsession.metadata.ntr-1
    indITI_current=behEventFrameIndex.start(indTrial):behEventFrameIndex.stimOnset(indTrial);
    indITI=[indITI,indITI_current];
    indSound_current=behEventFrameIndex.stimOnset(indTrial):behEventFrameIndex.stimOffset(indTrial);
    indSound=[indSound,indSound_current];
    indDelay_current=behEventFrameIndex.stimOffset(indTrial):behEventFrameIndex.go(indTrial);
    indDelay=[indDelay,indDelay_current];
    indLick_current=behEventFrameIndex.go(indTrial):30+behEventFrameIndex.go(indTrial);
    indLick=[indLick,indLick_current];
end
%}

%% calculate the probability of significant event or their start
ind_all=1:size(significant_matrix,2);
chosenFrame_ITI=ismember(ind_all,indITI);
chosenFrame_sound=ismember(ind_all,indSound);
chosenFrame_delay=ismember(ind_all,indDelay);
chosenFrame_lick=ismember(ind_all,indLick);

frac_event=zeros(size(significant_matrix,1),4);%2d, ITI, sound, delay, lick
frac_event_start=zeros(size(start_sig_mat,1),4);%2d, ITI, sound, delay, lick
for iROI=1:size(significant_matrix,1)
    frac_event(iROI,1)=sum(significant_matrix(iROI,chosenFrame_ITI))/sum(chosenFrame_ITI);
    frac_event(iROI,2)=sum(significant_matrix(iROI,chosenFrame_sound))/sum(chosenFrame_sound);
    frac_event(iROI,3)=sum(significant_matrix(iROI,chosenFrame_delay))/sum(chosenFrame_delay);
    frac_event(iROI,4)=sum(significant_matrix(iROI,chosenFrame_lick))/sum(chosenFrame_lick);
    frac_event_start(iROI,1)=sum(start_sig_mat(iROI,chosenFrame_ITI))/sum(chosenFrame_ITI);
    frac_event_start(iROI,2)=sum(start_sig_mat(iROI,chosenFrame_sound))/sum(chosenFrame_sound);
    frac_event_start(iROI,3)=sum(start_sig_mat(iROI,chosenFrame_delay))/sum(chosenFrame_delay);
    frac_event_start(iROI,4)=sum(start_sig_mat(iROI,chosenFrame_lick))/sum(chosenFrame_lick);
end
%}
%% summarize across sessions and compare the proportion of significant activities happened
fig_SigProb=figure;
set(gcf,'Position',[100,100,600,300]);
ax1=subplot(1,2,1);
ax1=fScatterStat(ax1, frac_event_all,'Probability of significant activity');
ax2=subplot(1,2,2);
ax2=fScatterStat(ax2, frac_event_start_all,'Probability of significant activity start');
%}

%% help function
function ax=fScatterStat(ax, frac_event,ylabelstr)
axes(ax)
hold on;
color_epoch={[0,0,0],[1,0,0],[0,1,0],[0,0,1]};
for i=1:4
    scatter(ones(size(frac_event,1),1)*i,frac_event(:,i),10,color_epoch{i});
end
set(gca,'Xlim',[0,5],'XTick',1:4,'XTickLabel',{'ITI','sound','delay','lick'});
ylabel(ylabelstr);
y_lim=get(gca,'Ylim');
p2=signrank(frac_event(:,1),frac_event(:,2));
p3=signrank(frac_event(:,1),frac_event(:,3));
p4=signrank(frac_event(:,1),frac_event(:,4));
text(1.5,y_lim(end)*0.85,plabelsymbol(p2));
plot([1,2],[y_lim(end)*0.83,y_lim(end)*0.83],'k');
text(2,y_lim(end)*0.9,plabelsymbol(p3));
plot([1,3],[y_lim(end)*0.88,y_lim(end)*0.88],'k');
text(2.5,y_lim(end)*0.95,plabelsymbol(p4));
plot([1,4],[y_lim(end)*0.93,y_lim(end)*0.93],'k');
set(gca,'FontSize',10);
end

