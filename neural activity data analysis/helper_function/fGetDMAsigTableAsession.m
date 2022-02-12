function [tbSigTrial4Cell_out,figDMAsig,figSigCellNum,savingNameStr] = fGetDMAsigTableAsession(file_path,session,animal,date,field,celltype,trialTypeStr,activity_form,prct_delay2calculate_pResp,threshold_dur,fig_n_event,ax_event_time)
%fGetDMAsigTableAsession summary a session data to a table
%Input-
%Output-
%   Detailed explanation goes here
trial2include='all';
trial2exclude=nan;
datestr=strrep(date,'/','-');
%naming the temp file by method of AUC calculation correction
fileNameT=[file_path,filesep,animal,'-',datestr,'trialType',trialTypeStr,'TsigTrialNum4cell_prctDelay2CalPResp-',prct_delay2calculate_pResp,'.mat'];
disp(['fGetDMAsigTableAsession is processing ',session]);
%     file_path=Tchoose.file_path{isession};
%     session=[Tchoose.session{isession},'_',Tchoose.field{isession}];
%     savepath=[file_path,filesep,session];
objsession=Session2P(session,file_path,trial2include,trial2exclude);
if strcmp(activity_form,'spkr')
    activity_data=objsession.spkr;
    significant_criteria='3STD';
    sig_criteria_off='3STD';
    threshold_dur_range=100:50:400;
elseif strcmp(activity_form,'dff')
    activity_data=objsession.dff;
    significant_criteria='2STD';
    sig_criteria_off='0.5STD';
    threshold_dur_range=100:100:800;
end
frT = objsession.metadata.frT;
frameNumTime=[1,1.5];%from 0.5 before stimOnset to end of delay
frameNum=double(round(frameNumTime*1000/frT));

[behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( objsession.Data_extract, objsession.metadata.ind_1stFrame, frT ,objsession.metadata.ind_tr_1);%get behavior event time


%find index of baseline
% indBaseline = zeros(1,size(activity_data,2));
ind_stim_delay=zeros(1,size(activity_data,2));
for iTrial=1:length(behEventFrameIndex.start)
    %indBaseline(behEventFrameIndex.start(iTrial):behEventFrameIndex.stimOnset(iTrial))=1;
    ind_stim_delay(behEventFrameIndex.stimOnset(iTrial):behEventFrameIndex.delayOffset(iTrial))=1;
end
% indBaseline=logical(indBaseline);
indBaseline = objsession.mIndexBaselineAroundMode(activity_form);
ind_stim_delay=logical(ind_stim_delay);
%     %examine the creteria of folds of STD
%     inROI=1:3:objsession.metadata.nROI;%1:103;
%     indTrial=[];%in 11:20
%     figExample=fPlotF_ROI(file_path,'spkr','2P',inROI,indTrial,'sigThreshSTD',2,'sigShowingStyle','patch','BaselineIndex',indBaseline);%example for 2P contralateral preference during delay
%test which criteria for duration is suitable
n_sigevent=[];
n_sigevent_SD=[];%summary only the significant event during sim. and delay period
binSize=1;
binStep=1;
ts=(-frameNum(1):binStep:frameNum(2))*frT/1000;
frameNumTime_go=[0.5,1];%from 0.5 before go cue to 1s after go cue
frameNum_go=double(round(frameNumTime_go*1000/frT));
ts_go=(-frameNum_go(1):binStep:frameNum_go(2))*frT/1000;
%see how different threshold of continuing significant activities affect
%the significant event characteristics
%
for threshold_dur_temp=threshold_dur_range
%     [significant_matrix,start_sig_mat] = fSignificant(activity_data,'2STD',round(threshold_dur_temp/frT),'BaselineIndex',indBaseline);
    [significant_matrix,start_sig_mat] = fSignificant(activity_data,significant_criteria,round(threshold_dur_temp/frT),'BaselineIndex',indBaseline,'EventOffCriteria',sig_criteria_off);
    n_sigevent=[n_sigevent,sum(start_sig_mat,'all')/size(start_sig_mat,1)];
    n_sigevent_SD=[n_sigevent_SD,sum(start_sig_mat(:,ind_stim_delay),'all')/size(start_sig_mat,1)];
    start_sig_rate=fEventRate(start_sig_mat,binSize,binStep);
    sig_rate=fEventRate(significant_matrix,binSize,binStep);
    start_sig_rate_aligned = fAlignDelaySigal( start_sig_rate, behEventFrameIndex, frameNum,'raw' );
    axes(ax_event_time(1));
    fPlotMean_SE(ts,start_sig_rate_aligned,[1,1,1]*(threshold_dur_range(end)-threshold_dur_temp)/threshold_dur_range(end));
    hold on;
    start_sig_rate_aligned_go = fAlignSigBehEvent( start_sig_rate, behEventFrameIndex, lickingFrameIndex, 'go cue',frameNum_go,'raw' );
    axes(ax_event_time(2));
    fPlotMean_SE(ts_go,start_sig_rate_aligned_go,[1,1,1]*(threshold_dur_range(end)-threshold_dur_temp)/threshold_dur_range(end));
    hold on;
    sig_rate_aligned = fAlignDelaySigal( sig_rate, behEventFrameIndex, frameNum,'raw' );
    axes(ax_event_time(3));
    fPlotMean_SE(ts,sig_rate_aligned,[1,1,1]*(threshold_dur_range(end)-threshold_dur_temp)/threshold_dur_range(end));
    hold on;
    sig_rate_aligned_go = fAlignSigBehEvent( sig_rate, behEventFrameIndex, lickingFrameIndex, 'go cue',frameNum_go,'raw' );
    axes(ax_event_time(4));
    fPlotMean_SE(ts_go,sig_rate_aligned_go,[1,1,1]*(threshold_dur_range(end)-threshold_dur_temp)/threshold_dur_range(end));
    hold on;
end
for i_ax=1:length(ax_event_time)
    axes(ax_event_time(i_ax));
    y_lim=get(gca,'YLim');
    plot([0,0],[0,y_lim(2)],'k-');
    if mod(i_ax,2)==1    
        set(ax_event_time(i_ax),'Xlim',[-frameNumTime(1),frameNumTime(2)]); 
        plot([-0.5,-0.5],[0,y_lim(2)],'k-');
        xlabel('Time (s) from stim. onset');
        set(gca,'XTick',[-1,0,1.5],'XTickLabel',{'-1','0','1.5'});
    else
        set(ax_event_time(i_ax),'Xlim',[-frameNumTime_go(1),frameNumTime_go(2)]);
        xlabel('Time (s) from go cue');
        set(gca,'XTick',[-0.5,0,1],'XTickLabel',{'-0.5','0','1'});
    end
    if i_ax==1
        ylabel('Freqency (Significant response start)');
    elseif i_ax==3
        ylabel('Freqency (Significant response)');
    end
end


figure(fig_n_event);
subplot(1,2,1);
plot(threshold_dur_range,n_sigevent);hold on;
% title('All significant events');
subplot(1,2,2);
plot(threshold_dur_range,n_sigevent_SD);hold on;
% title('Significant events during stim. and delay period');
%}
savingNameStr=['_BinSize',num2str(binSize),'_BinStep',num2str(binStep)];

%choose the final critera
[significant_matrix,start_sig_mat] = fSignificant(activity_data,significant_criteria,round(threshold_dur/frT),'BaselineIndex',indBaseline,'EventOffCriteria',sig_criteria_off);

if isnan(trial2include)
    indTrial2include=fExcludeTrials(trial2exclude,objsession.metadata.ind_1stFrame,'logical');
    trial2include=[1,length(objsession.metadata.ind_1stFrame)];
elseif strcmp(trial2include,'all')
    trial2include=[1,length(objsession.metadata.ind_1stFrame)];
    indTrial2include=fIncludeTrials(trial2include,objsession.metadata.ind_1stFrame,'logical');
else
    indTrial2include=fIncludeTrials(trial2include,objsession.metadata.ind_1stFrame,'logical');
end

n_ROI=size(activity_data,1);
n_trial=length(objsession.metadata.ind_1stFrame);
sig_trial_matrix=zeros(n_ROI,n_trial);
activities_aligned_cell=cell(n_ROI,1);
for iROI=1:n_ROI
    for iTrial=1:n_trial
        sig_start_ind=find(start_sig_mat(iROI,:));
        if sum((behEventFrameIndex.stimOnset(iTrial) < sig_start_ind).*(behEventFrameIndex.go(iTrial) > sig_start_ind))>0
            sig_trial_matrix(iROI,iTrial)=1;
        end
    end
    [ activities_aligned_cell{iROI}, behEvent_aligned ] = fAlignDelaySigal( activity_data(iROI,:), behEventFrameIndex,  frameNum );
end
[trialType,rule] = fGetTrialType( objsession.Data_extract,[],3,'matrix','left','divideCorErr');%decide trial type, 1d cor/err/miss/vio, each 2d ipsi/contra, 3d trials

%---form table and concatenate across sessions---
[population_activities_cell,delay_length]=deal(cell(2,1));%ispi, contra
%
population_activities_1trial_cell=cell(2,1);
x_lim1=0;
x_lim2=0;
resultstr={'correct','error'};
for nStim=1:2 %ipsi, contra
    switch trialTypeStr
        case 'cor'
            selectedTrialInd=trialType(1,nStim,:);
        case 'err'
            selectedTrialInd=trialType(2,nStim,:);
        case 'cor and err'
            selectedTrialInd=trialType(1:2,nStim,:);
    end
    selectedTrialInd=reshape(logical(squeeze(selectedTrialInd)),[],1);
    indTrial2include=reshape(indTrial2include,[],1);
    selectedTrialInd=logical(selectedTrialInd.*indTrial2include);
    sig_data_selected=sig_trial_matrix(:,selectedTrialInd);
    sot=behEvent_aligned.stimOnset(selectedTrialInd);%stim onset
    go=behEvent_aligned.go(selectedTrialInd);%go cue
    delay_length{nStim}=(go-sot)*frT/1000;%s
    [B,I]=sort(delay_length{nStim},'descend');%sort trial based on delay length
    delay_length{nStim}=delay_length{nStim}(I);
    population_activities_cell{nStim}=sig_data_selected(:,I);
    nTrial_selected=sum(selectedTrialInd);
    population_activities_1trial_cell{nStim}=cell(1,nTrial_selected);
    activities_selected_cell=cellfun(@(x) x(selectedTrialInd,:),activities_aligned_cell,'UniformOutput',false);
    activities_sorted_cell=cellfun(@(x) x(I,:), activities_selected_cell,'UniformOutput',false);
    population_activities_1trial_cell{nStim}=activities_sorted_cell;
end
significantTrialNum4cell_cell=cellfun(@(x) sum(x,2),population_activities_cell,'UniformOutput',false);
delay_length_rep=cellfun(@(x) repmat(x,[n_ROI,1]), delay_length,'UniformOutput',false);
varROI=(1:n_ROI);
varROI=reshape(varROI,[],1);
[varanimal{1:n_ROI}]=deal(animal);
varanimal=reshape(varanimal,[],1);
[vardate{1:n_ROI}]=deal(date);
vardate=reshape(vardate,[],1);
[varfield{1:n_ROI}]=deal(field);
varfield=reshape(varfield,[],1);
[varcelltype{1:n_ROI}]=deal(celltype);
varcelltype=reshape(varcelltype,[],1);
tbSigTrial4Cell=table(varanimal,vardate,varfield,varcelltype,varROI,...
    significantTrialNum4cell_cell{1},significantTrialNum4cell_cell{2},...
    significantTrialNum4cell_cell{2}-significantTrialNum4cell_cell{1},...
    population_activities_cell{1},population_activities_cell{2},...
    delay_length_rep{1},delay_length_rep{2},...
    population_activities_1trial_cell{1},population_activities_1trial_cell{2},...
    'VariableNames',{'animal','date','field','celltype','nROI',...
    'ipsi_sum','contra_sum','contraMipsi',...
    'ipsi_sig_data','contra_sig_data',...
    'ipsi_delay','contra_delay','ipsi_raw_data','contra_raw_data'});

prct_delay=prct_delay2calculate_pResp*ones(n_ROI,1);
tTemp=addvars(tbSigTrial4Cell,prct_delay);
tTemp=removevars(tTemp,{'animal','date','field','celltype','nROI','ipsi_sum','contra_sum','contraMipsi'});
tbSigTrialPerct=rowfun(@fSigResponse,tTemp,'OutputVariableNames',{'respSide','pSigResponse'});
tbSigTrial4Cell=[tbSigTrial4Cell,tbSigTrialPerct];
tbSigTrial4Cell=addvars(tbSigTrial4Cell,prct_delay);
%to be concatenate with other table, change variables to cell array so the
%size is consistent
tbSigTrial4Cell_out=tbSigTrial4Cell;
rowDist=ones(size(tbSigTrial4Cell_out,1),1);
tbSigTrial4Cell_out.ipsi_sig_data=mat2cell(tbSigTrial4Cell_out.ipsi_sig_data,rowDist);
tbSigTrial4Cell_out.contra_sig_data=mat2cell(tbSigTrial4Cell_out.contra_sig_data,rowDist);
tbSigTrial4Cell_out.ipsi_delay=mat2cell(tbSigTrial4Cell_out.ipsi_delay,rowDist);
tbSigTrial4Cell_out.contra_delay=mat2cell(tbSigTrial4Cell_out.contra_delay,rowDist);

%---plot cell profile showing their preference of ipsi/contra side---
[figDMAsig,figSigCellNum]=deal([]);
%{
figDMAsig=figure;
ax_raster_ipsi = axes('Position',[0.1,0.6,0.7,0.35],'Box','off');
ax_barh_ipsi = axes('Position',[0.85,0.6,0.1,0.35],'Box','off');
ax_raster_contra = axes('Position',[0.1,0.2,0.7,0.35],'Box','off');
ax_barh_contra = axes('Position',[0.85,0.2,0.1,0.35],'Box','off');
ax_stem = axes('Position',[0.1,0.05,0.7,0.1],'Box','off');
set(gcf,'Position',[100,100,600,900]);
titlestr={'Ipsi','Contra'};
axes(ax_stem);
tbSigTrial4Cell=sortrows(tbSigTrial4Cell,{'contraMipsi','contra_sum'},{'descend','descend'});
x=1:size(tbSigTrial4Cell,1);
h=stem(x,[tbSigTrial4Cell.ipsi_sum,-tbSigTrial4Cell.contra_sum],'filled','LineWidth',1);
set(gca,'Xlim',[0.5,n_ROI-0.5]);
set(h(1),'Color','blue');
set(h(2),'Color','red');
xlabel('Cell ID');
ylabel('# trials with significant responses');
for nStim=1:2 %ipsi, contra
    if nStim==1
        axes(ax_raster_ipsi);
        sig_data_selected=tbSigTrial4Cell.ipsi_sig_data;
    elseif nStim==2
        axes(ax_raster_contra);
        sig_data_selected=tbSigTrial4Cell.contra_sig_data;
    end
    hold on;
    imagesc(sig_data_selected');%transpose, so row represent trials and column represent ROIs
    nROI=size(sig_data_selected,1);
    set(gca,'xtick',nROI,'xticklabel',num2str(nROI),'xdir','normal');
    ts=1:size(sig_data_selected,2);
    set(gca,'Ylim',[1,size(sig_data_selected,2)]);
    set(gca,'Xlim',[1,size(sig_data_selected,1)]);
    ylabel(['# ',titlestr{nStim},' trials']);
    
    %barh plot of each trial
    y_data=sum(sig_data_selected);
    if nStim==1
        axes(ax_barh_ipsi);
        barh(ts,y_data,'EdgeColor','b','FaceColor','b');
    elseif nStim==2
        axes(ax_barh_contra);
        barh(ts,y_data,'EdgeColor','r','FaceColor','r');
        xlabel('# Cells');
    end
    hold on;
    set(gca,'Ylim',[1,size(sig_data_selected,2)]);
    set(gca,'Xlim',[ts(1),ts(end)]);
    %             ylabel('# Trial ID with significant response');
    x_lim=[min(y_data),max(y_data)];
    x_lim1=min(x_lim1,x_lim(1));
    x_lim2=max(x_lim2,x_lim(2));
end
temp_population_activities=[population_activities_cell{1},population_activities_cell{2}];
ColLimit = prctile(reshape(temp_population_activities,1,[]),100);
ColLimit_min = prctile(reshape(temp_population_activities,1,[]),1);
%     for i_subplot=[1,3]
%         subplot(4,1,i_subplot);
%         set(gca,'clim',[0 ColLimit]);
%     end
if x_lim2>x_lim1
    axes(ax_barh_ipsi);
    set(gca,'xlim',[x_lim1,x_lim2]);
    plot([0,0],[x_lim1,x_lim2],'k-');
    axes(ax_barh_contra);
    set(gca,'xlim',[x_lim1,x_lim2]);
    plot([0,0],[x_lim1,x_lim2],'k-');
end
ColBar = colorbar;
set(ColBar,'position',[0.91 0.05 0.02 0.1]);
ColBar.Label.String='significant response';
%figure for histogram of significant cell number
significantNum=cellfun(@(x) sum(x,1),population_activities_cell,'UniformOutput',false);
[h,p]=ttest2(significantNum{1},significantNum{2});
figSigCellNum=figure;
set(gcf,'Position',[100,200,400,200]);
subplot(1,2,1);%scatter plot comparing ipsi and contra
scatter(1*ones(length(significantNum{1}),1),significantNum{1},10,'k');hold on;
scatter(2*ones(length(significantNum{2}),1),significantNum{2},10,'k');
set(gca,'Xlim',[0,3],'XTick',[1,2],'XTickLabel',{'ipsi','contra'});
text(0.5,0.9,plabelsymbol(p),'Unit','Normalized');
subplot(1,2,2);%histogram comparing significant cell number/proportion of ipsi and contra trials
histogram(significantNum{1},'BinWidth',1,'DisplayStyle','stairs','EdgeColor',[0,0,1]);
hold on;
histogram(significantNum{2},'BinWidth',1,'DisplayStyle','stairs','EdgeColor',[1,0,0]);
%     xlabel('Trial ID with significant responses');
ylabel('# Cells');

%}
end

%------- helper functions----------------------
function [respSide, pSigResponse]= fSigResponse(ipsi_sig_data,contra_sig_data,ipsi_delay,contra_delay,ipsi_raw_data,contra_raw_data,prct_delay)
%find preferred side
mean_ipsi=nanmean(ipsi_raw_data{1},'all');%ipsi_raw_data is a 1-by-1 cell
mean_contra=nanmean(contra_raw_data{1},'all');
if mean_ipsi>mean_contra
    respSide='I';
    delay_preferred=ipsi_delay;
    data_preferred=ipsi_sig_data;
else
    respSide='C';%note here, rowfun output neet to be vercat, so should with same size if using string. 'ipsi' and 'contra' do not meet the criteria
    delay_preferred=contra_delay;
    data_preferred = contra_sig_data;
end
ind_trial=(delay_preferred>prctile(delay_preferred,prct_delay));
pSigResponse=mean(data_preferred(ind_trial));
end
%from event raster(in a 0-1 matrix) to event happen rate of vector
function eventRate=fEventRate(eventMat,binSize,binStep)
sampleNum=size(eventMat,2);
nTrial=size(eventMat,1);
eventRate=zeros(1,sampleNum);
for j=1:sampleNum
    ind1=max(1,binStep*(j-1)+1-ceil((binSize-1)/2));%in binSize=2n, then left bin longer than right bin
    ind2=min(sampleNum,binStep*(j-1)+1+floor((binSize-1)/2));
    eventRate(j+ceil(binSize/binStep))=nansum(nansum(eventMat(:,ind1:ind2)))*1000/sum(sum(~isnan(eventMat(:,ind1:ind2))));%denominator only include not nan data
end
end