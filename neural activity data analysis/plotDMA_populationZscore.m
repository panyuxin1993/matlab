%plot dynamic memory assembly, how neurons dynamicly fire during different
%trials. 
close all;
behEventAlign='delay onset';
frameNumTime=[0.5,2];

masklick='yes';
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
ind_session=logical(strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*strcmp(T.cell_type,'M2').*contains(T.ROI_type,'soma'));%not probe session
Tchoose=T(ind_session,:);
% file_path='H:\2P\pyx349_20210418\im_data_reg\result_save';
% session='pyx349_20210418';

trial2include='all';
trial2exclude=[];
activity_form='spkr';%{'dff','spkr'};
chosen_result={'Correct','Error'};%{'Correct','Error','Miss','Violation'};
chosen_stim={'ipsi','contra'};%{'ipsi','contra'};
n_trial_show=50;
ind_ROI=1;

%plot
%method 1, define the index
%{
i=1;
file_path=Tchoose.file_path{i};
session=[Tchoose.session{i},'_',Tchoose.field{i}];
%}

%method 2, directly give path and session name
%{
file_path='E:\2P\CD058_20180127\im_data_reg\result_save';
session='CD058_20180127';
%}

neuralActivityMean=cell(2,size(Tchoose,1));
for isession=1:size(Tchoose,1)
    file_path=Tchoose.file_path{isession};
    session=[Tchoose.session{isession},'_',Tchoose.field{isession}];
    savepath=[file_path,filesep,session];
    objsession=Session2P(session,file_path,trial2include,trial2exclude);
%     zscored_spkr=objsession.zscored_spkr;
    zscored_spkr=objsession.dff;
    frT = objsession.metadata.frT;
    frameNum=double(round(frameNumTime*1000/frT));
    [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( objsession.Data_extract, objsession.metadata.ind_1stFrame, frT ,objsession.metadata.ind_tr_1);%get behavior event time
    
    if isnan(trial2include)
        indTrial2include=fExcludeTrials(trial2exclude,objsession.metadata.ind_1stFrame,'logical');
        trial2include=[1,length(objsession.metadata.ind_1stFrame)];
    elseif strcmp(trial2include,'all')
        trial2include=[1,length(objsession.metadata.ind_1stFrame)];
        indTrial2include=fIncludeTrials(trial2include,objsession.metadata.ind_1stFrame,'logical');
    else
        indTrial2include=fIncludeTrials(trial2include,objsession.metadata.ind_1stFrame,'logical');
    end
    
    n_ROI=size(zscored_spkr,1);
    activities_aligned_cell=cell(n_ROI,1);%each cell represent a ROI
    for iROI=1:n_ROI
        %if strcmp(behEventAlign,'delay onset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
        [ activities_aligned_cell{iROI}, behEvent_aligned ] = fAlignDelaySigal( zscored_spkr(iROI,:), behEventFrameIndex,  frameNum );
        %[activities_aligned_cell{iROI}, behEvent_aligned] = fAlignSigBehEvent( zscored_spkr(iROI,:), behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
    end
    [trialType,rule] = fGetTrialType( objsession.Data_extract,[],3,'matrix','left','divideCorErr');%decide trial type, 1d cor/err/miss/vio, each 2d ipsi/contra, 3d trials
    figure;
    set(gcf,'Position',[100,100,1200,800]);
    titlestr={'ipsi','contra'};
    population_activities_cell=cell(2,1);
    population_activities_1trial_cell=cell(2,1);
    y_lim1=0;
    y_lim2=0;
    resultstr={'correct','error'};
    for nStim=1:2 %ipsi, contra
        for  nResult=1 %4 column(correct/error/miss/violation),only show correct trials
            selectedTrialInd=trialType(nResult,nStim,:);
            selectedTrialInd=reshape(logical(squeeze(selectedTrialInd)),[],1);
            indTrial2include=reshape(indTrial2include,[],1);
            selectedTrialInd=logical(selectedTrialInd.*indTrial2include);
            sot=behEvent_aligned.stimOnset(selectedTrialInd);% stim onset, white
            flt=behEvent_aligned.lickFirst(selectedTrialInd);% first lick, black
            flt_l=behEvent_aligned.lickFirst_left(selectedTrialInd);%left lick, magenta dots
            flt_r=behEvent_aligned.lickFirst_right(selectedTrialInd);%right lick, red dots
            rwt=behEvent_aligned.rewTime(selectedTrialInd);%reward, cyan dots
            go=behEvent_aligned.go(selectedTrialInd);%go cue,white
            temp=find(selectedTrialInd);
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
            activities_selected_cell=cellfun(@(x) x(selectedTrialInd,:),activities_aligned_cell,'UniformOutput',false);
            activities_sorted_cell=cellfun(@(x) x(I,:), activities_selected_cell,'UniformOutput',false);
            nTrial=sum(selectedTrialInd);
            trial_length=size(activities_sorted_cell{1},2);
            nROI=size(activities_sorted_cell,1);
            population_activities=[];
            population_activities_1trial_cell{nStim}=cell(1,nTrial);
            for i=1:nTrial
                population_activities_1trial=cellfun(@(x) x(i,:),activities_sorted_cell,'UniformOutput',false);
                population_activities_1trial=cell2mat(population_activities_1trial);
                population_activities_1trial_cell{nStim}{i}=population_activities_1trial;
                population_activities=[population_activities,population_activities_1trial];
                go(i)=go(i)+trial_length*(i-1);
                sot(i)=sot(i)+trial_length*(i-1);
                flt(i)=flt(i)+trial_length*(i-1);
                flt_l(i)=flt_l(i)+trial_length*(i-1);
                flt_r(i)=flt_r(i)+trial_length*(i-1);
                rwt(i)=rwt(i)+trial_length*(i-1);
            end
            population_activities_cell{nStim}=population_activities;
            subplot(4,1,nStim*2-1);
            hold on;
            imagesc(population_activities);
            
            for i=1:nTrial
                text(trial_length*(i-0.5),nROI*1.08,num2str(i));
                plot([sot(i),sot(i)],[1,nROI],'w-','LineWidth',1);
                hold on;
                plot([go(i),go(i)],[1,nROI],'w-','LineWidth',1);
                %             plot([flt(i),flt(i)],[1,nROI],'k-','LineWidth',1);
                %             plot([flt_l(i),flt_l(i)],[1,nROI],'m-','LineWidth',1);
                %             plot([flt_r(i),flt_r(i)],[1,nROI],'r-','LineWidth',1);
                %             plot([rwt(i),rwt(i)],[1,nROI],'c-','LineWidth',1);
            end
            
            set(gca,'ytick',nROI,'yticklabel',num2str(nROI),'xtick',0:20/frT*1000:size(population_activities,2),'xticklabel',0:20:size(population_activities,2)*frT/1000,'ydir','normal');
            ts=(1:size(population_activities,2))*frT/1000;
            set(gca,'Xlim',[1,size(population_activities,2)]);
            title([titlestr{nStim},' trial']);
            xlabel('Time (s)');
            ylabel('#ROI');
                        
            subplot(4,1,nStim*2);%mean trace
            fPlotSum(ts,population_activities,'k');hold on;
            
            set(gca,'Xlim',[ts(1),ts(end)]);
            xlabel('Time (s)');
            y_lim=get(gca,'Ylim');
            y_lim1=min(y_lim1,y_lim(1));
            y_lim2=max(y_lim2,y_lim(2));
        end
    end
    temp_population_activities=[population_activities_cell{1},population_activities_cell{2}];
    ColLimit = prctile(reshape(temp_population_activities,1,[]),98);
    ColLimit_min = prctile(reshape(temp_population_activities,1,[]),5);
    for i_subplot=[1,3]
        subplot(4,1,i_subplot);
        set(gca,'clim',[0 ColLimit]);
    end
    for i_subplot=[2,4]
        subplot(4,1,i_subplot);
        set(gca,'ylim',[y_lim1,y_lim2]);
        plot([y_lim1,y_lim2],[0,0],'k-');
    end
    ColBar = colorbar;
    set(ColBar,'position',[0.91 0.76 0.02 0.14]);
    ColBar.Label.String='z-scored firing rate';
    saveas(gcf,['E:\2P\summary\population_activities\',objsession.name,resultstr{nResult},'.pdf'],'pdf');
    
    [ neuralActivityMean_cell{1}, neuralActivitySE ] = cellfun(@fMean_SE, population_activities_1trial_cell{1},'UniformOutput',false);
    [ neuralActivityMean_cell{2}, neuralActivitySE ] = cellfun(@fMean_SE, population_activities_1trial_cell{2},'UniformOutput',false);
    
    neuralActivityMean{1,isession}=cellfun(@(x) nanmean(x),neuralActivityMean_cell{1});
    neuralActivityMean{2,isession}=cellfun(@(x) nanmean(x),neuralActivityMean_cell{2});
    
    [h,p]=ttest2(neuralActivityMean{1,isession},neuralActivityMean{2,isession});
    figure;
    scatter(1*ones(length(neuralActivityMean{1,isession}),1),neuralActivityMean{1,isession},10,'k');hold on;
    scatter(2*ones(length(neuralActivityMean{2,isession}),1),neuralActivityMean{2,isession},10,'k');
    set(gca,'Xlim',[0,3],'XTick',[1,2],'XTickLabel',{'ipsi','contra'});
    text(0.5,0.9,plabelsymbol(p),'Unit','Normalized');
end
%% group neurons together
neuralActivityMean_all{1}=cell2mat(neuralActivityMean(1,:));
neuralActivityMean_all{2}=cell2mat(neuralActivityMean(2,:));
[h,p]=ttest2(neuralActivityMean_all{1},neuralActivityMean_all{2});
figure;
x=[neuralActivityMean_all{1}';neuralActivityMean_all{2}'];
g=[repmat({'ipsi'},length(neuralActivityMean_all{1}),1);repmat({'contra'},length(neuralActivityMean_all{2}),1)];
boxplot(x,g,'Notch','on');hold on;
set(gca,'Xlim',[0,3]);
text(0.5,0.9,plabelsymbol(p),'Unit','Normalized');





