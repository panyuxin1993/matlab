%plot PCA of each session, dividing different trial types
%plot dynamic memory assembly, how neurons dynamicly fire during different
%trials.
close all;
%align to stim and delay
% behEventAlign='delay onset';
% frameNumTime=[1,1.5];
% masklick='yes';
% behEventSort='go cue';
%align to go cue
behEventAlign='go cue';
frameNumTime=[1,2];
masklick='no';
behEventSort='first lick';

trial2include='all';
trial2exclude=[];
activity_form='spkr';%{'dff','spkr'};
activity_smooth_binsize=400;%ms, 200for dff, 400 for spkr
chosen_result={'Correct','Error'};%{'Correct','Error','Miss','Violation'};
chosen_stim={'ipsi','contra'};%{'ipsi','contra'};
n_trial_show=50;
ind_ROI=1;
celltype={'M2','vglut2','vgat'};
binsize=7;%used for PCA population analysis
binstep=7;

for i_celltype=1:length(celltype)
    [num,txt,raw] =xlsread('C:\Users\PYX\Documents\DataSummary\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
    T=cell2table(raw(2:end,1:15));
    T.Properties.VariableNames=strrep(raw(1,1:15),' ','_');%table variable name can't have ' ',so replace them
    ind_session=logical(strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control')...
        .*(strcmp(T.cell_type,celltype{i_celltype})+strcmp(T.cell_type,[celltype{i_celltype},'-flpo']))...
        .*contains(T.ROI_type,'soma').*(~contains(T.field,'soma')));%not probe session
    Tchoose=T(ind_session,:);
    
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
    PCA_explained_varience=cell(1,size(Tchoose,1));
    CD_explained_varience=cell(1,size(Tchoose,1));
    session_name=cell(1,size(Tchoose,1));
    dot_CD_PCA_cell=cell(1,size(Tchoose,1));
    for isession=1:size(Tchoose,1)
        file_path=Tchoose.file_path{isession};
        session=[Tchoose.session{isession},'_',Tchoose.field{isession}];
        savepath=[file_path,filesep,session];
        trial2include='all';
        trial2exclude=[];
        objsession=Session2P(session,file_path,trial2include,trial2exclude);
        neuralActivityRaw=objsession.mSmoothFR(activity_form,activity_smooth_binsize);%{objsession.zscored_spkr,objsession.dff,objsession.spkr}%choose among several acitivities form
        switch activity_form
            case 'dff'
                str_colorBar='\it\DeltaF/F';
            case 'spkr'
                str_colorBar='Firing rate';
            case 'zscored_spkr'
                str_colorBar='Z-scored firing rate';
        end
        
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
        
        n_ROI=size(neuralActivityRaw,1);
        activities_aligned_cell=cell(n_ROI,1);%each cell represent a ROI
        for iROI=1:n_ROI
            if strcmp(behEventAlign,'delay onset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
                [ activities_aligned_cell{iROI}, behEvent_aligned ] = fAlignDelaySigal( neuralActivityRaw(iROI,:), behEventFrameIndex,  frameNum );
            else
                [activities_aligned_cell{iROI}, behEvent_aligned] = fAlignSigBehEvent( neuralActivityRaw(iROI,:), behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
            end
        end
        [trialType,rule] = fGetTrialType( objsession.Data_extract,[],3,'matrix','left','divideCorErr');%decide trial type, 1d cor/err/miss/vio, each 2d ipsi/contra, 3d trials
        figSigActivity(isession)=figure();
        set(gcf,'Position',[100,100,1200,800]);
        titlestr={'ipsi','contra'};
        population_activities_cell=cell(2,1);
        population_activities_1trial_cell=cell(2,1);
        population_mean_trajectory=cell(2,1);
        y_lim1=0;
        y_lim2=0;
        resultstr={'correct','error'};
        dataPCA=[];%merge the activities from different trial types, nROI-by-nTotalFrame
        ind_trialtype=[];%n-by-m 1/0 matrix, n- total trial number, m- trial type number
        ind_singletrial=[];%a-by-n 1/0 matrix, a- total frame number, n- total trial number
        
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
                population_mean_trajectory{nStim}=cellfun(@nanmean, activities_selected_cell,'UniformOutput',false);
                population_mean_trajectory{nStim}=cell2mat(population_mean_trajectory{nStim});
                activities_sorted_cell=cellfun(@(x) x(I,:), activities_selected_cell,'UniformOutput',false);
                nTrial=sum(selectedTrialInd);
                trial_length=size(activities_sorted_cell{1},2);
                nROI=size(activities_sorted_cell,1);
                population_activities=[];
                population_activities_1trial_cell{nStim}=cell(1,nTrial);
                indTrial_currentTotal=[];
                for i=1:nTrial
                    population_activities_1trial=cellfun(@(x) x(i,:),activities_sorted_cell,'UniformOutput',false);
                    population_activities_1trial=cell2mat(population_activities_1trial);
                    population_activities_1trial_cell{nStim}{i}=population_activities_1trial;
                    %some binning to reduce noise and perform population analysis
                    ts_raw=(-frameNum(1):frameNum(2))*frT/1000;
                    [population_activities_1trial_binned,ts_binned] = fBinActivity(population_activities_1trial,binsize,binstep,ts_raw);
                    population_activities=[population_activities,population_activities_1trial_binned];
                    indTrial_current1trial=ones(size(population_activities_1trial_binned,2),1);
                    indTrial_currentTotal(size(indTrial_currentTotal,1)+1:size(indTrial_currentTotal,1)+size(population_activities_1trial_binned,2),size(indTrial_currentTotal,2)+1)=indTrial_current1trial;
                    go(i)=go(i)+trial_length*(i-1);
                    sot(i)=sot(i)+trial_length*(i-1);
                    flt(i)=flt(i)+trial_length*(i-1);
                    flt_l(i)=flt_l(i)+trial_length*(i-1);
                    flt_r(i)=flt_r(i)+trial_length*(i-1);
                    rwt(i)=rwt(i)+trial_length*(i-1);
                end
                dataPCA=[dataPCA,population_activities];
                tmp_sz=size(ind_trialtype);
                ind_trialtype(tmp_sz(1)+1:tmp_sz(1)+nTrial,tmp_sz(2)+1)=ones(nTrial,1);
                tmp_sz1=size(ind_singletrial);
                ind_singletrial(tmp_sz1(1)+1:tmp_sz1(1)+size(population_activities,2),tmp_sz1(2)+1:tmp_sz1(2)+nTrial)=indTrial_currentTotal;
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
        [nTotalTrial,nTrialType]=size(ind_trialtype);
        [coeff_PCA,score,latent,tsquared,explained,mu] = pca(dataPCA');
        PCA_explained_varience{isession}=explained;
        
        %plot conponent change and trajectory for top 3 PC
        figPCA(isession)=figure();
        set(gcf,'Position',[100,100,600,300]);
        color_trialType={'b','r'};%ipsi, contra
        lineStyle_trialType={'-','-'};%{'-','--'};%correct, error
        for i_trial=1:size(ind_singletrial,2)
            indFrame_currentTrial=logical(ind_singletrial(:,i_trial));
            ind_trialType=logical(ind_trialtype(i_trial,:));
            for iPC=1:3%left panels, plot each PCs change
                subplot(3,2,2*iPC-1);
                hold on;
                plot(score(indFrame_currentTrial,iPC),'LineStyle',lineStyle_trialType{ind_trialType},'Color',color_trialType{ind_trialType});
                title(['PC',num2str(iPC)]);
            end
            subplot(3,2,[2,4,6]);%plot 3D trajectory
            plot3(score(indFrame_currentTrial,1),score(indFrame_currentTrial,2),score(indFrame_currentTrial,3),...
                'LineStyle',lineStyle_trialType{ind_trialType},'Color',color_trialType{ind_trialType});
            hold on;
        end
        %saveas(figPCA(isession),['E:\2P\summary\population_activities\cases\',objsession.name,'_',activity_form,'_',resultstr{nResult},'_align_',behEventAlign,'population_trajectory.png'],'png');
        session_name{isession}=objsession.name;
        
        %calculate coding direction
        ts=-frameNumTime(1):frT/1000:frameNumTime(2);
        n_timepoint=min(size(population_mean_trajectory{1},2),size(population_mean_trajectory{1},2));
        coding_direction=population_mean_trajectory{1}(:,1:n_timepoint)-population_mean_trajectory{2}(:,1:n_timepoint);%when aligned to delay 
        cd_corrcoef=corrcoef(coding_direction);%each column one variable,row neurons
        figCD(isession)=figure();
        set(gcf,'Position',[100,100,600,400]);
        subplot(2,2,1);%correlation of coding direction, to see if cd keep stable
        imagesc(ts,ts,flipud(cd_corrcoef));
        hold on;
        chandle=colorbar;
        chandle.Label.String= 'CD corrcoef';
        plot(frameNum(1)*[1,1],[0,sum(frameNum)],'k--');
        plot([0,sum(frameNum)],frameNum(1)*[1,1],'k--');
        if strcmp(behEventAlign,'delay onset')
            plot((frameNum(1)-0.5*1000/frT)*[1,1],[0,sum(frameNum)],'k--');
            plot([0,sum(frameNum)],(frameNum(1)-0.5*1000/frT)*[1,1],'k--');
            frameIndStimOn=ceil(frameNum(1)-0.5*1000/frT+1);
            frameIndDelayOff=frameIndStimOn+2*1000/frT;
            frameIndLick=min(frameIndDelayOff+1000/frT, sum(frameNum));
            xlabelstr='Time (s) from delay onset';
        elseif strcmp(behEventAlign,'go cue')
            frameIndStimOn=frameNum(1)+1;%in this case, choose mean of cd during licking as coding direction
            frameIndDelayOff=frameIndStimOn+1*1000/frT;
            frameIndLick=min(frameIndDelayOff+1000/frT, sum(frameNum));
            xlabelstr='Time (s) from go cue';
        end
        ts_show=get(gca,'XTick');
        ind_show=zeros(size(ts_show));
        for i=1:length(ts_show)
            ind_show(i)=find(abs(ts_show(i)-ts)==min(abs(ts_show(i)-ts)));
        end
        ind_show_y=length(ts)+1-ind_show;
        set(gca,'YTick',ts(fliplr(ind_show_y)),'YTickLabel',fliplr(ts_show));
%         set(gca,'XTick',frameNum(1)-floor(frameNum(1)):1000/frT:sum(frameNum),'XTickLabel',-floor(frameNumTime(1)):sum(frameNumTime)-floor(frameNumTime(1)));
%         set(gca,'YTick',frameNum(1)-floor(frameNum(1)):1000/frT:sum(frameNum),'YTickLabel',-floor(frameNumTime(1)):sum(frameNumTime)-floor(frameNumTime(1)));
        ind0=find(min(abs(ts))==abs(ts));
        if abs(ts(ind0))<0.1%which mean is not epoch SVM, zero point exist
            ind0y=length(ts)+1-ind0(1);
            plot([0,0],[ts(1),ts(end)],'k--');
            plot([ts(1),ts(end)],[ts(ind0y),ts(ind0y)],'k--');
        end
        
        subplot(2,2,2);%see different trial's activities projection on coding direction, show individual cases
        hold on;
        %%%%%%%%%%%%-----------------choose one method------------%%%%%%%%%%%%
        cd_unit=nanmean(coding_direction(:,frameIndStimOn:frameIndDelayOff),2);%mean from stim to delay period
        cd_str='delay';
%         cd_unit=nanmean(coding_direction(:,frameIndDelayOff:frameIndLick),2);%mean from go to lick
%         cd_str='lick';
        cd_unit=cd_unit/sqrt(dot(cd_unit,cd_unit));
        cat_activity_raw=[];
        cat_activity_cd=[];
        for nStim=1:2 %ipsi, contra
            nTrial=length(population_activities_1trial_cell{nStim});
            if exist('cd_unit_cell','var')
                clear cd_unit_cell;
            end
            [cd_unit_cell{1:nTrial}]=deal(cd_unit);
            population_activities_cd=cellfun(@fCDproject, population_activities_1trial_cell{nStim},cd_unit_cell,'UniformOutput',false);
            for iTrial=1:nTrial
                plot(ts,population_activities_cd{iTrial},'Color',color_trialType{nStim});
                cat_activity_raw=[cat_activity_raw,population_activities_1trial_cell{nStim}{iTrial}];
                cat_activity_cd=[cat_activity_cd,population_activities_cd{iTrial}];
            end
        end
        %see projection of cd on PCs, length, no sign
        dot_cd_pc=zeros(size(coeff_PCA,2),1);
        for i_pc=1:size(coeff_PCA,2)
            dot_cd_pc(i_pc)=abs(dot(coeff_PCA(:,i_pc),cd_unit));
        end
        xlabel('Time (s)');
        ylabel('CD projection');
        %calculate varience explained by the CD projection
        CD_explained_varience{isession}=var(cat_activity_cd,0,2,'omitnan')/sum(var(cat_activity_raw,0,2,'omitnan'));
        
        subplot(2,2,3);%see different coding direction projection on different PCs
        plot(dot_cd_pc);
        xlabel('PC ID');
        ylabel('Projection of coding direction vector on PC');
        dot_CD_PCA_cell{isession}=dot_cd_pc;
        
        subplot(2,2,4);%see different trial's activities projection on coding direction, show mean+se
        for nStim=1:2 %ipsi, contra
            nTrial=length(population_activities_1trial_cell{nStim});
            if exist('cd_unit_cell','var')
                clear cd_unit_cell;
            end
            [cd_unit_cell{1:nTrial}]=deal(cd_unit);
            population_activities_cd=cellfun(@fCDproject, population_activities_1trial_cell{nStim},cd_unit_cell,'UniformOutput',false);
            population_activities_cd=reshape(population_activities_cd,[],1);
            population_activities_cd=cell2mat(population_activities_cd);
            fPlotMean_SE( ts,population_activities_cd,color_trialType{nStim} );
        end
        xlabel(xlabelstr);
        ylabel('CD projection');
        
        figure(figSigActivity(isession));
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
        ColBar.Label.String=str_colorBar;
        saveas(figSigActivity(isession),['E:\2P\summary\population_activities\cases\',objsession.name,'_',activity_form,'_',resultstr{nResult},'_align_',behEventAlign,'.pdf'],'pdf');
        
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
    %save figure to PPT using exportToPPT
    %%
    pptFileName=['E:\2P\summary\population_activities\population_activities_cases_',activity_form,'_align-',behEventAlign,'_project-',cd_str,'.pptx'];
    isOpen  = exportToPPTX();
    if ~isempty(isOpen)
        % If PowerPoint already started, then close first and then open a new one
        exportToPPTX('close');
    end
    slidesize=[16 9];
    if exist(pptFileName,'file')
        exportToPPTX('open',pptFileName);
    else
        exportToPPTX('new','Dimensions',slidesize, ...
            'Title','Population_activities_cases', ...
            'Author','PYX', ...
            'Comments','This file has been automatically generated by exportToPPTX');
    end
    
    
    for isession=1:size(Tchoose,1)
        slideId = exportToPPTX('addslide');
        exportToPPTX('addpicture',figPCA(isession),'Position',[1 1 6 3]);
        exportToPPTX('addpicture',figCD(isession),'Position',[8 1 6 4]);
        exportToPPTX('addtext',['CD explaining varience=',num2str(CD_explained_varience{isession})],'Position',[8 0 3 0.5]);
        exportToPPTX('addtext',session_name{isession},'Position',[1 0 3 0.5],'Vert','top');
        exportToPPTX('addtext',char(celltype{i_celltype}),'Position',[slidesize(1)-1 0 3 0.5],'Vert','top');
    end
    
    
    % Save and close (in one command)
    newFile = exportToPPTX('saveandclose',pptFileName);
    fprintf('New file has been saved: <a href="matlab:open(''%s'')">%s</a>\n',newFile,newFile);
    close all;
    
    %}
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
    %plot PCA explained varience percentage
    %method1, PC1-10
    PCAexplained_pc1_10=cellfun(@(x) x(1:10), PCA_explained_varience,'UniformOutput',false);
    PCAexplained_pc1_10=cell2mat(PCAexplained_pc1_10);
    dot_CD_PCA_pc1_10=cellfun(@(x) x(1:10),dot_CD_PCA_cell,'UniformOutput',false);
    dot_CD_PCA_pc1_10=cell2mat(dot_CD_PCA_pc1_10);
    %method2, PC rank percentage, 0-100%
    PCAexplained_pc0_perct100=cellfun(@(x) interp1(1:length(x),x,linspace(1,length(x))'), PCA_explained_varience,'UniformOutput',false);
    PCAexplained_pc0_perct100=cellfun(@(x) x/sum(x)*100,PCAexplained_pc0_perct100,'UniformOutput',false);
    PCAexplained_pc0_perct100=cell2mat(PCAexplained_pc0_perct100);
    dot_CD_PCA_pc0_perct100=cellfun(@(x) interp1(1:length(x),x,linspace(1,length(x))'), dot_CD_PCA_cell,'UniformOutput',false);
    dot_CD_PCA_pc0_perct100=cellfun(@(x) x/sum(x)*100,dot_CD_PCA_pc0_perct100,'UniformOutput',false);
    dot_CD_PCA_pc0_perct100=cell2mat(dot_CD_PCA_pc0_perct100);
    
    save(['E:\2P\summary\population_activities\PCA_explained_varience_',celltype{i_celltype},'_',resultstr{nResult},'_using_',activity_form,'.mat'],...
        'PCA_explained_varience','CD_explained_varience','PCAexplained_pc1_10','PCAexplained_pc0_perct100','dot_CD_PCA_cell','dot_CD_PCA_pc1_10','dot_CD_PCA_pc0_perct100');
end

%% compare PCA explained varience for different cell types
figPCAexplained=figure();
set(gcf,'Position',[100,100,500,500]);
color_celltype={'k','b','r'};%M2,SC vglut2, SC vgat
celltype_str={'M2','SC E','SC I'};
for i_celltype=1:length(celltype)
    load(['E:\2P\summary\population_activities\PCA_explained_varience_',celltype{i_celltype},'_',resultstr{nResult},'_using_',activity_form,'.mat']);
    subplot(2,2,1);%plot PC1-10
    curve_PCAexplained_pc1_10(i_celltype) = fPlotMean_SE( 1:10,PCAexplained_pc1_10',color_celltype{i_celltype} );
    xlabel('PCs');
    ylabel('Percent varience explained');
    hold on;
    set(gca,'Xlim',[0,5],'FontSize',10);
    subplot(2,2,2);%plot PC rank percentage, 0-100%
    curve_PCAexplained_pc0_perct100(i_celltype) = fPlotMean_SE( 1:100,PCAexplained_pc0_perct100',color_celltype{i_celltype} );
    xlabel('Percentile of PCs rank');
    set(gca,'Xlim',[0,20],'FontSize',10);
    hold on;
    subplot(2,2,3);%plot projection of CD on PCs, PC1-10
    curve_PCAexplained_pc1_10(i_celltype) = fPlotMean_SE( 1:10,dot_CD_PCA_pc1_10',color_celltype{i_celltype} );
    xlabel('PCs');
    ylabel('CD projection');
    hold on;
    set(gca,'Xlim',[0,10],'FontSize',10);
    subplot(2,2,4);%plot projection of CD on PCs, PC rank percentage, 0-100%
    curve_PCAexplained_pc0_perct100(i_celltype) = fPlotMean_SE( 1:100,dot_CD_PCA_pc0_perct100',color_celltype{i_celltype} );
    xlabel('Percentile of PCs rank');
    set(gca,'Xlim',[0,30],'FontSize',10);
    hold on;
end
subplot(2,2,4);
hl=legend(curve_PCAexplained_pc0_perct100,celltype_str);
set(hl,'box','off');

saveas(gcf,['E:\2P\summary\population_activities\PCA_varience_explain_using_',activity_form,'.pdf'],'pdf');


%% helper function
function projection = fCDproject(x,y)
yx=repmat(y,1,size(x,2));
projection=dot(x,yx,1);
end