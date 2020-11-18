%This code choose some sessions and plot PSTH of each cases;
%both show individual sessions and animals
% dbstop if error
clear;
close all;
%% import .csv to decide which session to include
TepochAUC=readtable('H:\FP\summary\EpochAUC.csv');
% TepochAUCsig=TepochAUC(true(size(TepochAUC.pAUClate)),:);%include all data
TepochAUCsig=TepochAUC(TepochAUC.pAUClate<0.025 | TepochAUC.pAUClate>0.975,:);%include only data with significant late AUC
nameSig=TepochAUCsig.DatapointName;
%% import data summay file
savepath='H:\FP\summary';
[num,txt,raw] =xlsread('D:\xulab\project\fiber photometry data summary.xlsx');%criteria to choose sessions come from this file
T=cell2table(raw(2:end,1:13));
T.Properties.VariableNames=strrep(raw(1,1:13),' ','_');%table variable name can't have ' ',so replace them
T.delay_min=zeros(size(T,1),1);%convert delay variable to numeric
T.delay_max=zeros(size(T,1),1);
T.retract=zeros(size(T,1),1);
for i=1:size(T,1)
    if ischar(T.delay_length{i})
        delay_temp=split(T.delay_length{i},'-');
        T.delay_min(i)=str2double(delay_temp{1});
        T.delay_max(i)=str2double(delay_temp{end});
    elseif isnumeric(T.delay_length{i})%usually only one number, means fixed delay
        T.delay_min(i)=T.delay_length{i};
        T.delay_max(i)=T.delay_length{i};
    else
        T.delay_min(i)=nan;
        T.delay_max(i)=nan;
    end
    if ischar(T.note{i})
        if contains(T.note{i},'retract')
            T.retract(i)=1;
        end
    end
end
%settings of label and legend
i_selectivity=4;%*********variable**************
selectivitystr={'stimuli','difficulty','sensory','choice'};%sensory means grouping difficulties;
trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
fiberstr={'Soma'};%fiberstr={'Soma','Terminal'};
dffstr={'dff baseline correction first','dff motion correction first','dff470','dff410'};
behEventAlign='delay onset';%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
masklick='yes';%'yes' if mask lickings or 'no' when want to see all activity
stimDur=0.5;%0.5s click duration
frameNumTime=[1,1.5];%from 2s before align point to 5s after align point
% yrange=[-0.01,0.06];
% const settings of fiber photometry
Session470LEDstartFrame=1;% usually 205
Session410LEDstartFrame=2;
FrameRate=40;
FrameTime=1000/FrameRate;
frT=FrameTime*2;%2 channel, so framerate and frametime should be half
frameNum=double(round(frameNumTime*1000/frT));
%parrameters for calculating AUC etc.
binsize=3;
binstep=1;
nshuffle=10;
p_sepOnset=0.05;%the p-value used to determine when separation happen
%set the criteria to choose sessions
i_region=3;%1-bilateral, 2-unilateral*********variable**************
region={'bilateral','t SC','SC'};%{'bilateral SC','SC'}
regionstr={'bilateral','unilateral','mixed'};
whichAnimal='pyx241';
celltype={'SC vglut2','SC vgat'};
% celltype={'ALM terminal'};
n_siteCell=cell(1,2);%each cell store a cell type
n_datapointCell=cell(1,2);%each cell store a cell type
for i_celltype=1:length(celltype)
    %select sessions for summary analysis
    experiment=celltype{i_celltype};
    if strcmp(selectivitystr{i_selectivity},'stimuli')|| strcmp(selectivitystr{i_selectivity},'difficulty')%here only include sessions with probe trials
        ind_session=(T.delay_max==1.5).*(T.delay_min==0.3).*strcmp(T.experiment,experiment).*strcmp(T.used_as_data,'yes').*(T.retract==0).*((contains(T.brain_region,region{i_region}))).*(T.probe_fraction>0);%
    else
        ind_session=(T.delay_max==1.5).*(T.delay_min==0.3).*strcmp(T.experiment,experiment).*strcmp(T.used_as_data,'yes').*(T.retract==0).*((contains(T.brain_region,region{i_region})));%.*(T.probe_fraction>0)
    end
    %     ind_session=strcmp(T.animal,whichAnimal).*strcmp(T.used_as_data,'yes').*(strcmp(T.date,'2019/11/30'));
    
    ind_session=find(ind_session);
    n_session=length(ind_session);
    animal_unique=unique(T.animal(ind_session));
    n_animal=length(animal_unique);
    n_site=zeros(2,1);%store n sites for soma and terminal
    n_datapoint=zeros(2,1);%store n datapoints for soma and terminal
    for i=1:n_session%this loop to know how many data points are included
        if strcmp( T.brain_region(ind_session(i)),'bilateral SC')
            n_datapoint(1)=n_datapoint(1)+2;%soma
            n_datapoint(2)=n_datapoint(2);%terminal
        else
            n_datapoint(1)=n_datapoint(1)+1;
            n_datapoint(2)=n_datapoint(2)+1;
        end
    end
    for i=1:n_animal%this loop is to find how many sites are included in the dataset
        ind_row=find(strcmp(T.animal,animal_unique{i}));
        if strcmp( T.brain_region(ind_row),'bilateral SC')
            n_site(1)=n_site(1)+2;%soma
            n_site(2)=n_site(2);%terminal
        else
            n_site(1)=n_site(1)+1;
            n_site(2)=n_site(2)+1;
        end
    end
    [dff_aligned_cat_animal,trialType_cat_animal,behEvent_aligned_cat_animal,licking_aligned_cat_animal]=deal(cell(2,n_site(1)));%two fibers so two cells to store,for multiple sessions, store 1D--1-soma,2-terminal data,2D--sites; note soma sites are more than terminal sites
    [dff_aligned_by_session,trialType_by_session,behEvent_aligned_by_session,licking_aligned_by_session]=deal(cell(2,n_datapoint(1)));%2D-data points, each means one site one session
    for nfiber=1:length(fiberstr)%here means fiber bilaterally implanted,1-soma, 2-terminal
        %store AUC and p-value
        fileNameAUCrecent=strcat(savepath,filesep,experiment,'-recentAUCfile');
        dataInfoStr=strcat(experiment,'-',regionstr{i_region},'-',fiberstr{nfiber});
        dataProcessStr=strcat('-grouped by ',selectivitystr{i_selectivity},'-binCentered-size',num2str(binsize),'-shuffle',num2str(nshuffle));
        fileNameAUC=strcat(savepath,filesep,dataInfoStr,dataProcessStr,'-n',num2str(n_datapoint(nfiber)),'-sites',num2str(n_site(nfiber)),'-AUC.mat');
        datapointNameCell=cell(n_datapoint(nfiber),1);
        siteNameCell=[];%later will be struct with field 'nameSite','nameCases'
        fileNameTtest=strcat(savepath,filesep,experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-binCentered-size',num2str(binsize),'-n',num2str(n_datapoint(nfiber)),'-sites',num2str(n_site(nfiber)),'-Ttest.mat');
        pTtestcell=cell(n_datapoint(nfiber),4);
        pTtestcellSite=cell(n_site(nfiber),4);
        %% go through each session, for each data point/site, calculate AUC and p
        ind_datapoint=1;
        ind_site=1;
        for i=1:n_session%遍历满足条件的session并合并
            if nfiber==1%soma
                if strcmp(T.brain_region(ind_session(i)),'left SC')
                    fiberSide={'left'};
                    ind_fiber=[1];
                elseif strcmp(T.brain_region(ind_session(i)),'right SC')
                    fiberSide={'right'};
                    ind_fiber=[2];
                elseif strcmp(T.brain_region(ind_session(i)),'bilateral SC')
                    fiberSide={'left','right'};
                    ind_fiber=[1,2];
                end
            elseif nfiber==2%terminal
                if strcmp(T.brain_region(ind_session(i)),'left SC')
                    fiberSide={'right'};
                    ind_fiber=[2];
                elseif strcmp(T.brain_region(ind_session(i)),'right SC')
                    fiberSide={'left'};
                    ind_fiber=[1];
                elseif strcmp(T.brain_region(ind_session(i)),'bilateral SC')
                    fiberSide=[];
                    ind_fiber=[];
                end
            end
            sessiondate=datevec(T.date(ind_session(i)),'yyyy/mm/dd');
            formatOut = 'yyyymmdd';
            rootpath=strcat('H:\FP\',T.animal(ind_session(i)),'_',datestr(sessiondate,formatOut));%根据总结文件找到对应session的文件夹
            rootpath=rootpath{1};%change from cell array to char array
            files = dir(strcat(rootpath,'\*Virables.mat'));
            if length(files)==1
                behaviorFile=files(1).name;
            else
                warning('Multiple beh files');
            end
            cd(rootpath);
            FrameInfo = dir('*.log');
            fileID = fopen(FrameInfo.name);
            C=textscan(fileID,'%d %d','HeaderLines',16);
            fclose(fileID);
            TrialCount=C{1,1};
            TrialStart_FrameCount=C{1,2};
            nTrial = length(TrialCount);
            %     ImagingSetup
            ImagingSetup=FrameInfo.name(1:end-4);
            disp(ImagingSetup);%in command window display a value
            ImagingSetup(ImagingSetup=='_')='-';
            ImagingSetup_L=sprintf('%s_Left', ImagingSetup);
            ImagingSetup_R=sprintf('%s_Right', ImagingSetup);
            ImagingSetup_Ctrl=sprintf('%s_410Ctrl', ImagingSetup);
            
            %% compute dff or just load
            if exist('dff_temp.mat','file')
                load('dff_temp');
            else
                dff=fGetFPdff(rootpath);
            end
            %% load Beh mat data and extract behavior events
            load(behaviorFile);%load behavior data
            [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime(Data_extract,round(double(TrialStart_FrameCount')/2),frT);
            for i_fiber=1:length(fiberSide)
                %name of this data point
                tempdatapointNameCell=strsplit(rootpath,filesep);
                datapointNameCell{ind_datapoint,1}=strcat(tempdatapointNameCell{end},'_',fiberSide{i_fiber});
                %name this site
                tempSiteName=strcat(T.animal(ind_session(i)),'-',fiberSide{i_fiber});
                tempNSite=length(siteNameCell);
                if isfield(siteNameCell,'nameSite')
                    tempIndSite=arrayfun(@(x) strcmp(x.nameSite,tempSiteName),siteNameCell);
                else
                    tempIndSite=0;
                end
                if sum(tempIndSite)>0
                    ind_site=find(tempIndSite);
                else
                    ind_site=tempNSite+1;%create a new name in a new empty cell
                    siteNameCell(ind_site).nameSite=tempSiteName;
                end
                %add case name to site name
                if isfield(siteNameCell,'nameCases')
                    siteNameCell(ind_site).nameCases{length(siteNameCell(ind_site).nameCases)+1}=datapointNameCell{ind_datapoint,1};
                else
                    siteNameCell(ind_site).nameCases=cell(1,1);
                    siteNameCell(ind_site).nameCases{1}=datapointNameCell{ind_datapoint,1};
                end
                %get trial type
                [trialType,behrule,trialTypeStr] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix',fiberSide{i_fiber});%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
                trialType_cat_animal{nfiber,ind_site}=cat(3,trialType_cat_animal{nfiber,ind_site},trialType);
                trialType_by_session{nfiber,ind_datapoint}=trialType;
                %% concatenate dff, trial type maxtrix and behavior event, etc.
                for ndff=[1] %plot for each dff, see which is better and test whether 410 signals have similar trend
                    if strcmp(behEventAlign,'delay onset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
                        [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff{1,ndff}(ind_fiber(i_fiber),:), behEventFrameIndex,  frameNum );
                    else
                        [ dff_aligned, behEvent_aligned, licking_aligned ] = fAlignSigBehEvent( dff{1,ndff}(ind_fiber(i_fiber),:), behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
                    end
                    dff_aligned_cat_animal{nfiber,ind_site}=cat(1,dff_aligned_cat_animal{nfiber,ind_site},dff_aligned);
                    dff_aligned_by_session{nfiber,ind_datapoint}=dff_aligned;
                    if strcmp(fiberSide{i_fiber},'right')
                        [ behEvent_aligned.lickFirst_left, behEvent_aligned.lickFirst_right]=deal(behEvent_aligned.lickFirst_right,behEvent_aligned.lickFirst_left);%swarp value
                        [licking_aligned.leftLick,licking_aligned.rightLick]=deal(licking_aligned.rightLick,licking_aligned.leftLick);
                    end
                    behEvent_aligned_cat_animal{nfiber,ind_site}=fMergeStruct(behEvent_aligned_cat_animal{nfiber,ind_site},behEvent_aligned);
                    licking_aligned_cat_animal{nfiber,ind_site}=fMergeStruct(licking_aligned_cat_animal{nfiber,ind_site},licking_aligned);
                    behEvent_aligned_by_session{nfiber,ind_datapoint}=behEvent_aligned;
                    licking_aligned_by_session{nfiber,ind_datapoint}=licking_aligned;
                end
                ind_datapoint=ind_datapoint+1;
            end
        end
        indSig=false(size(datapointNameCell));
        for i=1:length(indSig)
            temp=cellfun(@(x) strcmp(x,datapointNameCell{i}),nameSig);
            if sum(temp)>0
                indSig(i)=true;
            end
        end
        
        %         [neuralActivity_by_site,neuralActivity_by_site_raw]=fDFFbyTrialType(dff_aligned_cat_animal,trialType_cat_animal,'z-score');
        %         [neuralActivity_by_datapoint,neuralActivity_by_datapoint_raw]=fDFFbyTrialType(dff_aligned_by_session,trialType_by_session,'z-score');
        [neuralActivity_by_site,neuralActivity_by_site_raw,celldata_by_site]=fDFFbyTrialType(dff_aligned_cat_animal,trialType_cat_animal,'minmax');
        [neuralActivity_by_datapoint,neuralActivity_by_datapoint_raw,celldata_by_datapoint]=fDFFbyTrialType(dff_aligned_by_session(:,indSig),trialType_by_session(:,indSig),'minmax');
        neuralActivity2Plot=cell(1,2);
        neuralActivity2Plot{1}=neuralActivity_by_site;
        neuralActivity2Plot{2}=neuralActivity_by_datapoint;
        
        %% plot raster and mean trace for large session
        for ndff=[1]         %  titlestr=strcat(animal_name,'-',fiberstr{nfiber},'-',titlestr);
            if isempty(neuralActivity2Plot)%if no data, then end this loop
                break;
            end
            %         fig=figure;%plot raster
            %         set(gcf, 'position', [0 400 1400 300]);
            figMeanTrace=figure;%plot mean trace
            set(gcf, 'position', [0 0 1400 600]);
            if contains(trialTypeStr,'combineCorErr')
                titlestr={'do','Miss','Violation'};
            else
                titlestr={'Correct','Error','Miss','Violation'};
            end
            titlestr_row={'-by site','-by session'};
            %             ColLimit = prctile(dff_aligned_cat{nfiber}',98,'all');%here data is a matrix rather than a vector, so use 'all' to find the prctile of all data（ver2018b)
            for figrow=1:2
                for nResult=1:size(neuralActivity2Plot{figrow},1) %4 column(correct/error/miss/violation),companied with 4 lick raster
                    %% plot raster plot and PSTH
                    for nStim=1:size(neuralActivity2Plot{figrow},2) %for each stimulus%[1,6]%just 2 end trials
                        if nResult==(size(neuralActivity2Plot{figrow},1)-1) && size(neuralActivity2Plot{figrow},2)==2%miss trials && grouped with choice rather than stimuli
                            color_mean_trace={[0,0,0],[0,0,0]};
                        elseif size(neuralActivity2Plot{figrow},2)==6%using stimuli to group trials
                            color_mean_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};%left-low--right-high;since cat sessions together, color don't care rule, just using the standard condition
                        elseif size(neuralActivity2Plot{figrow},2)==4%using difficulty to group trials
                            color_mean_trace={[0 0 1],[0.3 0 0.7],[0.7 0 0.3],[1 0 0]};%left-low--right-high;since cat sessions together, color don't care rule, just using the standard condition
                        elseif size(neuralActivity2Plot{figrow},2)==2%using difficulty to group trials
                            color_mean_trace={[0 0 1],[1 0 0]};
                        else
                            color_mean_trace=fMeanTraceColor(behrule,size(trialType,2));
                        end
                        ts= -frameNumTime(1):frT/1000:frameNumTime(2);
                        %plot mean trace together
                        figure(figMeanTrace);%save mean trace
                        subplot(2,size(neuralActivity2Plot{figrow},1),nResult+(figrow-1)*size(neuralActivity2Plot{figrow},1));
                        neuralActivity=neuralActivity2Plot{figrow}{nResult,nStim};
                        curve_meanTrace(nStim)=fPlotMean_SE(ts,neuralActivity,color_mean_trace{nStim});
                        if ~exist('yrange','var')
                            y_lim=get(gca,'Ylim');
                        else
                            y_lim=yrange;
                        end
                        %                     set(gca,'Ylim',yrange);%for comparing different cell type
                    end
                    plot([0,0],[y_lim(1),y_lim(2)],'k-');%align to a behavior event
                    if strcmp(behEventAlign,'stim onset')
                        plot(stimDur*ones(1,2),[y_lim(1),y_lim(2)],'k-');%if align to stimOnset, also plot delay onset
                    elseif strcmp(behEventAlign,'delay onset')
                        plot(-stimDur*ones(1,2),[y_lim(1),y_lim(2)],'k-');%if align to stimOnset, also plot delay onset
                    end
                    xlabel(['Time (s) from ',behEventAlign]);
                    title(strcat(titlestr{nResult},titlestr_row{figrow}));
                    set(gca,'xtick',[-floor(frameNumTime(1)):0.5:frameNumTime(2)],'xticklabel',[-floor(frameNumTime(1)):0.5:frameNumTime(2)]);
                    %   ylim([nanmin(nanmin(neuralActivityGrouped)),nanmax(nanmax(neuralActivityGrouped))]);
                    set(gca,'FontName','Arial','FontSize',14);
                    box off;
                    if nResult==1
                        ylabel('normalized \it\DeltaF/F');
                    end
                end
            end
            figure(figMeanTrace);
            %text label
            subplot(2,size(neuralActivity2Plot{figrow},1),3);%miss plot has most space
            if contains(trialTypeStr,'stimuli')
                h=legend(curve_meanTrace(:),{'ipsi easy','','ipsi hard','contra hard','','contra easy'},'Location','best');
            elseif contains(trialTypeStr,'difficulty')
                h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','contra hard','contra easy'},'Location','best');
            elseif contains(trialTypeStr,'first lick')
                h=legend(curve_meanTrace(:),{'ipsi lick first','contra lick first'},'Location','best');
            elseif contains(trialTypeStr,'sensory')
                h=legend(curve_meanTrace(:),{'low click rate','high click rate'},'Location','best');
            end
            set(h,'box','off');
            text(0,1,strcat('n=',num2str(n_datapoint),' datapoints',num2str(n_session),' sessions,',num2str(n_site(nfiber)),' site,',num2str(n_animal),' animal'),'FontName','Arial','FontSize',14,'Unit','Normalized');
            
            %plot cases
            if strcmp(selectivitystr{i_selectivity},'choice')
                color_cases_trace={[0 0 1],[1 0 0];[0.3 0 0.7],[0.7 0 0.3]};%1st row correct, 2nd row error
                color_diff={[0,0,0];[0,0,0]};
            elseif strcmp(selectivitystr{i_selectivity},'stimuli')
                color_cases_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0];[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};
                color_diff={[0,0,0],[0.4,0.4,0.4],[0.7,0.7,0.7];[0,0,0],[0.4,0.4,0.4],[0.7,0.7,0.7]};
            elseif strcmp(selectivitystr{i_selectivity},'difficulty')
                color_cases_trace={[0 0 1],[0.3 0 0.7],[0.7 0 0.3],[1 0 0];[0 0 1],[0.3 0 0.7],[0.7 0 0.3],[1 0 0]};
                color_diff={[0,0,0],[0.5,0.5,0.5];[0,0,0],[0.5,0.5,0.5]};
            end
            resultStr={'cor','err','miss','vio'};
            %             figDiffMeanTrace=figure;
            
            %plot by datapoint
            col_figcase=6;
            iResult=1;
            %neuralActivity{nResult,nStim}(ncase,timepoint)
            [figCases,curve_caseTrace] = fPlotCases(ts,neuralActivity_by_datapoint_raw(iResult,:),col_figcase,color_cases_trace,behEventAlign,frameNumTime,selectivitystr{i_selectivity},stimDur);
            neuralActivity_by_datapoint_diff=fContraIpsiDiff(neuralActivity_by_datapoint_raw(iResult,:),2);
            [figCasesDiff,curve_caseDiffTrace] = fPlotCases(ts,neuralActivity_by_datapoint_diff,col_figcase,color_diff,behEventAlign,frameNumTime,selectivitystr{i_selectivity},stimDur);
            %plot by site
            col_figsite=4;
            [figSite,curve_siteTrace] = fPlotCases(ts,neuralActivity_by_site_raw(iResult,:),col_figsite,color_cases_trace,behEventAlign,frameNumTime,selectivitystr{i_selectivity},stimDur);
            neuralActivity_by_site_diff=fContraIpsiDiff(neuralActivity_by_site_raw(iResult,:),2);
            [figSiteDiff,curve_siteDiffTrace] = fPlotCases(ts,neuralActivity_by_site_diff,col_figsite,color_diff,behEventAlign,frameNumTime,selectivitystr{i_selectivity},stimDur);
            %                 normalizedDFF_by_datapoint_diff=cellfun(@(x) fLocalNormalized(x, 'minmax'),neuralActivity_by_datapoint_diff,'UniformOutput',0);
            %                 normalizedDFF_by_site_diff=cellfun(@(x) fLocalNormalized(x, 'minmax'),neuralActivity_by_site_diff,'UniformOutput',0);
            normalizedDFF_by_datapoint_diff=cellfun(@(x) fLocalNormalized(x, 'z-score'),neuralActivity_by_datapoint_diff,'UniformOutput',0);
            normalizedDFF_by_site_diff=cellfun(@(x) fLocalNormalized(x, 'z-score'),neuralActivity_by_site_diff,'UniformOutput',0);
            %                 %plot diff mean trace
            %                 figure(figDiffMeanTrace);
            %                 subplot(2,length(resultStr),iResult);%by session
            %                 for nStim=1:size(normalizedDFF_by_datapoint_diff,2)
            %                     curve_diffMeanTrace(nStim)=fPlotMean_SE(ts,normalizedDFF_by_datapoint_diff{nStim},color_diff{1,nStim});hold on;
            %                 end
            %                 title(strcat(titlestr{iResult},'by session'));
            %                 subplot(2,length(resultStr),iResult+length(resultStr));%by site
            %                 for nStim=1:size(normalizedDFF_by_site_diff,2)
            %                     curve_diffMeanTrace(nStim)=fPlotMean_SE(ts,normalizedDFF_by_site_diff{nStim},color_diff{1,nStim});hold on;
            %                 end
            %                 title(strcat(titlestr{iResult},'by site'));
            %                 for irow=1:2
            %                     subplot(2,length(resultStr),(irow-1)*length(resultStr)+iResult);
            %                     plot([0,0],[y_lim(1),y_lim(2)],'k-');%align to a behavior event
            %                     if strcmp(behEventAlign,'stim onset')
            %                         plot(stimDur*ones(1,2),[y_lim(1),y_lim(2)],'k-');%if align to stimOnset, also plot delay onset
            %                     elseif strcmp(behEventAlign,'delay onset')
            %                         plot(-stimDur*ones(1,2),[y_lim(1),y_lim(2)],'k-');%if align to stimOnset, also plot delay onset
            %                     end
            %                     xlabel(['Time (s) from ',behEventAlign]);
            %                     set(gca,'xtick',[-floor(frameNumTime(1)):0.5:frameNumTime(2)],'xticklabel',[-floor(frameNumTime(1)):0.5:frameNumTime(2)]);
            %                     %   ylim([nanmin(nanmin(neuralActivityGrouped)),nanmax(nanmax(neuralActivityGrouped))]);
            %                     set(gca,'FontName','Arial','FontSize',14);
            %                     box off;
            %                     if nResult==1
            %                         ylabel('normalized \it\DeltaF/F (contra-ipsi)');
            %                     end
            %                 end
            
            
        end
        %save figure
        if n_animal==1%only one animal
            if strcmp(masklick,'yes')
                suptitle(strcat(animal_name,'-',regionstr{i_region},'-',fiberstr{nfiber},'-',dffstr{ndff},'-masklick'));
            else
                suptitle(strcat(animal_name,'-',regionstr{i_region},'-',fiberstr{nfiber},'-',dffstr{ndff}));
            end
            saveas(figMeanTrace,[savepath,filesep,animal_name,'-',experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-sort ',behEventSort,'-masklick-',masklick,'-n',num2str(n_datapoint(nfiber)),'-mean trace with end point dff.pdf'],'pdf');
            %             saveas(fig,[savepath,filesep,animal_name,'-',experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-sort ',behEventSort,'-masklick-',masklick,'-raster.pdf'],'pdf');
        else%use common feature, e.g. experiment/group name
            if strcmp(masklick,'yes')
                suptitle(strcat(experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-',dffstr{ndff},'-masklick'));
            else
                suptitle(strcat(experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-',dffstr{ndff}));
            end
            saveas(figMeanTrace,[savepath,filesep,experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-sort ',behEventSort,'-masklick-',masklick,'-n',num2str(n_datapoint(nfiber)),'-mean trace with end point dff.pdf'],'pdf');
            %                 saveas(figCases,[savepath,filesep,experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-masklick-',masklick,'-individual cases.pdf'],'pdf');
            %                 saveas(figSite,[savepath,filesep,experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-masklick-',masklick,'-individual cases by site.pdf'],'pdf');
            %                 saveas(figCasesDiff,[savepath,filesep,experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-masklick-',masklick,'-individual cases difference.pdf'],'pdf');
            %                 saveas(figSiteDiff,[savepath,filesep,experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-masklick-',masklick,'-individual cases difference by site.pdf'],'pdf');
        end
    end
    
    n_datapointCell{i_celltype}=n_datapoint;
    n_siteCell{i_celltype}=n_site;
    
end

%}

%% auxiliary function
%to merge data points together to a larger session
function s=fMergeStruct(a,b)
if isempty(a)
    s=b;
else
    f = fieldnames(b);
    for j=1:length(f)
        s.(f{j})=cat(2,a.(f{j}),b.(f{j}));
    end
end
end

%to transform cell to mat, and fill empty cell with nan
function out=fCell2Mat(in)
sizein=size(in);%usually should be 2-d
out=zeros(sizein);
for i=1:sizein(1)
    for j=1:sizein(2)
        if isempty(in{i,j})
            out(i,j)=nan;
        else
            out(i,j)=in{i,j};
        end
    end
end
if (length(size(in))-sum(size(in)==1))==1%1d- vector
    out(isnan(out))=[];
end
end

%to use trialtype and dff_aligned to get normalized plotting data
function [matdata,matdataraw,celldata]=fDFFbyTrialType(dff,trialtypecell,method)
trialtype=trialtypecell{1};
matdata=cell(size(trialtype,1),size(trialtype,2));
matdataraw=cell(size(trialtype,1),size(trialtype,2));
celldata=cell(size(trialtype,1),size(trialtype,2),length(dff));%1d-nResult,2d-nStim,3d-nsession
for nResult=1:size(trialtype,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
    for nStim=1:size(trialtype,2) %for each stimulus%[1,6]%just 2 end trials
        matdataraw{nResult,nStim}=zeros(length(dff),size(dff{1},2));
        for idatapoint=1:length(dff)
            selectedTrialInd=trialtypecell{1,idatapoint}(nResult,nStim,:);
            selectedTrialInd=logical(squeeze(selectedTrialInd))';
            rawdata=dff{1,idatapoint}(selectedTrialInd,:);
            matdataraw{nResult,nStim}(idatapoint,:)=nanmean(rawdata,1);
            celldata{nResult,nStim,idatapoint}=rawdata;
        end
    end
end
%normalize for each datapoint
if strcmp(method,'minmax')
    disp('normalized by (x-min)/(max-min)');
elseif strcmp(method,'z-score')
    disp('normalized by z-score');
end
maxDataCell=cellfun(@(x) max(x,[],2), matdataraw,'UniformOutput',false);
minDataCell=cellfun(@(x) min(x,[],2), matdataraw,'UniformOutput',false);
for idatapoint=1:length(dff)
    maxDataMat=cellfun(@(x) x(idatapoint),maxDataCell);
    minDataMat=cellfun(@(x) x(idatapoint),minDataCell);
    maxDataVal=max(max(maxDataMat));
    minDataVal=min(min(minDataMat));
    
    for nResult=1:size(trialtype,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
        for nStim=1:size(trialtype,2) %for each stimulus%[1,6]%just 2 end trials
            if strcmp(method,'minmax')
                matdata{nResult,nStim}(idatapoint,:)=(matdataraw{nResult,nStim}(idatapoint,:)-minDataVal)/(maxDataVal-minDataVal);%method 1
            elseif strcmp(method,'z-score')
                matdata{nResult,nStim}(idatapoint,:)=(matdataraw{nResult,nStim}(idatapoint,:)-nanmean(nanmean(dff{idatapoint})))/nanstd(dff{idatapoint},0,'all');%method 2, z-score
            end
        end
    end
end
end

%from cases to get diff
function [datadiff]= fContraIpsiDiff(dataRaw,d)
%d indicates which dimesion to substract
if d ==2
    datadiff=cell(size(dataRaw,1),size(dataRaw,2)/2);
    for nResult=1:size(dataRaw,1)
        for nStim=1:size(dataRaw,2)/2
            datadiff{nResult,nStim}=dataRaw{nResult,size(dataRaw,2)+1-nStim}-dataRaw{nResult,nStim};
        end
    end
elseif d ==1
    datadiff=cell(size(dataRaw,1)/2,size(dataRaw,2));
    for nResult=1:size(dataRaw,1)/2
        for nStim=1:size(dataRaw,2)
            datadiff{nResult,nStim}=dataRaw{size(dataRaw,2)+1-nResult,nStim}-dataRaw{nResult,nStim};
        end
    end
end
end

%use to plot individual cases
function [figCases,curve_caseTrace] = fPlotCases(ts,neuralActivity,col_figcase,color_cases_trace,behEventAlign,frameNumTime,selectivitystr,stimDur)
n_cases=size(neuralActivity{1,1},1);
row_figcase=ceil(n_cases/col_figcase);
figCases=figure;
set(gcf, 'paperPosition', [0 0 2*col_figcase 2*row_figcase]);
figure(figCases);
for i_case=1:n_cases
    for nResult=1:size(neuralActivity,1)
        for nStim=1:size(neuralActivity,2)
            subplot(row_figcase,col_figcase,i_case);hold on;
            curve_caseTrace(nStim)=plot(ts,neuralActivity{nResult,nStim}(i_case,:),'color',color_cases_trace{nResult,nStim});
            if ~exist('yrange','var')
                y_lim=get(gca,'Ylim');
            else
                y_lim=yrange;
            end
            %   set(gca,'Ylim',yrange);%for comparing different cell type
            plot([0,0],[y_lim(1),y_lim(2)],'k-');%align to a behavior event
            if strcmp(behEventAlign,'stim onset')
                plot(stimDur*ones(1,2),[y_lim(1),y_lim(2)],'k-');%if align to stimOnset, also plot delay onset
            elseif strcmp(behEventAlign,'delay onset')
                plot(-stimDur*ones(1,2),[y_lim(1),y_lim(2)],'k-');%if align to stimOnset, also plot delay onset
            end
            if i_case>(n_cases-col_figcase)% the last row
                xlabel(['Time (s) from ',behEventAlign]);
            end
            set(gca,'xtick',[-floor(frameNumTime(1)):0.5:frameNumTime(2)],'xticklabel',[-floor(frameNumTime(1)):0.5:frameNumTime(2)]);
            xlim([-frameNumTime(1),frameNumTime(2)]);
            %   ylim([nanmin(nanmin(neuralActivityGrouped)),nanmax(nanmax(neuralActivityGrouped))]);
            set(gca,'FontName','Arial','FontSize',12);
            box off;
            if mod(i_case,col_figcase)==1 % the 1st col
                ylabel('\it\DeltaF/F');
            end
            if i_case==n_cases
                switch selectivitystr
                    case 'choice'
                        legend(curve_caseTrace(:),'ipsi correct','contra correct');
                        %                             legend(curve_caseTrace(:),'ipsi correct','contra correct','ipsi error','contra error');
                    case 'difficulty'
                        legend(curve_caseTrace(:),'ipsi easy','ipsi hard','contra hard','contra easy');
                    case 'stimuli'
                        legend(curve_caseTrace(:),'','ipsi','','','contra','');
                end
            end
            
        end
    end
end
end