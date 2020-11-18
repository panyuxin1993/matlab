%This code choose some sessions and plot PSTH of each cases; 
%both show individual sessions and animals
% calculate prefered/prefered-nonprefered values and compare this value of
% different difficulties using 1)AUC 2)ANOVA
% dbstop if error
clear;
close all;
% %% import .csv to decide which session to include
% TepochAUC=readtable('F:\FP\summary\EpochAUC.csv');
% TepochAUCsig=TepochAUC(TepochAUC.pAUClate<0.025 | TepochAUC.pAUClate>0.975,:);
% nameSig=TepochAUCsig.DatapointName;
%% import data summay file
savepath='F:\FP\summary';
[num,txt,raw] =xlsread('F:\FP\fiber photometry data summary.xlsx');%criteria to choose sessions come from this file
T=cell2table(raw(2:end,1:12));
T.Properties.VariableNames=strrep(raw(1,1:12),' ','_');%table variable name can't have ' ',so replace them
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
iPrefered=2;
preferedStr={'prefered','prefered-non'};
iAUCstr=1;
AUCstr={'hard-easy','hardest-easiest'};
% yrange=[-0.01,0.06];
% const settings of fiber photometry
Session470LEDstartFrame=1;% usually 205
Session410LEDstartFrame=2;
FrameRate=40;
FrameTime=1000/FrameRate;
frT=FrameTime*2;%2 channel, so framerate and frametime should be half
frameNum=double(round(frameNumTime*1000/frT));
ts= -frameNumTime(1):frT/1000:frameNumTime(2);
%parrameters for calculating AUC etc.
binsize=3;
binstep=1;
nshuffle=1000;
p_sepOnset=0.05;%the p-value used to determine when separation happen
%set the criteria to choose sessions
i_region=3;%1-bilateral, 2-unilateral*********variable**************
region={'bilateral','t SC','SC'};%{'bilateral SC','SC'}
regionstr={'bilateral','unilateral','mixed'};
whichAnimal='pyx241';
celltype={'SC vgat'};
% celltype={'ALM terminal'};
n_siteCell=cell(1,2);%each cell store a cell type
n_datapointCell=cell(1,2);%each cell store a cell type
for i_celltype=1:length(celltype)
    %select sessions for summary analysis
    experiment=celltype{i_celltype};
    ind_session=strcmp(T.experiment,experiment).*contains(T.used_as_data,'yes').*(T.retract==0).*((contains(T.brain_region,region{i_region}))).*(contains(T.manipulation,'opto'));
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
        dataProcessStr=strcat('-grouped by ',selectivitystr{i_selectivity},'-binCentered-size',num2str(binsize),'-shuffle',num2str(nshuffle),preferedStr{iPrefered},'-AUC',AUCstr{iAUCstr});
        fileNameDifficultyAUC=strcat(savepath,filesep,dataInfoStr,dataProcessStr,'-n',num2str(n_datapoint(nfiber)),'-sites',num2str(n_site(nfiber)),'-difficulty.mat');
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
            rootpath=strcat('F:\FP\',T.animal(ind_session(i)),'_',datestr(sessiondate,formatOut));%根据总结文件找到对应session的文件夹
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
                tempdatapointNameCell=strsplit(rootpath,filesep);
%                 tempName= strcat(tempdatapointNameCell{end},'_',fiberSide{i_fiber});
%                 % include only significant late activities sites
%                 temp=cellfun(@(x) strcmp(x,tempName),nameSig);
%                 if sum(temp)==0
%                     continue;
%                 end
                %name of this data point
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
                [trialType,behrule,trialTypeStr] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix',fiberSide{i_fiber},'divideCorErr','divideOpto');%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
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
%         indSig=false(size(datapointNameCell));
%         for i=1:length(indSig)
%             temp=cellfun(@(x) strcmp(x,datapointNameCell{i}),nameSig);
%             if sum(temp)>0
%                 indSig(i)=true;
%             end
%         end
%         datapointNameCellSig=datapointNameCell(indSig);

%         [neuralActivity_by_site,neuralActivity_by_site_raw]=fDFFbyTrialType(dff_aligned_cat_animal,trialType_cat_animal,'z-score');
%         [neuralActivity_by_datapoint,neuralActivity_by_datapoint_raw]=fDFFbyTrialType(dff_aligned_by_session,trialType_by_session,'z-score');
%         %rule out empty cell that may be due to data exclusion
%         dff_aligned_cat_animal=fRefine(dff_aligned_cat_animal);
%         dff_aligned_by_session=fRefine(dff_aligned_by_session);
%         trialType_cat_animal=fR
        [neuralActivity_by_site,neuralActivity_by_site_raw,celldata_by_site]=fDFFbyTrialType(dff_aligned_cat_animal,trialType_cat_animal,'minmax');
        [neuralActivity_by_datapoint,neuralActivity_by_datapoint_raw,celldata_by_datapoint]=fDFFbyTrialType(dff_aligned_by_session,trialType_by_session,'minmax');
        
        if size(celldata_by_datapoint,4)==2  && size(celldata_by_datapoint,2)==2
            temp_datapoint=cell(size(celldata_by_datapoint,1),size(celldata_by_datapoint,2)*2,size(celldata_by_datapoint,3));
            temp_site =cell(size(celldata_by_site,1),size(celldata_by_site,2)*2,size(celldata_by_site,3));
            temp_datapoint(:,1,:)=celldata_by_datapoint(:,1,:,1);
            temp_datapoint(:,2,:)=celldata_by_datapoint(:,2,:,1);
            temp_datapoint(:,3,:)=celldata_by_datapoint(:,2,:,2);
            temp_datapoint(:,4,:)=celldata_by_datapoint(:,1,:,2);
            celldata_by_datapoint=temp_datapoint;
            temp_site(:,1,:)=celldata_by_site(:,1,:,1);
            temp_site(:,2,:)=celldata_by_site(:,2,:,1);
            temp_site(:,3,:)=celldata_by_site(:,2,:,2);
            temp_site(:,4,:)=celldata_by_site(:,1,:,2);
            celldata_by_site=temp_site;
        end
        
        neuralActivity2Plot=cell(1,2);
        neuralActivity2Plot{1}=celldata_by_site;
        neuralActivity2Plot{2}=celldata_by_datapoint;
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
                        elseif size(neuralActivity2Plot{figrow},2)==4%with opto trials
                            color_mean_trace={[0 0 1],[0 0.5 1],[1 0.5 0],[1 0 0]};%left-low--right-high;since cat sessions together, color don't care rule, just using the standard condition
                        elseif size(neuralActivity2Plot{figrow},2)==2%using difficulty to group trials
                            color_mean_trace={[0 0 1],[1 0 0]};
                        else
                            color_mean_trace=fMeanTraceColor(behrule,size(trialType,2));
                        end
                        
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
                h=legend(curve_meanTrace(:),{'ipsi lick first ctrl','ipsi lick first opto','contra lick first opto','contra lick first ctrl'},'Location','best');
            elseif contains(trialTypeStr,'sensory')
                h=legend(curve_meanTrace(:),{'low click rate','high click rate'},'Location','best');
            end
            set(h,'box','off');
            text(0,1,strcat('n=',num2str(n_datapoint),' datapoints',num2str(n_session),' sessions,',num2str(n_site(nfiber)),' site,',num2str(n_animal),' animal'),'FontName','Arial','FontSize',14,'Unit','Normalized');
            
            %plot cases
            if strcmp(selectivitystr{i_selectivity},'choice')
                color_cases_trace={[0 0 1],[1 0 0];[0.3 0 0.7],[0.7 0 0.3]};%1st row correct, 2nd row error
                color_diff={[0,0,0];[0,0,0]};
            elseif size(celldata_by_datapoint,2)==6 %strcmp(selectivitystr{i_selectivity},'stimuli')
                color_cases_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0];[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};
                color_diff={[0,0,0],[0.4,0.4,0.4],[0.7,0.7,0.7];[0,0,0],[0.4,0.4,0.4],[0.7,0.7,0.7]};
            elseif size(celldata_by_datapoint,2)==4 %strcmp(selectivitystr{i_selectivity},'difficulty')
                color_cases_trace={[0 0 1],[0.3 0 0.7],[0.7 0 0.3],[1 0 0];[0 0 1],[0.3 0 0.7],[0.7 0 0.3],[1 0 0]};
                color_diff={[0,0,0],[0.5,0.5,0.5];[0,0,0],[0.5,0.5,0.5]};
            end
            resultStr={'cor','err','miss','vio'};
%             figDiffMeanTrace=figure;
            
            for iResult=1%:length(resultStr) 
                %find prefer direction
                diffdata_by_datapoint=cell(size(celldata_by_datapoint,1),size(celldata_by_datapoint,2)/2,size(celldata_by_datapoint,3));
                diffdata_by_site=cell(size(celldata_by_site,1),size(celldata_by_site,2)/2,size(celldata_by_site,3));
                for iData=1:size(celldata_by_datapoint,3)%for each session
                    cell_by_datapoint=celldata_by_datapoint(iResult,:,iData);%each session, every stimuli
                    meandff=cellfun(@(x) nanmean(nanmean(x)),cell_by_datapoint);
                    ipsimeandff=mean(meandff(1:length(meandff)/2));
                    contrameandff=mean(meandff(length(meandff)/2+1:end));
                    if ipsimeandff>contrameandff
                        nonprefered=cell_by_datapoint(length(cell_by_datapoint)/2+1:end);
                        nonprefered=fliplr(nonprefered);
                        prefered=cell_by_datapoint(1:length(cell_by_datapoint)/2);
                    else
                        nonprefered=cell_by_datapoint(1:length(cell_by_datapoint)/2);
                        prefered=cell_by_datapoint(length(cell_by_datapoint)/2+1:end);
                        prefered=fliplr(prefered);
                    end
                    if strcmp(preferedStr{iPrefered},'prefered-non')
                        datadiff=cellfun(@(x,y) y-nanmean(x,1),nonprefered,prefered,'UniformOutput',false);
                    elseif strcmp(preferedStr{iPrefered},'prefered')
                        datadiff=prefered;
                    end

                    poslabel=2;
                    Structout= fCell2CategoryStruct(datadiff);
                    [StructDiff(iData).label,StructDiff(iData).data]=deal(Structout.label,Structout.data);
                    StructDiff(iData).earlyData = nanmean(Structout.data(:,(frameNumTime(1)-0.5)*1000/frT:(frameNumTime(1))*1000/frT),2);
                    StructDiff(iData).lateData = nanmean(Structout.data(:,(frameNumTime(1)+1)*1000/frT:(frameNumTime(1)+1.5)*1000/frT),2);
                    StructDiff(iData).datapointName=datapointNameCell{iData};
                    fileNameDifficultyAUC_datapoint=strcat(savepath,filesep,dataInfoStr,dataProcessStr,'-',StructDiff(iData).datapointName,dataProcessStr,'-difficultyAUC.mat');%adding datapoint name to the AUC file name
                    if exist(fileNameDifficultyAUC_datapoint,'file')
                        load(fileNameDifficultyAUC_datapoint);
                    end
                    %calculate choiceAUC for different difficulty
                    if ~exist('choiceAUC_diffopto','var') 
                        n_difficulty= size(celldata_by_datapoint,2)/2;
                        for i=1:n_difficulty
                            datachoiceAUC=celldata_by_datapoint(iResult,[i,size(celldata_by_datapoint,2)-i+1],iData);
                            tempStruct=fCell2CategoryStruct(datachoiceAUC);
                            [choiceAUC_diffopto(i,:),pchoiceAUC_diffopto(i,:)] = fMovingAUC(tempStruct.label,tempStruct.data,poslabel,nshuffle,binsize,binstep);
                        end
                    end
                    [StructDiff(iData).choiceAUC_diffopto,StructDiff(iData).pchoiceAUC_diffopto] = deal(choiceAUC_diffopto,pchoiceAUC_diffopto);
                    if length(unique(StructDiff(iData).label))==2
                        if ~exist('optoctrlAUCcase','var')
                            [optoctrlAUCcase,poptoctrlAUCcase] = fMovingAUC(StructDiff(iData).label,StructDiff(iData).data,poslabel,nshuffle,binsize,binstep);
                        end
                        [StructDiff(iData).AUC,StructDiff(iData).pAUC] = deal(optoctrlAUCcase,poptoctrlAUCcase);
                        if ~exist('AUCearly','var')
                            [AUCearly,pAUCearly]=fAUC(StructDiff(iData).label,StructDiff(iData).earlyData,poslabel,nshuffle);
                            [AUClate,pAUClate] = fAUC(StructDiff(iData).label,StructDiff(iData).lateData,poslabel,nshuffle);                 
                        end
                        [StructDiff(iData).AUCearly,StructDiff(iData).pAUCearly,StructDiff(iData).AUClate,StructDiff(iData).pAUClate] = deal(AUCearly,pAUCearly,AUClate,pAUClate);
                         save(fileNameDifficultyAUC_datapoint,'optoctrlAUCcase','poptoctrlAUCcase','AUCearly','pAUCearly','AUClate','pAUClate','choiceAUC_diffopto','pchoiceAUC_diffopto');
                        clearvars  'AUCearly' 'pAUCearly' 'AUClate' 'pAUClate' 'optoctrlAUCcase' 'poptoctrlAUCcase' 'choiceAUC_diffopto' 'pchoiceAUC_diffopto'
                    elseif length(unique(StructDiff(iData).label))>2 %using anova
                        StructDiff(iData).pANOVA1=fMovingANOVA1(StructDiff(iData).data,StructDiff(iData).label,binsize,binstep);
                        save(fileNameDifficultyAUC_datapoint,'choiceAUC_diffopto','pchoiceAUC_diffopto');
                        clearvars 'choiceAUC_diffopto' 'pchoiceAUC_diffopto'
                    else
                        warning([num2str(length(unique(StructDiff(iData).label))),',not 2 category, unable to calculate AUC']);
                    end

                end

                col_figcase=6;
                %plot AUC betweeen difficulty by datapoint
                if length(unique(StructDiff(iData).label))==2
                    [AUCarray{1:length(StructDiff)}]=deal(StructDiff.AUC);
                    [casename{1:length(StructDiff)}]=deal(StructDiff.datapointName);
                    [pAUCcell{1:length(StructDiff)}]=deal(StructDiff.pAUC);
                    sigX=cellfun(@(x) (x<p_sepOnset/2)|(x>1-p_sepOnset/2),pAUCcell,'UniformOutput',0);
                    fPlotAUCCases(ts,AUCarray,col_figcase,casename,sigX);
                elseif length(unique(StructDiff(iData).label))==3
                    [pANOVA{1:length(StructDiff)}]=deal(StructDiff.pANOVA1);
                    sigX=cellfun(@(x) (x<p_sepOnset),pANOVA,'UniformOutput',0);
                end
                clearvars AUCarray casename pAUCcell
                %plot choice AUC of different diffiulty levels
                tempcellAUC=cell(1,size(celldata_by_datapoint,2)/2,size(celldata_by_datapoint,3));
                for i_difficulty=1:size(tempcellAUC,2)
                    for i_case= 1:size(tempcellAUC,3)
                      tempcellAUC{1,i_difficulty, i_case}= StructDiff(i_case).choiceAUC_diffopto(i_difficulty,:);
                    end
                end
                [figAUCCases,curve_AUCcaseTrace] = fPlotCases(ts,tempcellAUC(1,:,:),col_figcase,color_diff,behEventAlign,'choice AUC',frameNumTime,selectivitystr{i_selectivity},stimDur,sigX);
                %plot by datapoint
                %neuralActivity{nResult,nStim}(ncase,timepoint)
                [figCases,curve_caseTrace] = fPlotCases(ts,celldata_by_datapoint(iResult,:,:),col_figcase,color_cases_trace,behEventAlign,'\it\DeltaF/F',frameNumTime,selectivitystr{i_selectivity},stimDur,sigX);
                neuralActivity_by_datapoint_diff=fContraIpsiDiff(celldata_by_datapoint(iResult,:,:),preferedStr{iPrefered});
                [figCasesDiff,curve_caseDiffTrace] = fPlotCases(ts,neuralActivity_by_datapoint_diff,col_figcase,color_diff,behEventAlign,'\it\DeltaF/F(contra-ipsi)',frameNumTime,selectivitystr{i_selectivity},stimDur,sigX);
                %% similar analysis for site
                for iData=1:size(celldata_by_site,3)%for each session
                    cell_by_site=celldata_by_site(iResult,:,iData);%each session, every stimuli
                    meandff=cellfun(@(x) nanmean(nanmean(x)),cell_by_site);
                    ipsimeandff=mean(meandff(1:length(meandff)/2));
                    contrameandff=mean(meandff(length(meandff)/2+1:end));
                    if ipsimeandff>contrameandff
                        nonprefered=cell_by_site(length(cell_by_site)/2+1:end);
                        nonprefered=fliplr(nonprefered);
                        prefered=cell_by_site(1:length(cell_by_site)/2);
                    else
                        nonprefered=cell_by_site(1:length(cell_by_site)/2);
                        prefered=cell_by_site(length(cell_by_site)/2+1:end);
                        prefered=fliplr(prefered);
                    end
                    if strcmp(preferedStr{iPrefered},'prefered-non')
                        datadiff=cellfun(@(x,y) y-nanmean(x,1),nonprefered,prefered,'UniformOutput',false);
                    elseif strcmp(preferedStr{iPrefered},'prefered')
                        datadiff=prefered;
                    end
                    
 
                    StructoutSite= fCell2CategoryStruct(datadiff);
                    [StructDiffSite(iData).label,StructDiffSite(iData).data]=deal(StructoutSite.label,StructoutSite.data);
                    StructDiffSite(iData).earlyData = nanmean(StructoutSite.data(:,(frameNumTime(1)-0.5)*1000/frT:(frameNumTime(1))*1000/frT),2);
                    StructDiffSite(iData).lateData = nanmean(StructoutSite.data(:,(frameNumTime(1)+1)*1000/frT:(frameNumTime(1)+1.5)*1000/frT),2);
                    StructDiffSite(iData).siteName=siteNameCell(iData).nameSite;
                    fileNameDifficultyAUC_site=strcat(savepath,filesep,dataInfoStr,dataProcessStr,'-',StructDiffSite(iData).siteName{1},dataProcessStr,'-difficultyAUC.mat');%adding datapoint name to the AUC file name
                    if exist(fileNameDifficultyAUC_site,'file')
                        load(fileNameDifficultyAUC_site);
                    end
                    %calculate choiceAUC for different difficulty
                    if ~exist('choiceAUC_diffopto','var') 
                        n_difficulty= size(celldata_by_site,2)/2;
                        for i=1:n_difficulty
                            datachoiceAUC=celldata_by_site(iResult,[i,size(celldata_by_datapoint,2)-i+1],iData);
                            tempStruct=fCell2CategoryStruct(datachoiceAUC);
                            [choiceAUC_diffopto(i,:),pchoiceAUC_diffopto(i,:)] = fMovingAUC(tempStruct.label,tempStruct.data,poslabel,nshuffle,binsize,binstep);
                        end
                    end              
                    [StructDiffSite(iData).choiceAUC_diffopto,StructDiffSite(iData).pchoiceAUC_diffopto] = deal(choiceAUC_diffopto,pchoiceAUC_diffopto);
                    if length(unique(StructDiffSite(iData).label))==2
                        if ~exist('optoctrlAUCcase','var')
                            [optoctrlAUCcase,poptoctrlAUCcase] = fMovingAUC(StructDiffSite(iData).label,StructDiffSite(iData).data,poslabel,nshuffle,binsize,binstep);
                        end
                        [StructDiffSite(iData).AUC,StructDiffSite(iData).pAUC] = deal(optoctrlAUCcase,poptoctrlAUCcase);
                        if ~exist('AUCearly','var')
                            [AUCearly,pAUCearly]=fAUC(StructDiffSite(iData).label,StructDiffSite(iData).earlyData,poslabel,nshuffle);
                            [AUClate,pAUClate] = fAUC(StructDiffSite(iData).label,StructDiffSite(iData).lateData,poslabel,nshuffle);                 
                        end
                        [StructDiffSite(iData).AUCearly,StructDiffSite(iData).pAUCearly,StructDiffSite(iData).AUClate,StructDiffSite(iData).pAUClate] = deal(AUCearly,pAUCearly,AUClate,pAUClate);
                         save(fileNameDifficultyAUC_site,'optoctrlAUCcase','poptoctrlAUCcase','AUCearly','pAUCearly','AUClate','pAUClate','choiceAUC_diffopto','pchoiceAUC_diffopto');
                        clearvars  'AUCearly' 'pAUCearly' 'AUClate' 'pAUClate' 'choiceAUC_diffopto' 'pchoiceAUC_diffopto'
                    elseif length(unique(StructDiffSite(iData).label))>2 %using anova
                        StructDiffSite(iData).pANOVA1=fMovingANOVA1(StructDiffSite(iData).data,StructDiffSite(iData).label,binsize,binstep);
                        save(fileNameDifficultyAUC_site,'choiceAUC_diffopto','pchoiceAUC_diffopto');
                        clearvars 'choiceAUC_diffopto' 'pchoiceAUC_diffopto'
                    else
                        warning([num2str(length(unique(StructDiffSite(iData).label))),',not 2 category, unable to calculate AUC']);
                    end

                end
                save(fileNameDifficultyAUC,'StructDiff','StructDiffSite');
                col_figsite=4;
                %plot AUC by site
                if length(unique(StructDiffSite(iData).label))==2
                    [AUCarray{1:length(StructDiffSite)}]=deal(StructDiffSite.AUC);
                    [casename{1:length(StructDiffSite)}]=deal(StructDiffSite.siteName);
                    [pAUCcell{1:length(StructDiffSite)}]=deal(StructDiffSite.pAUC);
                    sigX=cellfun(@(x) (x<p_sepOnset/2)|(x>1-p_sepOnset/2),pAUCcell,'UniformOutput',0);
                    fPlotAUCCases(ts,AUCarray,col_figsite,casename,sigX);
                elseif length(unique(StructDiffSite(iData).label))==3
                    [pANOVA{1:length(StructDiffSite)}]=deal(StructDiffSite.pANOVA1);
                    sigX=cellfun(@(x) (x<p_sepOnset),pANOVA,'UniformOutput',0);
                end
                clearvars AUCarray casename pAUCcell 
                %plot choice AUC of different diffiulty levels by site
                tempcellAUC=cell(1,size(celldata_by_site,2)/2,size(celldata_by_site,3));
                for i_difficulty=1:size(tempcellAUC,2)
                    for i_case= 1:size(tempcellAUC,3)
                      tempcellAUC{1,i_difficulty, i_case}= StructDiffSite(i_case).choiceAUC_diffopto(i_difficulty,:);
                    end
                end
                [figAUCCases,curve_AUCcaseTrace] = fPlotCases(ts,tempcellAUC(1,:,:),col_figcase,color_diff,behEventAlign,'choice AUC',frameNumTime,selectivitystr{i_selectivity},stimDur,sigX);

                %plot by site
                [figSite,curve_siteTrace] = fPlotCases(ts,celldata_by_site(iResult,:,:),col_figsite,color_cases_trace,behEventAlign,'\it\DeltaF/F',frameNumTime,selectivitystr{i_selectivity},stimDur,sigX);
                neuralActivity_by_site_diff=fContraIpsiDiff(celldata_by_site(iResult,:,:),preferedStr{iPrefered});
                [figSiteDiff,curve_siteDiffTrace] = fPlotCases(ts,neuralActivity_by_site_diff,col_figsite,color_diff,behEventAlign,'\it\DeltaF/F(contra-ipsi)',frameNumTime,selectivitystr{i_selectivity},stimDur,sigX);

%                 normalizedDFF_by_datapoint_diff=cellfun(@(x) fLocalNormalized(x, 'minmax'),neuralActivity_by_datapoint_diff,'UniformOutput',0);
%                 normalizedDFF_by_site_diff=cellfun(@(x) fLocalNormalized(x, 'minmax'),neuralActivity_by_site_diff,'UniformOutput',0);
%                 normalizedDFF_by_datapoint_diff=cellfun(@(x) fLocalNormalized(x, 'z-score'),neuralActivity_by_datapoint_diff,'UniformOutput',0);
%                 normalizedDFF_by_site_diff=cellfun(@(x) fLocalNormalized(x, 'z-score'),neuralActivity_by_site_diff,'UniformOutput',0);

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
                saveas(figMeanTrace,[savepath,filesep,animal_name,'-',experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-sort ',behEventSort,'-masklick-',masklick,'-n',num2str(n_datapoint(nfiber)),'-mean trace.pdf'],'pdf');
                %             saveas(fig,[savepath,filesep,animal_name,'-',experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-sort ',behEventSort,'-masklick-',masklick,'-raster.pdf'],'pdf');
            else%use common feature, e.g. experiment/group name
                if strcmp(masklick,'yes')
                    suptitle(strcat(experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-',dffstr{ndff},'-masklick'));
                else
                    suptitle(strcat(experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-',dffstr{ndff}));
                end
                if length(unique(StructDiff(iData).label))==2
                    sigLabel='AUC significance';
                elseif length(unique(StructDiff(iData).label))==3
                    sigLabel='ANOVA significance';
                end
                saveas(figMeanTrace,[savepath,filesep,experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-sort ',behEventSort,'-masklick-',masklick,'-n',num2str(n_datapoint(nfiber)),'-mean trace.pdf'],'pdf');
                saveas(figCases,[savepath,filesep,experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-masklick-',masklick,'-individual cases-',sigLabel,'.pdf'],'pdf');
                saveas(figSite,[savepath,filesep,experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-masklick-',masklick,'-individual cases by site-',sigLabel,'.pdf'],'pdf');
                saveas(figCasesDiff,[savepath,filesep,experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-masklick-',masklick,'-individual cases difference-',sigLabel,'.pdf'],'pdf');
                saveas(figSiteDiff,[savepath,filesep,experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-masklick-',masklick,'-individual cases difference by site-',sigLabel,'.pdf'],'pdf');
            end
        end
        
        n_datapointCell{i_celltype}=n_datapoint;
        n_siteCell{i_celltype}=n_site;  
    end
    [celltypeVar{1:length(StructDiff)}]=deal(celltype{i_celltype});
    [answerVar{1:length(StructDiff)}]=deal('correct');
    celltypeVar=categorical(celltypeVar);
    answerVar=categorical(answerVar);
    tempT=struct2table(StructDiff);
    tempT.celltype=celltypeVar';
    tempT.Answer=answerVar';
    if i_celltype==1 
        Tdifficulty=tempT;
    else
        Tdifficulty=vertcat(Tdifficulty,tempT);
    end
    clearvars 'celltypeVar' 'answerVar'
    %for site
    [celltypeVar{1:length(StructDiffSite)}]=deal(celltype{i_celltype});
    [answerVar{1:length(StructDiffSite)}]=deal('correct');
    celltypeVar=categorical(celltypeVar);
    answerVar=categorical(answerVar);
    tempTsite=struct2table(StructDiffSite);
    tempTsite.celltype=celltypeVar';
    tempTsite.Answer=answerVar';
    if i_celltype==1 
        TdifficultySite=tempTsite;
    else
        TdifficultySite=vertcat(TdifficultySite,tempTsite);
    end
    clearvars 'celltypeVar' 'answerVar'
end
%}

color_celltype={'F16820','646464'};%vglut2,vgat,ALM
color_celltype=fHex2RGB(color_celltype);
%     color_celltype={[1 0.5 0.5],[0.5 0.5 1],[0.5,0.5,0.5]};%blue-vgat,red-vglut2
color_celltype_mean=color_celltype;%blue-vglut2,red-vgat,black-ALM
experiment={'SC vglut2','SC vgat'};
celltypestr={'vglut2','vgat'};
titlestr={'Correct','Error','Miss','Violation'};
%plot AUC together
figDiffAUC=figure;
set(gcf,'PaperPosition',[1,1,2*length(experiment),4]);
for ncelltype=1:length(experiment)
    dataInfoStr=strcat(experiment{ncelltype},'-',regionstr{i_region},'-',fiberstr{1});
    fileNameAUC=strcat(savepath,filesep,dataInfoStr,dataProcessStr,'-n',num2str(n_datapointCell{ncelltype}(1)),'-sites',num2str(n_siteCell{ncelltype}(1)),'-difficulty.mat');
    load(fileNameAUC);
    [cellAUC{1:length(StructDiff)}]=deal(StructDiff.AUC);
    matAUC=cell2mat(cellAUC');
    [cellAUCsite{1:length(StructDiffSite)}]=deal(StructDiffSite.AUC);
    matAUCsite=cell2mat(cellAUCsite');
    figure(figDiffAUC);
    subplot(2,2,ncelltype);
    plot([ts(1),ts(end)],[0.5,0.5],'k-');hold on;
    plot(ts,matAUC,'color',color_celltype_mean{ncelltype});
    ylim([0,1]);
    xlim([ts(1),ts(end)]);
    ylabel('AUC distinguishing easy/hard');
    xlabel('Time (s) from delay onset');
    title([experiment{ncelltype},'-by sesseion']);
    subplot(2,2,ncelltype+2);
    plot([ts(1),ts(end)],[0.5,0.5],'k-');hold on;
    plot(ts,matAUCsite,'color',color_celltype_mean{ncelltype});
    ylim([0,1]);
    xlim([ts(1),ts(end)]);
    ylabel('AUC distinguishing easy/hard');
    xlabel('Time (s) from delay onset');
    title([experiment{ncelltype},'-by site']);
end
saveas(figDiffAUC,[savepath,filesep,'difficultyMovingAUC-algin to ',behEventAlign,'-masklick-',masklick,'.pdf'],'pdf');

%plot early and late AUC
Tvar=Tdifficulty.Properties.VariableNames;
tempAUCvar=cellfun(@(x) contains(x,'AUC'),Tvar);
tempANOVAvar=cellfun(@(x) contains(x,'ANOVA'),Tvar);
if sum(tempAUCvar)>0 % have AUC value
    pSig=0.05;
    figEpochAUC=figure;
    set(gcf,'PaperPosition',[1,1,4,4]);
    subplot(2,2,1);
    fScatterAUCCmpOneEpoch(Tdifficulty,pSig,'early');
    title('early AUC by session');
    subplot(2,2,2);
    fScatterAUCCmpOneEpoch(Tdifficulty,pSig,'late');
    title('late AUC by session');
    subplot(2,2,3);
    fScatterAUCCmpOneEpoch(TdifficultySite,pSig,'early');
    title('early AUC by site');
    subplot(2,2,4);
    fScatterAUCCmpOneEpoch(TdifficultySite,pSig,'late');
    title('late AUC by site');
    saveas(figEpochAUC,[savepath,filesep,'difficultyEpochAUC-algin to ',behEventAlign,'-masklick-',masklick,'.pdf'],'pdf');
elseif sum(tempANOVAvar)>0 % have ANOVA p value
    
end
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
matdata=cell(size(trialtype,1),size(trialtype,2),size(trialtype,4));
matdataraw=cell(size(trialtype,1),size(trialtype,2),size(trialtype,4));
indempty=cellfun(@isempty, dff);
dff=dff(~indempty);
trialtypecell=trialtypecell(~indempty);
celldata=cell(size(trialtype,1),size(trialtype,2),length(dff),size(trialtype,4));%1d-nResult,2d-nStim,3d-nsession,4d-opto/ctrl
for nResult=1:size(trialtype,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
    for nStim=1:size(trialtype,2) %for each stimulus%[1,6]%just 2 end trials
        for iOpto=1:size(trialtype,4)
            matdataraw{nResult,nStim,iOpto}=zeros(length(dff),size(dff{1},2));
            for idatapoint=1:length(dff)            
                selectedTrialInd=trialtypecell{idatapoint}(nResult,nStim,:,iOpto);
                selectedTrialInd=logical(squeeze(selectedTrialInd))';
                rawdata=dff{idatapoint}(selectedTrialInd,:);
                matdataraw{nResult,nStim,iOpto}(idatapoint,:)=nanmean(rawdata,1);
                celldata{nResult,nStim,idatapoint,iOpto}=rawdata;
            end
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
    maxDataVal=max(max(max(maxDataMat)));
    minDataVal=min(min(min(minDataMat)));
    
    for nResult=1:size(trialtype,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
        for nStim=1:size(trialtype,2) %for each stimulus%[1,6]%just 2 end trials
            if strcmp(method,'minmax')
                matdata{nResult,nStim,iOpto}(idatapoint,:)=(matdataraw{nResult,nStim,iOpto}(idatapoint,:)-minDataVal)/(maxDataVal-minDataVal);%method 1
            elseif strcmp(method,'z-score')
                matdata{nResult,nStim,iOpto}(idatapoint,:)=(matdataraw{nResult,nStim,iOpto}(idatapoint,:)-nanmean(nanmean(dff{idatapoint})))/nanstd(dff{idatapoint},0,'all');%method 2, z-score
            end
        end
    end
end
end

%from cell to struct for further AUC
function [Structout]= fCell2CategoryStruct(cellin)
ncategory=length(cellin);
data=[];
label=[];
for i=1:ncategory
    label=[label;ones(size(cellin{i},1),1)*i];
    data=[data;cellin{i}];
end
Structout.label=label;
Structout.data=data;
end

%from cases to get diff
function [datadiff]= fContraIpsiDiff(dataRaw,preferedStr)
datadiff=cell(size(dataRaw,1),size(dataRaw,2)/2,size(dataRaw,3));
for nResult=1:size(dataRaw,1)
    for iData=1:size(dataRaw,3)
        cell_by_datapoint=dataRaw(nResult,:,iData);
        nonprefered=cell_by_datapoint(1:length(cell_by_datapoint)/2);
        prefered=cell_by_datapoint(length(cell_by_datapoint)/2+1:end);
        prefered=fliplr(prefered);
        if strcmp(preferedStr,'prefered-non')
            datadiff(nResult,:,iData)=cellfun(@(x,y) y-nanmean(x,1),nonprefered,prefered,'UniformOutput',false);
        elseif strcmp(preferedStr,'prefered')
            datadiff(nResult,:,iData)=prefered;
        end
    end
end
end

%use to plot individual cases
function [figCases,curve_caseTrace] = fPlotCases(ts,neuralActivity,col_figcase,color_cases_trace,behEventAlign,ylabelstr,frameNumTime,selectivitystr,stimDur,SigX)
n_cases=size(neuralActivity,3);
row_figcase=ceil(n_cases/col_figcase);
figCases=figure;
set(gcf, 'paperPosition', [0 0 2*col_figcase 2*row_figcase]);
figure(figCases);
for i_case=1:n_cases
    for nResult=1:size(neuralActivity,1)
        for nStim=1:size(neuralActivity,2)
            subplot(row_figcase,col_figcase,i_case);hold on;
%             plot(ts,neuralActivity{nResult,nStim,i_case},'color',color_cases_trace{nResult,nStim});
            curve_caseTrace(nStim)=fPlotMean_SE( ts,neuralActivity{nResult,nStim,i_case},color_cases_trace{nResult,nStim});
        end
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
            ylabel(ylabelstr);
        end
        if i_case==n_cases
            switch selectivitystr
                case 'choice'
                    if length(curve_caseTrace)==2
                        legend(curve_caseTrace(:),'ctrl','opto','AutoUpdate','off');
                    elseif length(curve_caseTrace)==4
                        legend(curve_caseTrace(:),'ipsi ctrl','ipsi opto','contra opto','contra ctrl','AutoUpdate','off');
                    end
                    %                             legend(curve_caseTrace(:),'ipsi correct','contra correct','ipsi error','contra error');
                case 'difficulty'
                    if length(curve_caseTrace)==4
                        legend(curve_caseTrace(:),'ipsi easy','ipsi hard','contra hard','contra easy','AutoUpdate','off');
                    elseif length(curve_caseTrace)==2
                        legend(curve_caseTrace(:),'easy','hard','AutoUpdate','off');
                    end
                case 'stimuli'
                    if length(curve_caseTrace)==6
                        legend(curve_caseTrace(:),'','ipsi','','','contra','','AutoUpdate','off');
                    elseif length(curve_caseTrace)==3
                        legend(curve_caseTrace(:),'easy','medium','hard','AutoUpdate','off');
                    end
            end
        end
    end
    %plot bar indicating significance
    ySig=ones(size(ts))*0.9*y_lim(end);
    tempts=ts;
    tempts(~SigX{i_case})=nan;
    ySig(~SigX{i_case})=nan;
    plot(tempts,ySig,'k-','LineWidth',2);
end

end

%plot AUC individual cases
function [figCases] = fPlotAUCCases(ts,AUCarray,col_figcase,casename,SigX)
n_case=length(AUCarray);
nrow=ceil(n_case/col_figcase);
figCases=figure;
for i=1:n_case
    subplot(nrow,col_figcase,i);
    plot(ts,AUCarray{i},'k-');
    hold on;
    plot([ts(1),ts(end)],[0.5,0.5],'k-');
    plot([0,0,],[0,1],'k-');
    title(casename{i});
    ylim([0,1]);
    xlim([ts(1),ts(end)])
    if mod(i,col_figcase)==1
        ylabel('AUC');
    end
    ySig=ones(size(ts))*0.9;
    tempts=ts;
    tempts(~SigX{i})=nan;
    ySig(~SigX{i})=nan;
    plot(tempts,ySig,'k-','LineWidth',2);
    
end
end


%plot AUC for one epoch
function [Tout] = fScatterAUCCmpOneEpoch(Tin,pSig,epoch)
if strcmp(epoch,'early')
    ind_col=8:9;
    ylabelstr='stim AUC';
elseif strcmp(epoch,'late')
    ind_col=10:11;
    ylabelstr='late delay AUC';
end
TAUCVglut2=Tin(logical((Tin.celltype=='SC vglut2').*(Tin.Answer == 'correct')),ind_col);
TAUCVglut2=table2array(TAUCVglut2);
indSigVglut2=logical((TAUCVglut2(:,2)<pSig/2)+(TAUCVglut2(:,2)>1-pSig/2));
scatter(ones(length(TAUCVglut2(~indSigVglut2,1)),1),TAUCVglut2(~indSigVglut2,1),10,'k','filled')
hold on;
scatter(ones(length(TAUCVglut2(indSigVglut2,1)),1),TAUCVglut2(indSigVglut2,1),10,'r','filled')
TAUCVgat=Tin(logical((Tin.celltype=='SC vgat').*(Tin.Answer == 'correct')),ind_col);
TAUCVgat=table2array(TAUCVgat);
indSigVgat=logical((TAUCVgat(:,2)<pSig/2)+(TAUCVgat(:,2)>1-pSig/2));
scatter(ones(length(TAUCVgat(~indSigVgat,1)),1)*2,TAUCVgat(~indSigVgat,1),10,'k','filled')
scatter(ones(length(TAUCVgat(indSigVgat,1)),1)*2,TAUCVgat(indSigVgat,1),10,'r','filled')
% % TAUCALM=Tin(logical((Tin.celltype=='ALM terminal').*(Tin.Answer == 'correct')),ind_col);
% % TAUCALM=table2array(TAUCALM);
% % indSigALM=logical((TAUCALM(:,2)<pSig/2)+(TAUCALM(:,2)>1-pSig/2));
% % scatter(ones(length(TAUCALM(~indSigALM,1)),1)*3,TAUCALM(~indSigALM,1),10,'k','filled')
% % scatter(ones(length(TAUCALM(indSigALM,1)),1)*3,TAUCALM(indSigALM,1),10,'r','filled')
% [h,p1]=ttest2(TAUCVglut2(:,1),TAUCVgat(:,1));
% [h,p2]=ttest2(TAUCVglut2(:,1),TAUCALM(:,1));
% [h,p3]=ttest2(TAUCALM(:,1),TAUCVgat(:,1));
[h,varp1]=vartest2(TAUCVglut2(:,1),TAUCVgat(:,1));
% % [h,varp2]=vartest2(TAUCVglut2(:,1),TAUCALM(:,1));
% % [h,varp3]=vartest2(TAUCALM(:,1),TAUCVgat(:,1));
p1=ranksum(TAUCVglut2(:,1),TAUCVgat(:,1));
% % p2=ranksum(TAUCVglut2(:,1),TAUCALM(:,1));
% % p3=ranksum(TAUCALM(:,1),TAUCVgat(:,1));
xlim([0,3]);
ylim([0,1]);
plot([1,1.9],[0.8,0.8],'k-');
text(1,0.9,plabelsymbol(p1));
% % plot([2.1,3],[0.8,0.8],'k-');
% % text(2.5,0.9,plabelsymbol(p2));
% % plot([1,3],[0.9,0.9],'k-');
% % text(1.5,1,plabelsymbol(p3));
indSigContraVglut2=logical(TAUCVglut2(:,2)<pSig/2);
indSigContraVgat=logical(TAUCVgat(:,2)<pSig/2);
% % indSigContraALM=logical(TAUCALM(:,2)<pSig/2);
text(0.5,0.1,strcat(num2str(sum(indSigContraVglut2)),'/',num2str(length(indSigContraVglut2)),'hard'));
text(1.5,0.1,strcat(num2str(sum(indSigContraVgat)),'/',num2str(length(indSigContraVgat)),'hard'));
% % text(2.5,0.1,strcat(num2str(sum(indSigContraALM)),'/',num2str(length(indSigContraALM)),'contra'));
set(gca,'XTick',[1,2],'XTickLabel',{'vglut2','vgat'});
ylabel(ylabelstr);
plot([0,3],[0.5,0.5],'k--');
Tout=Tin;
set(gca,'FontName','Arial','FontSize',12);
end