%This code choose some sessions and group(but not concanetate) data together
% to plot PSTH and Ending point trajectories; both show individual sessions
% and animals
% dbstop if error
clear;
dbstop if error;
close all;
savepath='F:\FP\example';
[num,txt,raw] =xlsread('F:\FP\fiber photometry data summary.xlsx');%criteria to choose sessions come from this file
T=cell2table(raw(2:end,1:11));
T.Properties.VariableNames=strrep(raw(1,1:11),' ','_');%table variable name can't have ' ',so replace them
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
selectivitystr={'stimuli','sensory difficulty','sensory','choice'};%sensory means grouping difficulties;
trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
combineCorErr='divideCorErr';%{'combineCorErr','divideCorErr'}
if strcmp(combineCorErr,'combineCorErr') %used for deciding matrix size
    n_answer=1;
elseif strcmp(combineCorErr,'divideCorErr')
    n_answer=2;
end
fiberstr={'Soma'};%fiberstr={'Soma','Terminal'};
dffstr={'dff baseline correction first','dff motion correction first','dff470','dff410'};
behEventAlign='first lick';%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
masklick='no';%'yes' if mask lickings or 'no' when want to see all activity
stimDur=0.5;%0.5s click duration
if strcmp(behEventAlign,'delay onset')
    frameNumTime=[2.5,3.5];%from 2s before align point to 5s after align point
else
    frameNumTime=[0.5,1];%from 2s before align point to 5s after align point
end
yrange=[-0.03,0.15];
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
nshuffle=1000;
p_sepOnset=0.05;%the p-value used to determine when separation happen
%set the criteria to choose sessions
i_region=3;%1-bilateral, 2-unilateral*********variable**************
region={'bilateral','t SC','SC'};%{'bilateral SC','SC'}
regionstr={'bilateral','unilateral','mixed'};
whichAnimal='pyx241';%************variable******pyx237 for vgat, pyx241 for vglut2
whichSite='left';%if bilateral fiber implanted, can choose site here******variable******right for pyx237, left for pyx241
celltype={'SC vglut2','SC vgat'};
n_siteCell=cell(1,2);%each cell store a cell type
n_datapointCell=cell(1,2);%each cell store a cell type
time_before_delay=0.5;
time_before_go=0.5;
dataProcessStr=strcat('-grouped by ',selectivitystr{i_selectivity},'-',combineCorErr,'-binCenteredSize',num2str(binsize),'-alignTo-',behEventAlign,'-timeWindow',num2str(frameNumTime(1)),'-',num2str(frameNumTime(2)),'-EpochWindow stim-',num2str(time_before_delay),'delay-',num2str(time_before_go),'-shuffle',num2str(nshuffle));

%
i_celltype=1;%1 for vglut2, 2 for vgat*********variable**************
%select sessions for summary analysis
experiment=celltype{i_celltype};
if strcmp(whichAnimal,'pyx241') && strcmp(whichSite,'left') && i_celltype==1
    ind_session=strcmp(T.animal,whichAnimal).*strcmp(T.used_as_data,'yes').*(strcmp(T.date,'2019/11/30')).*strcmp(T.experiment,experiment);%for vglut2 example, pyx241
elseif strcmp(whichAnimal,'pyx237') && strcmp(whichSite,'right') && i_celltype==2
    ind_session=strcmp(T.animal,whichAnimal).*strcmp(T.used_as_data,'yes').*(strcmp(T.date,'2019/12/1')).*strcmp(T.experiment,experiment);%for vgat example, F:\FP\pyx237_20191201
end
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
    dataInfoStr=strcat(experiment,'-',whichAnimal,'-',regionstr{i_region},'-',fiberstr{nfiber});
    fileNameAUC=strcat(savepath,filesep,dataInfoStr,dataProcessStr,'-n',num2str(n_datapoint(nfiber)),'-sites',num2str(n_site(nfiber)),'-AUC.mat');
    [AUCcell,pAUCcell]=deal(cell(n_datapoint(nfiber),4));%each cell 1d-data points, 2d-one result(cor/err/miss/vio)
    datapointNameCell=cell(n_datapoint(nfiber),1);
    
    [AUCearly,pAUCearly,AUClate,pAUClate]=deal(nan(n_answer*n_datapoint(nfiber),1));%variables of same length, later will combine as tabel
    [DatapointName, AnswerSession]=deal(cell(n_datapoint(nfiber)*n_answer,1));
    [AUCcellSite,pAUCcellSite]=deal(cell(n_site(nfiber),4));%each cell 1d-site, 2d-one result(cor/err/miss/vio)
    siteNameCell=[];%later will be struct with field 'nameSite','nameCases'
    [AUCearlySite,pAUCearlySite,AUClateSite,pAUClateSite]=deal(nan(n_answer*n_site(nfiber),1));
    [SiteName,AnswerSite]=deal(cell(n_site(nfiber)*n_answer,1));
    %store ttest p-value
    fileNameTtest=strcat(savepath,filesep,dataInfoStr,dataProcessStr,'-n',num2str(n_datapoint(nfiber)),'-sites',num2str(n_site(nfiber)),'-Ttest.mat');
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
            if ~isfield(siteNameCell,'nameCases')
                siteNameCell(ind_site).nameCases=cell(1,1);
                siteNameCell(ind_site).nameCases{1}=datapointNameCell{ind_datapoint,1};
            elseif isempty(siteNameCell(ind_site).nameCases)
                siteNameCell(ind_site).nameCases{1}=datapointNameCell{ind_datapoint,1};
            else
                temp=cellfun(@(x) strcmp(x,datapointNameCell{ind_datapoint,1}),siteNameCell(ind_site).nameCases);
                if sum(temp)>0 %this session name exist
                    ind_session4site=find(temp);
                else
                    ind_session4site=length(siteNameCell(ind_site).nameCases)+1;
                end
                siteNameCell(ind_site).nameCases{ind_session4site(1)}=datapointNameCell{ind_datapoint,1};
            end
            %get trial type
            RTVector=double(Data_extract.Answer_time-Data_extract.Go_time);
            indRT=true(size(RTVector));%logical((RTVector>700).*(RTVector<900));
            [trialType,behrule,trialTypeStr] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix',fiberSide{i_fiber},combineCorErr);%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
            %                 %this block choose only easy trials
            %                 tempTrialType=zeros(size(trialType,1),2,size(trialType,3));
            %                 tempTrialType(:,1,:)=trialType(:,1,:);
            %                 tempTrialType(:,2,:)=trialType(:,end,:);
            %                 trialType=tempTrialType;
            %                 trialTypeStr=strcat(trialTypeStr,'-easy tirals');
            disp([trialTypeStr,',rule of this session is ',behrule]);
            %}
            trialType_cat_animal{nfiber,ind_site}=cat(3,trialType_cat_animal{nfiber,ind_site},trialType);
            trialType_by_session{nfiber,ind_datapoint}=trialType;
            %% concatenate dff, trial type maxtrix and behavior event, etc.
            for ndff=[1] %plot for each dff, see which is better and test whether 410 signals have similar trend
                %ITI baseline correction
                dff{1,ndff}(ind_fiber(i_fiber),:)=fDffITIcorrection( dff{1,ndff}(ind_fiber(i_fiber),:), behEventFrameIndex,FrameRate/2);
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
            %save corresponding datapoint name- which session and which sites
            %
            fileNameAUC_datapoint=strcat(savepath,filesep,dataInfoStr,'-',datapointNameCell{ind_datapoint,1},dataProcessStr,'-AUC.mat');%adding datapoint name to the AUC file name
            if n_answer==2
                stranswer={'correct','error'};
            elseif n_answer==1
                stranswer={'correct and error'};
            end
            if exist(fileNameAUC_datapoint,'file')
                load(fileNameAUC_datapoint);
                for nResult=1:size(trialType,1)-2
                    AUCcell{ind_datapoint,nResult} = AUCcase(nResult).AUC;
                    pAUCcell{ind_datapoint,nResult} = AUCcase(nResult).pAUC;
                end
                if ~strcmp(datapointNameCell{ind_datapoint,1},nameCase)
                    warning('Name of data point mismatch, double check');
                end
            else
                nameCase=datapointNameCell{ind_datapoint,1};
            end
            for nResult=1:size(trialType,1)
                %calculate moving p to decide when selectivity becoming significant
                label = fTrialType2Label(trialType,2);
                indTrial=trialType(nResult,:,:);
                indTrial=sum(squeeze(indTrial),1);
                indTrial=logical(squeeze(indTrial));
                %plot moving AUC, method 1
                if nResult<size(trialType,1)-1 && nfiber==1 %only calculate AUC from cor/err soma
                    if true%isempty(AUCcell{ind_datapoint,nResult})
                        [AUCcell{ind_datapoint,nResult},pAUCcell{ind_datapoint,nResult}] = fMovingAUC(label(indTrial),dff_aligned(indTrial,:),2,nshuffle,binsize,binstep);
                        AUCcase(nResult).AUC=AUCcell{ind_datapoint,nResult};
                        AUCcase(nResult).pAUC=pAUCcell{ind_datapoint,nResult};
                    end

                end
                %calculate moving t-test, method 2
                pTtestcell{ind_datapoint,nResult}=fMovingTtest(label(indTrial),dff_aligned(indTrial,:),binsize,binstep);
            end
            save(fileNameAUC_datapoint,'AUCcase','nameCase');
            %}
            ind_datapoint=ind_datapoint+1;
            
        end
    end
    %calculate AUC for large sessions grouped by sites, only calculated soma, not terminal
    %
    for ind_site=1:n_site(1)
        fileNameAUC_site=strcat(savepath,filesep,dataInfoStr,'-',siteNameCell(ind_site).nameSite,dataProcessStr,'-AUC.mat');%adding datapoint name to the AUC file name
        fileNameAUC_site=fileNameAUC_site{1};
        if n_answer==2
            stranswer={'correct','error'};
        elseif n_answer==1
            stranswer={'correct and error'};
        end
        if exist(fileNameAUC_site,'file')
            load(fileNameAUC_site);
            %check whether site include same dataset
            if fEqual(nameCases,siteNameCell(ind_site).nameCases)
                for nResult=1:size(trialType,1)-2
                    AUCcellSite{ind_site,nResult} = AUCsite(nResult).AUC;
                    pAUCcellSite{ind_site,nResult} = AUCsite(nResult).pAUC;
                end
            else
                nameSite=siteNameCell(ind_site).nameSite;
                nameCases=siteNameCell(ind_site).nameCases;
            end
        else
            nameSite=siteNameCell(ind_site).nameSite;
            nameCases=siteNameCell(ind_site).nameCases;
        end
        nameCases=unique(nameCases);%remove elements that replicate
        siteNameCell(ind_site).nameCases=unique(siteNameCell(ind_site).nameCases);
        for nResult=1:1:size(trialType,1)-2 %only calculated cor/err trials
            label = fTrialType2Label(trialType_cat_animal{1,ind_site},2);
            indTrial=trialType_cat_animal{1,ind_site}(nResult,:,:);
            indTrial=sum(squeeze(indTrial),1);
            indTrial=logical(squeeze(indTrial));
            if true%isempty(AUCcellSite{ind_site,nResult}) && nResult<=2
                [AUCcellSite{ind_site,nResult},pAUCcellSite{ind_site,nResult}] = fMovingAUC(label(indTrial),dff_aligned_cat_animal{1,ind_site}(indTrial,:),2,nshuffle,binsize,binstep);
                AUCsite(nResult).AUC=AUCcellSite{ind_site,nResult};
                AUCsite(nResult).pAUC=pAUCcellSite{ind_site,nResult};
            end
            %calculate moving t-test, method 2
            pTtestcellSite{ind_site,nResult}=fMovingTtest(label(indTrial),dff_aligned_cat_animal{1,ind_site}(indTrial,:),binsize,binstep);
        end
        save(fileNameAUC_site,'AUCsite','nameSite','nameCases');
    end
    %}
    %plot PSTH
    indempty=cellfun(@isempty,datapointNameCell);
    indCase=cellfun(@(x) contains(x,whichSite),datapointNameCell(~indempty));
    [~,neuralActivity_by_site]=fDFFbyTrialType(dff_aligned_cat_animal(1,indCase),trialType_cat_animal,indRT);
    [~,neuralActivity_by_datapoint]=fDFFbyTrialType(dff_aligned_by_session(1,indCase),trialType_by_session,indRT);
    neuralActivity2Plot=cell(1,2);
    neuralActivity2Plot{1}=neuralActivity_by_site;
    neuralActivity2Plot{2}=neuralActivity_by_datapoint;
    %% plot raster and mean trace for large session, and individual cases by trial
    for ndff=[1]         %  titlestr=strcat(animal_name,'-',fiberstr{nfiber},'-',titlestr);
        if isempty(neuralActivity2Plot)%if no data, then end this loop
            break;
        end
        figraster=figure;%plot raster
        set(gcf, 'position', [0 400 1400 300]);
        figMeanTrace=figure;%plot mean trace
        set(gcf, 'PaperPosition', [0 0 sum(frameNumTime)*2 4]);

        if contains(trialTypeStr,'combineCorErr')
            titlestr={'do','Miss','Violation'};
        else
            titlestr={'Correct','Error','Miss','Violation'};
        end
        titlestr_row={'-by site','-by session'};
        choicestr={'ipsi','contra'};
        ts=double((-frameNum(1):binstep:frameNum(2))*frT/1000);
        %             ColLimit = prctile(dff_aligned_cat{nfiber}',98,'all');%here data is a matrix rather than a vector, so use 'all' to find the prctile of all data（ver2018b)
        for figrow=1:2
            n_col=size(neuralActivity2Plot{figrow},1)-2;%4 column(correct/error/miss/violation),companied with 4 lick raster
            for nResult=1: n_col
                %% plot raster plot and PSTH
                for nStim=1:size(neuralActivity2Plot{figrow},2) %for each stimulus%[1,6]%just 2 end trials
                    if nResult==(size(neuralActivity2Plot{figrow},1)-1) && size(neuralActivity2Plot{figrow},2)==2%miss trials && grouped with choice rather than stimuli
                        color_mean_trace={[0,0,0],[0,0,0]};
                    elseif size(neuralActivity2Plot{figrow},2)==6%using stimuli to group trials
                        color_mean_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};%left-low--right-high;since cat sessions together, color don't care rule, just using the standard condition
                    elseif size(neuralActivity2Plot{figrow},2)==4%using difficulty to group trials
                        color_mean_trace={[0 0 1],[0.3 0 0.7],[0.7 0 0.3],[1 0 0]};%left-low--right-high;since cat sessions together, color don't care rule, just using the standard condition
                    elseif contains(trialTypeStr,'first lick')
                        color_mean_trace={[0 0 1],[1 0 0]};
                    else
                        color_mean_trace=fMeanTraceColor(behrule,size(trialType,2));
                    end
                    neuralActivity=neuralActivity2Plot{figrow}{nResult,nStim};
                    ColLimit=[0,0.1];
                    %raster plot of activity and corresponding licking
                    figure(figraster);
                    axActivity=subplot(size(trialType_cat_animal,2),n_col,nResult);
                    x_lim=[0,sum(frameNum)];
                    fRasterDffandLicking(x_lim,axActivity,neuralActivity,behEvent_aligned_cat_animal{1},[],behEventSort,ColLimit);
%                     %add labels to activity raster block
%                     subplot(axActivity);
%                     if nStim == 1
%                         title(titlestr{nResult});
%                     end
%                     if nResult==1&&contains(trialTypeStr,'stimuli')
%                         ylabel([num2str(Data_extract.Stimuli(nStim))]);% ' clicks/s'
%                     elseif nResult==1&&contains(trialTypeStr,'first lick')
%                         ylabelstr={'ipsi lick first','contra lick first'};
%                         ylabel(ylabelstr{nStim});
%                     elseif nResult==1&&contains(trialTypeStr,'sensory')
%                         ylabelstr={'low','high'};
%                         ylabel(ylabelstr{nStim});
%                     end
%                     xlabel(['time(s) from ',behEventAlign]);
                    %plot piled activities for each trials
                    if strcmp(whichAnimal,'pyx237')% for vgat case
                        indTrial=[];
%                         if nStim==2 %contra trial
%                             indTrial=[9,10,12];%[];
%                         else
%                             indTrial=[5,6,11];
%                         end
                    elseif strcmp(whichAnimal,'pyx241')
                        indTrial=[];
%                         if nStim==2 %contra trial
%                             indTrial=[12,13,15];%[];
%                         else
%                             indTrial=[4,10,12];
%                         end
                    end
                    xlabelstrPiled=['Time (s) from ',behEventAlign];
                    if figrow==1 && strcmp(titlestr{nResult},'Correct')%each row same data, since here only one example session
                        [figPiled] = fPlotPiledDff(neuralActivity,ts,indTrial,color_mean_trace{nStim});
                        figure(figPiled);
                        xlabel(xlabelstrPiled);
%                         title(titlestr{nResult});
                        set(figPiled,'PaperPosition',[1,1,1,1]);
                        box off;
                        saveas(figPiled,[savepath,filesep,datapointNameCell{indCase},'-',titlestr{nResult},'-',choicestr{nStim},'-single trial dff-align to-',behEventAlign,'.pdf'],'pdf');
                    end
                    %plot mean trace
                    figure(figMeanTrace);%save mean trace
                    subplot(2,n_col,nResult+(figrow-1)*n_col);
                    
%                     curve_meanTrace(nStim)=fPSTHBinDelayEndPoint(neuralActivity,color_mean_trace{nStim},'EndSmooth','smooth','EndSmoothMethod','bin','ParaEndSmoothMethod',0.3*1000/frT);
                    curve_meanTrace(nStim)=fPlotMean_SE(ts,neuralActivity,color_mean_trace{nStim});
                    if ~exist('yrange','var')
                        if ~exist('y_lim','var')
                            y_lim=get(gca,'Ylim');
                        else
                            y_lim_temp=get(gca,'Ylim');
                            y_lim=[min(y_lim(1), y_lim_temp(1)),max(y_lim(end), y_lim_temp(end))];
                        end
                    else
                        y_lim=yrange;
                    end
                    %                     set(gca,'Ylim',yrange);%for comparing different cell type

                end
                
                plot([0,0],[y_lim(1),y_lim(2)],'k-');
                if strcmp(behEventAlign,'stim onset')
                    plot([0.5,0.5],[y_lim(1),y_lim(2)],'k-');
                elseif strcmp(behEventAlign,'delay onset')
                    plot([-0.5,-0.5],[y_lim(1),y_lim(2)],'k');
                end
                xlabel(['Time (s) from ',behEventAlign]);
                title(strcat(titlestr{nResult},titlestr_row{figrow}));
                set(gca,'xtick',[-floor(frameNumTime(1)):0.5:frameNumTime(2)],'xticklabel',[-floor(frameNumTime(1)):0.5:frameNumTime(2)]);

                xlim([-frameNumTime(1),frameNumTime(2)]);
                ylim(y_lim);
                set(gca,'FontName','Arial','FontSize',10);
                box off;
                if strcmp(behEventAlign,'first lick')
                    ax=gca;
                    ax.YAxisLocation='right';
                end
%                 if nResult==1
%                     ylabel('\it\DeltaF/F');
%                 end
            end
        end
        figure(figMeanTrace);
        %text label
        subplot(2,size(neuralActivity2Plot{figrow},1),3);%miss plot has most space
        if contains(trialTypeStr,'stimuli')
            h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','','contra hard','','contra easy'},'Location','best');
        elseif contains(trialTypeStr,'difficulty')
            h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','contra hard','contra easy'},'Location','best');
        elseif contains(trialTypeStr,'first lick')
            h=legend(curve_meanTrace(:),{'ipsi lick first','contra lick first'},'Location','best');
        elseif contains(trialTypeStr,'sensory')
            h=legend(curve_meanTrace(:),{'low click rate','high click rate'},'Location','best');
        end
        set(h,'box','off');
%         text(x_lim(1),y_lim(end),strcat('n=',num2str(n_datapoint),' datapoints',num2str(n_session),' sessions,',num2str(n_site(nfiber)),' site,',num2str(n_animal),' animal'),'FontName','Arial','FontSize',14);
        %save figure
        if n_animal==1%only one animal
            if strcmp(masklick,'yes')
                suptitle(strcat(T.animal{ind_session(1)},'-',regionstr{i_region},'-',fiberstr{nfiber},'-',dffstr{ndff},'-masklick'));
            else
                suptitle(strcat(T.animal{ind_session(1)},'-',regionstr{i_region},'-',fiberstr{nfiber},'-',dffstr{ndff}));
            end
            saveas(figMeanTrace,[savepath,filesep,T.animal{ind_session(1)},'-',experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-sort ',behEventSort,'-masklick-',masklick,'-n',num2str(n_datapoint(nfiber)),'-mean trace.pdf'],'pdf');
            %             saveas(fig,[savepath,filesep,animal_name,'-',experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-sort ',behEventSort,'-masklick-',masklick,'-raster.pdf'],'pdf');
        else%use common feature, e.g. experiment/group name
            if strcmp(masklick,'yes')
                suptitle(strcat(experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-',dffstr{ndff},'-masklick'));
            else
                suptitle(strcat(experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-',dffstr{ndff}));
            end
            saveas(figMeanTrace,[savepath,filesep,experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-sort ',behEventSort,'-masklick-',masklick,'-n',num2str(n_datapoint(nfiber)),'-mean trace.pdf'],'pdf');
            %             saveas(fig,[savepath,filesep,experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-sort ',behEventSort,'-masklick-',masklick,'-raster.pdf'],'pdf');
        end
        %             saveas(figMovingAUC,[savepath,filesep,experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-n',num2str(n_datapoint(nfiber)),'-AUC.pdf'],'pdf');
        %             close all;
    end
    n_datapointCell{i_celltype}=n_datapoint;
    n_siteCell{i_celltype}=n_site;
    %save AUC, because calculation is time consuming
%     if strcmp(behEventAlign,'delay onset')
%         TAUCepochSession=table(string(DatapointName),AUCearly,pAUCearly, AUClate,pAUClate,categorical(AnswerSession),...
%             'VariableNames',{'DatapointName' 'AUCearly' 'pAUCearly' 'AUClate' 'pAUClate' 'Answer'});
%         TAUCepochSite=table(string(SiteName),AUCearlySite,pAUCearlySite,AUClateSite,pAUClateSite,categorical(AnswerSite),...
%             'VariableNames',{'SiteName' 'AUCearly' 'pAUCearly' 'AUClate' 'pAUClate' 'Answer'});
%         save(fileNameAUC,'AUCcell','pAUCcell','datapointNameCell','AUCcellSite','pAUCcellSite','siteNameCell','TAUCepochSite','TAUCepochSession');
%     else
        save(fileNameAUC,'AUCcell','pAUCcell','datapointNameCell','AUCcellSite','pAUCcellSite','siteNameCell');
%     end
    %         save(fileNameAUC,'AUCcell','pAUCcell','datapointNameCell','AUCcellSite','pAUCcellSite','siteNameCell');
    save(fileNameTtest,'pTtestcell','pTtestcellSite');
end

%}
%
%% compare different celltype p value(AUC, ttest)
%settings for which case to plot
ncelltype=i_celltype;%1-for vglut cell,while 2- vgat
iRes=1;%correct-iRes=1; error- iRes=2

if i_selectivity>=2%only when 2 group exist
    color_celltype={'F16820','646464'};%
    color_celltype=fHex2RGB(color_celltype);
    color_celltype_mean=color_celltype;%
    celltypestr={'vglut2','vgat'};
    ylabelstr={'AUC','p of AUC','p of t-test'};
    titlestr={'Correct','Error','Miss','Violation'};
    ytext=[1,0.1];
    
    
    AUCbySession{ncelltype}=cellfun(@(x) reshape(x(:),1,[]),AUCcell(:,iRes),'UniformOutput' ,false);
    AUCbySession{ncelltype}=cell2mat(AUCbySession{ncelltype}(indCase));
    AUCbySite{ncelltype}=cellfun(@(x) reshape(x(:),1,[]),AUCcellSite(:,iRes),'UniformOutput' ,false);
    AUCbySite{ncelltype}=cell2mat(AUCbySite{ncelltype}(indCase));
%     Ttemp = fExtendT(TAUCepochSession,celltypestr{ncelltype},'celltype');
%     if exist('TAUCepochSessionCombine','var')
%         TAUCepochSessionCombine=vertcat(TAUCepochSessionCombine,Ttemp);
%     else
%         TAUCepochSessionCombine=Ttemp;
%     end
%     Ttemp = fExtendT(TAUCepochSite,celltypestr{ncelltype},'celltype');
%     if exist('TAUCepochSiteCombine','var')
%         TAUCepochSiteCombine=vertcat(TAUCepochSiteCombine,Ttemp);
%     else
%         TAUCepochSiteCombine=Ttemp;
%     end
%     %refine table
%     TAUCepochSessionCombine=TAUCepochSessionCombine(contains(TAUCepochSessionCombine.DatapointName,whichSite),:);
%     TAUCepochSiteCombine=TAUCepochSiteCombine(contains(TAUCepochSiteCombine.SiteName,whichSite),:);

    %plot AUC separately for vglut2/vgat
    figAUC=figure;
    set(gcf,'paperPosition',[0 0 sum(frameNumTime)*2 4]);
    for icol=ncelltype%each col one cell type
        figure(figAUC);
        subplot(2,2,icol);%1st row
%         AUCbySession{icol}=abs(AUCbySession{icol}-0.5)+0.5;
%         plot(ts,AUCbySession{icol},'Color',color_celltype{icol},'LineWidth',0.5); hold on;
        curve(icol)=fPlotMean_CI(ts,AUCbySession{icol},0.05,color_celltype_mean{icol});
        titlestrAUC=strcat(celltypestr{icol},'-by-session');
        title(titlestrAUC);
        subplot(2,2,icol+2);%2nd row
%         AUCbySite{icol}=abs(AUCbySite{icol}-0.5)+0.5;
%         plot(ts,AUCbySite{icol},'Color',color_celltype{icol},'LineWidth',0.5);   hold on;
        curve(icol)=fPlotMean_CI(ts,AUCbySite{icol},0.05,color_celltype_mean{icol});
        titlestrAUC=strcat(celltypestr{icol},'-by-site');
        title(titlestrAUC);
        for irow=1:2
            subplot(2,2,icol+irow*2-2);%2nd row
            plot(ts,ones(length(ts),1)*0.5,'Color',[0 0 0]);%show 0.5
            plot([0,0],[0,1],'k-');
            if strcmp(behEventAlign,'stim onset')
                plot([0.5,0.5],[0,1],'k-');
            elseif strcmp(behEventAlign,'delay onset')
                plot([-0.5,-0.5],[0,1],'k');
            end
%             ylabel('AUC');
            xlabel(['Time (s) from ',behEventAlign]);
            xlim([-frameNumTime(1),frameNumTime(2)]);
            set(gca,'xtick',[-floor(frameNumTime(1)):0.5:frameNumTime(2)],'xticklabel',[-floor(frameNumTime(1)):0.5:frameNumTime(2)]);
            set(gca,'FontName','Arial','FontSize',10);
            box off;
            if strcmp(behEventAlign,'first lick')
                ax=gca;
                ax.YAxisLocation='right';
            end
        end
    end
    saveas(figAUC,[savepath,filesep,datapointNameCell{indCase},'AUC-algin to-',behEventAlign,'.pdf'],'pdf');
    
%     %scatter plot of epoch AUC
%     if strcmp(behEventAlign,'delay onset')
%         pSig=0.05;
%         figAUCepoch=figure;
%         set(gcf,'paperPosition',[0,0,8,4]);
%         subplot(2,4,1);%1d row by session,1st col late AUC comparison
%         [Tout] = fScatterAUCearlyCmp(TAUCepochSessionCombine,pSig);
%         title('by session');
%         subplot(2,4,2);%1d row by session,1st col late AUC comparison
%         [Tout] = fScatterAUClateCmp(TAUCepochSessionCombine,pSig);
%         title('by session');
%         subplot(2,4,5);
%         fScatterAUCearlyCmp(TAUCepochSiteCombine,pSig);
%         title('by site');
%         subplot(2,4,6);
%         fScatterAUClateCmp(TAUCepochSiteCombine,pSig);
%         title('by site');
%         subplot(2,4,3);
%         %     fScatterCmpEarlyLate(TAUCepochSessionCombine,'vglut2',pSig);
%         fScatterLineCmpEarlyLate(TAUCepochSessionCombine,'vglut2',pSig);
%         title('vglut2');
%         subplot(2,4,4);
%         %     fScatterCmpEarlyLate(TAUCepochSessionCombine,'vgat',pSig);
%         fScatterLineCmpEarlyLate(TAUCepochSessionCombine,'vgat',pSig);
%         title('vgat');
%         subplot(2,4,7);
%         %     fScatterCmpEarlyLate(TAUCepochSiteCombine,'vglut2',pSig);
%         fScatterLineCmpEarlyLate(TAUCepochSiteCombine,'vglut2',pSig);
%         title('vglut2');
%         subplot(2,4,8);
%         %     fScatterCmpEarlyLate(TAUCepochSiteCombine,'vgat',pSig);
%         fScatterLineCmpEarlyLate(TAUCepochSiteCombine,'vgat',pSig);
%         title('vgat');
%         saveas(figAUCepoch,[savepath,filesep,regionstr{i_region},datapointNameCell{indCase},'-epoch AUC comparison.pdf'],'pdf');
%     end
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
%to plot mean+-sem
function []=fPlotMeanSem(data,x_lim,color_data)
%data should be a vector
%x_lim indicate the ceter of x limit
%color_data indicated the color of lines
mean_data=mean(data);
sem=std(data)/sqrt(length(data));
plot([x_lim-0.2,x_lim+0.2],[mean_data, mean_data],'Color',color_data);hold on;
plot(x_lim,mean_data,'o','Color',color_data);
plot([x_lim-0.1,x_lim+0.1],[mean_data-sem, mean_data-sem],'Color',color_data);
plot([x_lim-0.1,x_lim+0.1],[mean_data+sem, mean_data+sem],'Color',color_data);
plot([x_lim,x_lim],[mean_data-sem,mean_data+sem],'Color',color_data);
end
%to use trialtype and dff_aligned to get normalized plotting data
function [matdata,matdataraw]=fDFFbyTrialType(dff,trialtypecell,varargin)
if isempty(varargin)
    trialtype=trialtypecell{1};
    matdata=cell(size(trialtype,1),size(trialtype,2));
    matdataraw=cell(size(trialtype,1),size(trialtype,2));
    for nResult=1:size(trialtype,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
        for nStim=1:size(trialtype,2) %for each stimulus%[1,6]%just 2 end trials
            matdata{nResult,nStim}=zeros(length(dff),size(dff{1},2));
            for idatapoint=1:length(dff)
                selectedTrialInd=trialtypecell{1,idatapoint}(nResult,nStim,:);
                selectedTrialInd=logical(squeeze(selectedTrialInd))';
                rawdata=dff{1,idatapoint}(selectedTrialInd,:);
                matdata{nResult,nStim}(idatapoint,:)=nanmean(rawdata,1);
                matdataraw{nResult,nStim}=rawdata;
            end
        end
    end
else
    indExtra=varargin{1};
    trialtype=trialtypecell{1};
    matdata=cell(size(trialtype,1),size(trialtype,2));
    matdataraw=cell(size(trialtype,1),size(trialtype,2));
    for nResult=1:size(trialtype,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
        for nStim=1:size(trialtype,2) %for each stimulus%[1,6]%just 2 end trials
            matdata{nResult,nStim}=zeros(length(dff),size(dff{1},2));
            for idatapoint=1:length(dff)
                selectedTrialInd=trialtypecell{1,idatapoint}(nResult,nStim,:);
                selectedTrialInd=logical(squeeze(selectedTrialInd))';
                selectedTrialInd=logical(selectedTrialInd.*indExtra);
                rawdata=dff{1,idatapoint}(selectedTrialInd,:);
                matdata{nResult,nStim}(idatapoint,:)=nanmean(rawdata,1);
                matdataraw{nResult,nStim}=rawdata;
            end
        end
    end
end
%normalize for each datapoint
maxDataCell=cellfun(@(x) max(x,[],2), matdata,'UniformOutput',false);
minDataCell=cellfun(@(x) min(x,[],2), matdata,'UniformOutput',false);
for idatapoint=1:length(dff)
    maxDataMat=cellfun(@(x) x(idatapoint),maxDataCell);
    minDataMat=cellfun(@(x) x(idatapoint),minDataCell);
    maxDataVal=max(max(maxDataMat));
    minDataVal=min(min(minDataMat));
    for nResult=1:size(trialtype,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
        for nStim=1:size(trialtype,2) %for each stimulus%[1,6]%just 2 end trials
            matdata{nResult,nStim}(idatapoint,:)=(matdata{nResult,nStim}(idatapoint,:)-minDataVal)/(maxDataVal-minDataVal);
        end
    end
end
end
%to plot mean trace and CI as patch
function [outputcurve]=fPlotMean_CI(ts,neuralActivity,PCI,color)
[neuralActivityMean, neuralActivityCI]=fMean_CI(neuralActivity,PCI);
% %plot individual traces
% color_case=(1+color)/2;
% plot(1:size(neuralActivity,2),neuralActivity,'Color',color_case,'linewidth',1);
% hold on;

% %shadow as ci
% xpatch=[ts(1:size(neuralActivity,2)), fliplr(ts(1:size(neuralActivity,2)))];
% ypatch=[neuralActivityCI(1,:),fliplr(neuralActivityCI(2,:))];
% p=patch(xpatch,ypatch,color);%plot confidence interval
% p.FaceAlpha=0.05;
% p.EdgeColor='none';%color;%'none';
% hold on;

outputcurve=plot(ts(1:size(neuralActivity,2)),neuralActivityMean,'Color',color,'linewidth',2);
hold on;
end
%to plot mean and sem
function [outputcurve]= fPlotMean_SE(ts,neuralActivity,color)
[ neuralActivityMean, neuralActivitySE ] = fMean_SE( neuralActivity );
% %plot individual traces
% color_case=(1+color)/2;
% plot(1:size(neuralActivity,2),neuralActivity,'Color',color_case,'linewidth',1);
% hold on;

%shadow as se
xpatch=[ts(1:size(neuralActivity,2)), fliplr(ts(1:size(neuralActivity,2)))];
ypatch=[neuralActivitySE(1,:),fliplr(neuralActivitySE(2,:))];
p=patch(xpatch,ypatch,color);%plot confidence interval
p.FaceAlpha=0.2;
p.EdgeColor='none';%color;%'none';
hold on;

outputcurve=plot(ts(1:size(neuralActivity,2)),neuralActivityMean,'Color',color,'linewidth',2);
hold on;
end
%to compare cell array whether they are same in value but can be in
%different order
function [flag] = fEqual(A,B)
%A,B should be cell array
if length(A)~=length(B)
    flag=false(1);
    return;
end
for i=1:length(A)
    temp=cellfun(@(x) strcmp(x,A{i}),B);
    if sum(temp)==0
        flag=false(1);
        return;
    end
end
flag=true(1);
end
%add a col of celltype
function [Tout]= fExtendT(Tin,element,nameVar)
nrowT=size(Tin,1);
[addvar{1:nrowT}]=deal(element);
addvar=reshape(addvar,[],1);
addvar=categorical(addvar);
Tout=table(addvar,'VariableNames',{nameVar});
Tout=[Tout Tin];
end
%table input, plot early AUC comparison
function [Tout] = fScatterAUCearlyCmp(Tin,pSig)
TAUCearlyVglut2=Tin(logical((Tin.celltype=='vglut2').*(Tin.Answer == 'correct')),3:4);
indSigVglut2=logical((TAUCearlyVglut2.pAUCearly<pSig/2)+(TAUCearlyVglut2.pAUCearly>1-pSig/2));
scatter(ones(length(TAUCearlyVglut2.AUCearly(~indSigVglut2)),1),TAUCearlyVglut2.AUCearly(~indSigVglut2),10,'k','filled')
hold on;
scatter(ones(length(TAUCearlyVglut2.AUCearly(indSigVglut2)),1),TAUCearlyVglut2.AUCearly(indSigVglut2),10,'r','filled')
TAUCearlyVgat=Tin(logical((Tin.celltype=='vgat').*(Tin.Answer == 'correct')),3:4);
indSigVgat=logical((TAUCearlyVgat.pAUCearly<pSig/2)+(TAUCearlyVgat.pAUCearly>1-pSig/2));
scatter(ones(length(TAUCearlyVgat.AUCearly(~indSigVgat)),1)*2,TAUCearlyVgat.AUCearly(~indSigVgat),10,'k','filled')
scatter(ones(length(TAUCearlyVgat.AUCearly(indSigVgat)),1)*2,TAUCearlyVgat.AUCearly(indSigVgat),10,'r','filled')
[h,p]=ttest2(TAUCearlyVglut2.AUCearly,TAUCearlyVgat.AUCearly);
xlim([0,3]);
ylim([0,1]);
plot([1,2],[0.8,0.8],'k-');
text(1,0.9,strcat('p=',num2str(p)));
indSigContraVglut2=logical(TAUCearlyVglut2.pAUCearly<pSig/2);
indSigContraVgat=logical(TAUCearlyVgat.pAUCearly<pSig/2);
text(0.5,1,strcat(num2str(sum(indSigContraVglut2)),'/',num2str(length(indSigContraVglut2)),'contra'));
text(1.5,1,strcat(num2str(sum(indSigContraVgat)),'/',num2str(length(indSigContraVgat)),'contra'));
set(gca,'XTick',[1,2],'XTickLabel',{'vglut2','vgat'});
ylabel('stim AUC');
plot([0,4],[0.5,0.5],'k--');
Tout=Tin;
set(gca,'FontName','Arial','FontSize',12);
end
%table input, plot early AUC comparison
function [Tout] = fScatterAUClateCmp(Tin,pSig)
TAUClateVglut2=Tin(logical((Tin.celltype=='vglut2').*(Tin.Answer == 'correct')),5:6);
indSigVglut2=logical((TAUClateVglut2.pAUClate<pSig/2)+(TAUClateVglut2.pAUClate>1-pSig/2));
scatter(ones(length(TAUClateVglut2.AUClate(~indSigVglut2)),1),TAUClateVglut2.AUClate(~indSigVglut2),10,'k','filled')
hold on;
scatter(ones(length(TAUClateVglut2.AUClate(indSigVglut2)),1),TAUClateVglut2.AUClate(indSigVglut2),10,'r','filled')
TAUClateVgat=Tin(logical((Tin.celltype=='vgat').*(Tin.Answer == 'correct')),5:6);
indSigVgat=logical((TAUClateVgat.pAUClate<pSig/2)+(TAUClateVgat.pAUClate>1-pSig/2));
scatter(ones(length(TAUClateVgat.AUClate(~indSigVgat)),1)*2,TAUClateVgat.AUClate(~indSigVgat),10,'k','filled')
scatter(ones(length(TAUClateVgat.AUClate(indSigVgat)),1)*2,TAUClateVgat.AUClate(indSigVgat),10,'r','filled')
[h,p]=ttest2(TAUClateVglut2.AUClate,TAUClateVgat.AUClate);
xlim([0,3]);
ylim([0,1]);
plot([1,2],[0.8,0.8],'k-');
text(1,0.9,strcat('p=',num2str(p)));
indSigContraVglut2=logical(TAUClateVglut2.pAUClate<pSig/2);
indSigContraVgat=logical(TAUClateVgat.pAUClate<pSig/2);
text(0.5,1,strcat(num2str(sum(indSigContraVglut2)),'/',num2str(length(indSigContraVglut2)),'contra'));
text(1.5,1,strcat(num2str(sum(indSigContraVgat)),'/',num2str(length(indSigContraVgat)),'contra'));
set(gca,'XTick',[1,2],'XTickLabel',{'vglut2','vgat'});
ylabel('late delay AUC');
plot([0,4],[0.5,0.5],'k--');
Tout=Tin;
set(gca,'FontName','Arial','FontSize',12);
end
%calculate correlation between early and late delay
function [Tout]=fScatterCmpEarlyLate(Tin,celltype,pSig)
TAUC=Tin(logical((Tin.celltype==celltype).*(Tin.Answer == 'correct')),3:6);
indSigEarly=logical((TAUC.pAUCearly<pSig/2)+(TAUC.pAUCearly>1-pSig/2));
indSigLate=logical((TAUC.pAUClate<pSig/2)+(TAUC.pAUClate>1-pSig/2));
indSig=logical(indSigEarly.*indSigLate);
curve1=scatter(TAUC.AUCearly(indSig),TAUC.AUClate(indSig),10,'r','filled');
hold on;
curve2=scatter(TAUC.AUCearly(~indSig),TAUC.AUClate(~indSig),10,'k','filled');
text(0,1,['pearson correlation correct=',num2str(corr(TAUC.AUCearly,TAUC.AUClate))]);
xlabel('stim AUC');
ylabel('late delay AUC');
set(gca,'Xlim',[0,1],'Ylim',[0,1],'FontSize',12,'FontName','Arial');
%for error
TAUCerr=Tin(logical((Tin.celltype==celltype).*(Tin.Answer == 'error')),3:6);
indSigEarlyErr=logical((TAUCerr.pAUCearly<pSig/2)+(TAUCerr.pAUCearly>1-pSig/2));
indSigLateErr=logical((TAUCerr.pAUClate<pSig/2)+(TAUCerr.pAUClate>1-pSig/2));
indSig=logical(indSigEarlyErr.*indSigLateErr);
curve3=scatter(TAUCerr.AUCearly(indSig),TAUCerr.AUClate(indSig),10,'r');
hold on;
curve4=scatter(TAUCerr.AUCearly(~indSig),TAUCerr.AUClate(~indSig),10,'k');
text(0,0.9,['pearson correlation error=',num2str(corr(TAUCerr.AUCearly,TAUCerr.AUClate))]);
xlabel('stim AUC');
ylabel('late delay AUC');
set(gca,'Xlim',[0,1],'Ylim',[0,1],'FontSize',12,'FontName','Arial');
% legend([curve1,curve2,curve3,curve4],'significant correct','n.s. correct','significant error','n.s. error');
end
%calculate correlation between early and late delay
function [Tout]=fScatterLineCmpEarlyLate(Tin,celltype,pSig)
TAUC=Tin(logical((Tin.celltype==celltype).*(Tin.Answer == 'correct')),3:6);
indSigEarly=logical((TAUC.pAUCearly<pSig/2)+(TAUC.pAUCearly>1-pSig/2));
indSigLate=logical((TAUC.pAUClate<pSig/2)+(TAUC.pAUClate>1-pSig/2));
curve2=scatter(ones(length(TAUC.AUCearly(~indSigEarly)),1),TAUC.AUCearly(~indSigEarly),10,'k','filled');hold on;
curve1=scatter(ones(length(TAUC.AUCearly(indSigEarly)),1),TAUC.AUCearly(indSigEarly),10,'r','filled');
curve4=scatter(ones(length(TAUC.AUClate(~indSigLate)),1)*2,TAUC.AUClate(~indSigLate),10,'k','filled');
curve3=scatter(ones(length(TAUC.AUClate(indSigLate)),1)*2,TAUC.AUClate(indSigLate),10,'r','filled');
AUCmat2=table2array(TAUC(~indSigEarly,[1,3]));
plot(AUCmat2','color',[0.5,0.5,0.5]);
AUCmat=table2array(TAUC(indSigEarly,[1,3]));
plot(AUCmat','color',[1,0.5,0.5]);
ylabel('AUC');
plot([0,3],[0.5,0.5],'k--');
set(gca,'Xlim',[0,3],'Ylim',[0,1],'FontSize',12,'FontName','Arial');
set(gca,'XTick',[1,2],'XTickLabel',{'stim','late delay'});

end
