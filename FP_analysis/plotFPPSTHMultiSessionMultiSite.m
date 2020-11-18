%This code choose some sessions and group(but not concanetate) data together 
% to plot PSTH and Ending point trajectories; both show individual sessions
% and animals
% dbstop if error
clear;
dbstop if error;
close all;
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
behEventAlign='delay onset';%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
masklick='yes';%'yes' if mask lickings or 'no' when want to see all activity
stimDur=0.5;%0.5s click duration
if strcmp(behEventAlign,'delay onset')
    frameNumTime=[1,1.5];%from 2s before align point to 5s after align point
else
    frameNumTime=[1,2];%from 2s before align point to 5s after align point
end
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
nshuffle=1000;
p_sepOnset=0.05;%the p-value used to determine when separation happen
%set the criteria to choose sessions
i_region=3;%1-bilateral, 2-unilateral*********variable**************
region={'bilateral','t SC','SC'};%{'bilateral SC','SC'}
regionstr={'bilateral','unilateral','mixed'};
% whichAnimal='pyx241';
celltype={'SC vglut2','SC vgat','ALM terminal'};
n_siteCell=cell(1,2);%each cell store a cell type
n_datapointCell=cell(1,2);%each cell store a cell type
time_before_delay=0.5;
time_before_go=0.5;
dataProcessStr=strcat('-grouped by ',selectivitystr{i_selectivity},'-',combineCorErr,'-binCenteredSize',num2str(binsize),'-alignTo-',behEventAlign,'-timeWindow',num2str(frameNumTime(1)),'-',num2str(frameNumTime(2)),'-EpochWindow stim-',num2str(time_before_delay),'delay-',num2str(time_before_go),'-shuffle',num2str(nshuffle));
% dataProcessStr=strcat('-grouped by ',selectivitystr{i_selectivity},'-binCentered-size',num2str(binsize),'alignTo-',behEventAlign,'-timeWindow',num2str(frameNumTime(1)),'-',num2str(frameNumTime(2)),'-shuffle',num2str(nshuffle));
%-grouped bychoice-binCentered-size3alignTo-delay onset-timeWindow1-1.5-shuffle10

%
for i_celltype=1:length(celltype)
    %select sessions for summary analysis
    experiment=celltype{i_celltype};
    if strcmp(selectivitystr{i_selectivity},'stimuli')%|| strcmp(selectivitystr{i_selectivity},'sensory difficulty')%here only include sessions with probe trials
        ind_session=(T.delay_max==1.5).*(T.delay_min==0.3).*strcmp(T.experiment,experiment).*strcmp(T.used_as_data,'yes').*(T.retract==0).*((contains(T.brain_region,region{i_region}))).*(T.probe_fraction>0);%
    else
        ind_session=(T.delay_max==1.5).*(T.delay_min==0.3).*strcmp(T.experiment,experiment).*strcmp(T.used_as_data,'yes').*(T.retract==0).*((contains(T.brain_region,region{i_region})));%.*(T.probe_fraction>0)
    end
%     ind_session=strcmp(T.animal,whichAnimal).*strcmp(T.used_as_data,'yes').*(strcmp(T.date,'2019/11/30')).*strcmp(T.experiment,experiment);
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
        dataInfoStr=strcat(experiment,'-',regionstr{i_region},'-',fiberstr{nfiber});
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
                [trialType,behrule,trialTypeStr] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix',fiberSide{i_fiber},combineCorErr);%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
%                 %this block choose only easy trials
%                 tempTrialType=zeros(size(trialType,1),2,size(trialType,3));
%                 tempTrialType(:,1,:)=trialType(:,1,:);
%                 tempTrialType(:,2,:)=trialType(:,end,:);
%                 trialType=tempTrialType;
%                 trialTypeStr=strcat(trialTypeStr,'-easy tirals');
                disp([trialTypeStr,',rule of this session is ',behrule]);

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
                %save corresponding datapoint name- which session and which sites
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
                        DatapointName{n_answer*ind_datapoint+nResult-n_answer} = datapointNameCell{ind_datapoint,1};
                        AUCearly(n_answer*ind_datapoint+nResult-n_answer) = AUCcase(nResult).AUCearly;
                        pAUCearly(n_answer*ind_datapoint+nResult-n_answer) = AUCcase(nResult).pAUCearly;
                        AUClate(n_answer*ind_datapoint+nResult-n_answer) = AUCcase(nResult).AUClate;
                        pAUClate(n_answer*ind_datapoint+nResult-n_answer) = AUCcase(nResult).pAUClate;
                        AnswerSession{n_answer*ind_datapoint+nResult-n_answer} = stranswer{nResult};
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
                        if isempty(AUCcell{ind_datapoint,nResult}) 
                            [AUCcell{ind_datapoint,nResult},pAUCcell{ind_datapoint,nResult}] = fMovingAUC(label(indTrial),dff_aligned(indTrial,:),2,nshuffle,binsize,binstep);
                            AUCcase(nResult).AUC=AUCcell{ind_datapoint,nResult};
                            AUCcase(nResult).pAUC=pAUCcell{ind_datapoint,nResult};
                        end
                        if  ~isfield(AUCcase(nResult),'AUCearly') || isnan(pAUClate(n_answer*ind_datapoint+nResult-n_answer))
                            %calculate AUC in different epoch
                            earlyActivity=nanmean(dff_aligned(indTrial,(frameNumTime(1)-time_before_delay)*1000/frT:(frameNumTime(1))*1000/frT),2);%stim and early delay
                            lateActivity=nanmean(dff_aligned(indTrial,(frameNumTime(1)+1.5-time_before_go)*1000/frT:(frameNumTime(1)+1.5)*1000/frT),2);%late delay
                            [AUCcase(nResult).AUCearly,AUCcase(nResult).pAUCearly] = fAUC(label(indTrial),earlyActivity,2,nshuffle);
                            [AUCcase(nResult).AUClate,AUCcase(nResult).pAUClate] = fAUC(label(indTrial),lateActivity,2,nshuffle);
                            AUCearly(n_answer*ind_datapoint+nResult-n_answer) = AUCcase(nResult).AUCearly;
                            pAUCearly(n_answer*ind_datapoint+nResult-n_answer) = AUCcase(nResult).pAUCearly;
                            AUClate(n_answer*ind_datapoint+nResult-n_answer) = AUCcase(nResult).AUClate;
                            pAUClate(n_answer*ind_datapoint+nResult-n_answer) = AUCcase(nResult).pAUClate;
                            DatapointName{n_answer*ind_datapoint+nResult-n_answer} = datapointNameCell{ind_datapoint,1};
                            AnswerSession{n_answer*ind_datapoint+nResult-n_answer} = stranswer{nResult};
                        end
                    end
                    %calculate moving t-test, method 2
                    pTtestcell{ind_datapoint,nResult}=fMovingTtest(label(indTrial),dff_aligned(indTrial,:),binsize,binstep);
                end
                save(fileNameAUC_datapoint,'AUCcase','nameCase');   
                ind_datapoint=ind_datapoint+1;
            end
        end
        %calculate AUC for large sessions grouped by sites, only calculated soma, not terminal
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
                        SiteName{n_answer*ind_site+nResult-n_answer} = siteNameCell(ind_site).nameSite;
                        AUCearlySite(n_answer*ind_site+nResult-n_answer) = AUCsite(nResult).AUCearly;
                        pAUCearlySite(n_answer*ind_site+nResult-n_answer) = AUCsite(nResult).pAUCearly;
                        AUClateSite(n_answer*ind_site+nResult-n_answer) = AUCsite(nResult).AUClate;
                        pAUClateSite(n_answer*ind_site+nResult-n_answer) = AUCsite(nResult).pAUClate;
                        AnswerSite{n_answer*ind_site+nResult-n_answer} = stranswer{nResult};
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
                if isempty(AUCcellSite{ind_site,nResult}) && nResult<=2
                    [AUCcellSite{ind_site,nResult},pAUCcellSite{ind_site,nResult}] = fMovingAUC(label(indTrial),dff_aligned_cat_animal{1,ind_site}(indTrial,:),2,nshuffle,binsize,binstep);
                    AUCsite(nResult).AUC=AUCcellSite{ind_site,nResult};
                    AUCsite(nResult).pAUC=pAUCcellSite{ind_site,nResult};
                end
                if  ~isfield(AUCsite(nResult),'AUCearly') || isnan(AUCearlySite(n_answer*ind_site+nResult-n_answer))
                    %calculate AUC in different epoch
                    earlyActivity=nanmean(dff_aligned_cat_animal{1,ind_site}(indTrial,(frameNumTime(1)-time_before_delay)*1000/frT:(frameNumTime(1))*1000/frT),2);%stim and early delay
                    lateActivity=nanmean(dff_aligned_cat_animal{1,ind_site}(indTrial,(frameNumTime(1)+1.5-time_before_go)*1000/frT:(frameNumTime(1)+1.5)*1000/frT),2);%late delay
                    [AUCsite(nResult).AUCearly,AUCsite(nResult).pAUCearly] = fAUC(label(indTrial),earlyActivity,2,nshuffle);
                    [AUCsite(nResult).AUClate,AUCsite(nResult).pAUClate] = fAUC(label(indTrial),lateActivity,2,nshuffle);
                    AUCearlySite(n_answer*ind_site+nResult-n_answer) = AUCsite(nResult).AUCearly;
                    pAUCearlySite(n_answer*ind_site+nResult-n_answer) = AUCsite(nResult).pAUCearly;
                    AUClateSite(n_answer*ind_site+nResult-n_answer) = AUCsite(nResult).AUClate;
                    pAUClateSite(n_answer*ind_site+nResult-n_answer) = AUCsite(nResult).pAUClate;
                    SiteName{n_answer*ind_site+nResult-n_answer} = siteNameCell(ind_site).nameSite;
                    AnswerSite{n_answer*ind_site+nResult-n_answer} = stranswer{nResult};
                end
                %calculate moving t-test, method 2
                pTtestcellSite{ind_site,nResult}=fMovingTtest(label(indTrial),dff_aligned_cat_animal{1,ind_site}(indTrial,:),binsize,binstep);
            end
            save(fileNameAUC_site,'AUCsite','nameSite','nameCases');   
        end
%         [neuralActivity_by_site,neuralActivity_by_site_raw]=fDFFbyTrialType(dff_aligned_cat_animal,trialType_cat_animal,'z-score');
%         [neuralActivity_by_datapoint,neuralActivity_by_datapoint_raw]=fDFFbyTrialType(dff_aligned_by_session,trialType_by_session,'z-score');
        [neuralActivity_by_site,neuralActivity_by_site_raw]=fDFFbyTrialType(dff_aligned_cat_animal,trialType_cat_animal,'minmax');
        [neuralActivity_by_datapoint,neuralActivity_by_datapoint_raw]=fDFFbyTrialType(dff_aligned_by_session,trialType_by_session,'minmax');

        neuralActivity2Plot=cell(1,2);
        neuralActivity2Plot{1}=neuralActivity_by_site;
        neuralActivity2Plot{2}=neuralActivity_by_datapoint;
        %% plot raster and mean trace for large session
        %
        for ndff=[1]         %  titlestr=strcat(animal_name,'-',fiberstr{nfiber},'-',titlestr);
            if isempty(neuralActivity2Plot)%if no data, then end this loop
                break;
            end
            %         fig=figure;%plot raster
            %         set(gcf, 'position', [0 400 1400 300]);
            figMeanTrace=figure;%plot mean trace
            set(gcf, 'position', [0 0 1400 600]);
            figCases=figure;%plot moving AUC
            set(gcf, 'position', [0 200 1400 600]);
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

                        %                 figure(fig);
                        %                 %raster plot of activity and corresponding licking
                        %                 axActivity=subplot(size(trialType_cat,2),2*size(trialType_cat,1),2*nResult-1+2*size(trialType_cat,1)*(nStim-1));
                        %                 axLicking=subplot(size(trialType,2),2*size(trialType,1),2*nResult+2*size(trialType,1)*(nStim-1));
                        %                 x_lim=[0,sum(frameNum)];
                        %                 fRasterDffandLicking(x_lim,axActivity,neuralActivity,behEvent_aligned_cat,selectedTrialInd,behEventSort,ColLimit,axLicking,licking_aligned_cat);
                        %                 %add labels to activity raster block
                        %                 subplot(axActivity);
                        %                 if nStim == 1
                        %                     title(titlestr{nResult});
                        %                 end
                        %                 if nResult==1&&contains(trialTypeStr,'stimuli')
                        %                     ylabel([num2str(Data_extract.Stimuli(nStim))]);% ' clicks/s'
                        %                 elseif nResult==1&&contains(trialTypeStr,'first lick')
                        %                     ylabelstr={'ipsi lick first','contra lick first'};
                        %                     ylabel(ylabelstr{nStim});
                        %                 elseif nResult==1&&contains(trialTypeStr,'sensory')
                        %                     ylabelstr={'low','high'};
                        %                     ylabel(ylabelstr{nStim});
                        %                 end
                        %                 xlabel(['time(s) from ',behEventAlign]);

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
                    xlabel(['time(s) from ',behEventAlign]);
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
                h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','','contra hard','','contra easy'},'Location','best');
            elseif contains(trialTypeStr,'difficulty')
                h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','contra hard','contra easy'},'Location','best');
            elseif contains(trialTypeStr,'first lick')
                h=legend(curve_meanTrace(:),{'ipsi lick first','contra lick first'},'Location','best');
            elseif contains(trialTypeStr,'sensory')
                h=legend(curve_meanTrace(:),{'low click rate','high click rate'},'Location','best');
            end
            set(h,'box','off');
            text(0.5,0.9,strcat('n=',num2str(n_datapoint),' datapoints',num2str(n_session),' sessions,',num2str(n_site(nfiber)),' site,',num2str(n_animal),' animal'),'FontName','Arial','FontSize',14,'Unit','Normalized');
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
        %}
        n_datapointCell{i_celltype}=n_datapoint;
        n_siteCell{i_celltype}=n_site;
        %save AUC, because calculation is time consuming
        TAUCepochSession=table(string(DatapointName),AUCearly,pAUCearly, AUClate,pAUClate,categorical(AnswerSession),...
            'VariableNames',{'DatapointName' 'AUCearly' 'pAUCearly' 'AUClate' 'pAUClate' 'Answer'});
        TAUCepochSite=table(string(SiteName),AUCearlySite,pAUCearlySite,AUClateSite,pAUClateSite,categorical(AnswerSite),...
            'VariableNames',{'SiteName' 'AUCearly' 'pAUCearly' 'AUClate' 'pAUClate' 'Answer'});
        save(fileNameAUC,'AUCcell','pAUCcell','datapointNameCell','AUCcellSite','pAUCcellSite','siteNameCell','TAUCepochSite','TAUCepochSession');
%         save(fileNameAUC,'AUCcell','pAUCcell','datapointNameCell','AUCcellSite','pAUCcellSite','siteNameCell');
        save(fileNameTtest,'pTtestcell','pTtestcellSite');
        clearvars DatapointName AUCearly pAUCearly AUClate pAUClate AnswerSession AUCearlySite pAUCearlySite AUClateSite pAUClateSite AnswerSite SiteName
    end
end
%}
% {
%% compare different celltype p value(AUC, ttest)
if i_selectivity>=2%only when 2 group exist
    color_celltype={'F16820','646464','46782D'};%vglut2,vgat,ALM
    color_celltype=fHex2RGB(color_celltype);
%     color_celltype={[1 0.5 0.5],[0.5 0.5 1],[0.5,0.5,0.5]};%blue-vgat,red-vglut2
    color_celltype_mean=color_celltype;%blue-vglut2,red-vgat,black-ALM
    experiment={'SC vglut2','SC vgat','ALM terminal'};
    celltypestr={'vglut2','vgat','ALM terminal'};
    ylabelstr={'AUC','p of AUC','p of t-test'};
    titlestr={'Correct','Error','Miss','Violation'};
    ytext=[1,0.1];
    n_celltype_toshow=2;%length(experiment)
%     figP=figure;%3*4,1st row-AUC,2nd row-pAUC,3rd row-pTtest
%     set(gcf, 'position', [0 0 1000 600]);
    ts=double((-frameNum(1):binstep:frameNum(2))*frT/1000);
    t_sep_AUC=cell(2,1);%
    t_sep_ttest=cell(2,1);

    for ncelltype=1:n_celltype_toshow
        dataInfoStr=strcat(experiment{ncelltype},'-',regionstr{i_region},'-',fiberstr{1});
        fileNameAUC=strcat(savepath,filesep,dataInfoStr,dataProcessStr,'-n',num2str(n_datapointCell{ncelltype}(1)),'-sites',num2str(n_siteCell{ncelltype}(1)),'-AUC.mat');
        load(fileNameAUC);
        fileNameTtest=strcat(savepath,filesep,dataInfoStr,dataProcessStr,'-n',num2str(n_datapointCell{ncelltype}(1)),'-sites',num2str(n_siteCell{ncelltype}(1)),'-Ttest.mat');
        load(fileNameTtest);    
        pAUCcell=pAUCcell;%%%%%%%%%%%%here decide to use pAUCcellSite/pAUCcell
        AUCcell=AUCcell;%%%%%%%%%%%%here decide to use AUCcellSite/AUCcell
        pTtestcell=pTtestcell;%%%%%%%%%%%%here decide to use pTtestcell/pTtestcellSite
        p_sig_AUC=cellfun(@(x) (x<p_sepOnset/2)+((1-x)<p_sepOnset/2),pAUCcell,'UniformOutput',0);
        p_sig_ttest=cellfun(@(x) x<p_sepOnset,pTtestcell,'UniformOutput',0);    
        threshold=(frameNumTime(1)-0.5)*1000/frT/binstep;%rule out t_sep before stimuli onset
        t_sep_AUC{ncelltype}=cell(size(p_sig_AUC));
        t_sep_ttest{ncelltype}=cell(size(p_sig_ttest));
        varShowCases='yes';%'no'
        
        if strcmp(varShowCases,'yes')
            figAUCcases=figure;
            n_col=6;
            n_row=ceil(2*size(AUCcell,1)/n_col);
            set(gcf,'PaperPosition',[0,0,1.5*n_col,1.5*n_row]);
            ind_subplot_cor=1:2:size(AUCcell,1)*2;
            ind_subplot_err=2:2:size(AUCcell,1)*2;
            ind_subplot=[ind_subplot_cor;ind_subplot_err];
            for ncase=1:size(AUCcell,1)
                for nResult=1:2
                    tempName=strsplit(datapointNameCell{ncase},'_');
                    tempstr=[celltypestr{ncelltype},'-',tempName{1},'-',tempName{3},'-',titlestr{nResult}];
                    subplot(n_row,n_col,ind_subplot(nResult,ncase));
                    t_sep_AUC{ncelltype}{ncase,nResult}=fOnsetPhaseChangeExamAUCCase(p_sig_AUC{ncase,nResult},AUCcell{ncase,nResult},ts,threshold,tempstr);
                    if ceil(ind_subplot(nResult,ncase)/n_col)==n_row
                        xlabel(['Time(s) from',behEventAlign]);
                    end
                    t_sep_ttest{ncelltype}{ncase,nResult}=fOnsetPhaseChange(p_sig_ttest{ncase,nResult},threshold);
%                     t_sep_ttest{ncelltype}{ncase,nResult}=fOnsetPhaseChangeExamAUCCase(p_sig_ttest{ncase,nResult},pTtestcell{ncase,nResult},ts,threshold,tempstr);
                end
            end
            saveas(figAUCcases,[savepath,filesep,regionstr{i_region},celltypestr{ncelltype},'-algin to ',behEventAlign,dataProcessStr,'-AUC cases.pdf'],'pdf');
        else
            for ncase=1:size(AUCcell,1)
                for nResult=1:2
                    t_sep_AUC{ncelltype}{ncase,nResult}=fOnsetPhaseChange(p_sig_AUC{ncase,nResult},threshold);
                    t_sep_ttest{ncelltype}{ncase,nResult}=fOnsetPhaseChange(p_sig_ttest{ncase,nResult},threshold);
                end
            end
        end
%         for i_data=1:size(AUCcell,1)
%             for icol=1:n_answer
%                 figure(figP);
%                 subplot(3,4,icol);%first row
%                 curve(ncelltype)=plot(ts,AUCcell{i_data,icol},'Color',color_celltype{ncelltype},'LineWidth',0.5);
%                 hold on;
%                 plot(ts,ones(length(ts),1)*0.5,'Color',[0 0 0]);%show 0.5
%                 subplot(3,4,icol+4);%2nd row
%                 curve(ncelltype)=semilogy(ts,pAUCcell{i_data,icol},'Color',color_celltype{ncelltype},'LineWidth',0.5);%here the true p is 1-pvalue
%                 hold on;
%                 semilogy(ts,ones(length(ts),1)*0.01,'Color',[0 0 0]);%show significant level
%                 %text(0,ytext(ncelltype),[experiment{ncelltype},'t_s_e_p=',num2str(ts(t_sep_AUC{ncelltype}{i_data,icol})),'s']);
%                 subplot(3,4,icol+8);%2nd row
%                 curve(ncelltype)=semilogy(ts,pTtestcell{i_data,icol},'Color',color_celltype{ncelltype},'LineWidth',0.5);
%                 hold on;
%                 semilogy(ts,ones(length(ts),1)*0.01,'Color',[0 0 0]);%show significant level
%                 %text(0,ytext(ncelltype),[experiment{ncelltype},'t_s_e_p=',num2str(ts(t_sep_ttest{ncelltype}{i_data,icol})),'s']);
%                 for irow=1:3%label
%                     subplot(3,4,irow*4-3);
%                     ylabel(ylabelstr{irow});
%                 end
%             end
%         end
        iRes=1;
        AUCbySession{ncelltype}=cellfun(@(x) reshape(x(:),1,[]),AUCcell(:,iRes),'UniformOutput' ,false);
        AUCbySession{ncelltype}=cell2mat(AUCbySession{ncelltype});
        AUCbySite{ncelltype}=cellfun(@(x) reshape(x(:),1,[]),AUCcellSite(:,iRes),'UniformOutput' ,false);
        AUCbySite{ncelltype}=cell2mat(AUCbySite{ncelltype});
        Ttemp = fExtendT(TAUCepochSession,celltypestr{ncelltype},'celltype');
        if exist('TAUCepochSessionCombine','var')
            TAUCepochSessionCombine=vertcat(TAUCepochSessionCombine,Ttemp);
        else
            TAUCepochSessionCombine=Ttemp;
        end
        Ttemp = fExtendT(TAUCepochSite,celltypestr{ncelltype},'celltype');
        if exist('TAUCepochSiteCombine','var')
            TAUCepochSiteCombine=vertcat(TAUCepochSiteCombine,Ttemp);
        else
            TAUCepochSiteCombine=Ttemp;
        end
    end
    %save combined table of AUC
    filename = [savepath, filesep,'AUC_table_by_session_summary.csv'];
    writetable(TAUCepochSessionCombine,filename,'WriteVariableNames',true);
    filename = [savepath, filesep,'AUC_table_by_site_summary.csv'];
    writetable(TAUCepochSiteCombine,filename,'WriteVariableNames',true);
%     % plot p value change of AUC/ttest during delay
%     for irow=1:3
%         for icol=1:4
%             figure(figP);
%             subplot(3,4,irow*4-4+icol);
%             set(gca,'FontSize',12,'FontName','Arial','Xlim',[-frameNumTime(1),frameNumTime(2)],'XTick',0:2,'XTickLabel',0:2);
%             plot([0,0],[0,1],'k-');
%             plot([0.5,0.5],[0,1],'k-');
%             if irow==2
%                 x_ttest=fCell2Mat(t_sep_AUC{1}(:,icol));
%                 y_ttest=fCell2Mat(t_sep_AUC{2}(:,icol));
%                 [~,p_t_sep]=ttest2(x_ttest,y_ttest);
%                 text(0,0.01,['p=',num2str(p_t_sep)]);
%             end
%             if irow==3
%                 x_ttest=fCell2Mat(t_sep_ttest{1}(:,icol));
%                 y_ttest=fCell2Mat(t_sep_ttest{2}(:,icol));
%                 [~,p_t_sep]=ttest2(x_ttest,y_ttest);
%                 text(0,0.01,['p=',num2str(p_t_sep)]);
%             end
%         end
%     end
%     subplot(3,4,4);
%     legend(curve_mean,experiment);
%     saveas(figP,[savepath,filesep,regionstr{i_region},'-algin to ',behEventAlign,'-pAUC_Ttest.pdf'],'pdf');
        
    %plot AUC separately for vglut2/vgat/AML terminal
    figAUC=figure;
    set(gcf,'paperPosition',[0 0 2*n_celltype_toshow 4]);
    for icol=1:n_celltype_toshow%each col one cell type
        figure(figAUC);
        subplot(2,n_celltype_toshow,icol);%1st row
%         AUCbySession{icol}=abs(AUCbySession{icol}-0.5)+0.5;
        plot(ts,AUCbySession{icol},'Color',color_celltype{icol},'LineWidth',0.5); hold on;
        curve(icol)=fPlotMean_CI(ts,AUCbySession{icol},0.05,color_celltype_mean{icol});        
        titlestrAUC=strcat(celltypestr{icol},'-by-session');
        title(titlestrAUC);
        box off;
        subplot(2,n_celltype_toshow,icol+n_celltype_toshow);%2nd row
%         AUCbySite{icol}=abs(AUCbySite{icol}-0.5)+0.5;
        plot(ts,AUCbySite{icol},'Color',color_celltype{icol},'LineWidth',0.5);   hold on;
        curve(icol)=fPlotMean_CI(ts,AUCbySite{icol},0.05,color_celltype_mean{icol});
        titlestrAUC=strcat(celltypestr{icol},'-by-site');
        title(titlestrAUC);
        box off;
        for irow=1:2
            subplot(2,n_celltype_toshow,icol+(irow-1)*n_celltype_toshow);%2nd row
            plot(ts,ones(length(ts),1)*0.5,'Color',[0 0 0]);%show 0.5
            plot([0,0],[0,1],'k-');
            if strcmp(behEventAlign,'stim onset')
                plot([0.5,0.5],[0,1],'k-');
            elseif strcmp(behEventAlign,'delay onset')
                plot([-0.5,-0.5],[0,1],'k');
            end
            ylabel('AUC');
            xlabel(['Time (s) from ',behEventAlign]);
            xlim([-frameNumTime(1),frameNumTime(2)]);
            set(gca,'FontName','Arial','FontSize',12);
        end  
    end
    saveas(figAUC,[savepath,filesep,'AUC-algin to-',behEventAlign,'.pdf'],'pdf');

    %scatter plot of epoch AUC
    mkr={'h','>','s'};
    pSig=0.05;
    figAUCepoch=figure;
    set(gcf,'paperPosition',[0,0,10,8]);
    subplot(4,5,1);%1d row by session,1st col late AUC comparison
    [Tout] = fScatterAUCCmpOneEpoch(TAUCepochSessionCombine,pSig,'early');
    title('by session');
    subplot(4,5,2);%1d row by session,1st col late AUC comparison
    [Tout] = fScatterAUCCmpOneEpoch(TAUCepochSessionCombine,pSig,'late');
    title('by session');
    subplot(4,5,11);
    [Tout] = fScatterAUCCmpOneEpoch(TAUCepochSiteCombine,pSig,'early');
    title('by site');
    subplot(4,5,12);
    [Tout] = fScatterAUCCmpOneEpoch(TAUCepochSiteCombine,pSig,'late');
    title('by site');
    subplot(4,5,3);
    fScatterLineCmpEarlyLate(TAUCepochSessionCombine,'vglut2',pSig);
    title('vglut2');
    subplot(4,5,4);
    fScatterLineCmpEarlyLate(TAUCepochSessionCombine,'vgat',pSig);
    title('vgat');
    subplot(4,5,5);
    fScatterLineCmpEarlyLate(TAUCepochSessionCombine,'ALM terminal',pSig);
    title('ALM terminal');
    subplot(4,5,6);
    fBarShiftedAUC(TAUCepochSessionCombine,celltypestr,color_celltype_mean,n_celltype_toshow);
    subplot(4,5,7);
    fScatterCmpEarlyLateCellType(TAUCepochSessionCombine,celltypestr,color_celltype_mean,mkr);
    subplot(4,5,8);
    fScatterCmpEarlyLate(TAUCepochSessionCombine,'vglut2',pSig);  
    title('vglut2');
    subplot(4,5,9);
    fScatterCmpEarlyLate(TAUCepochSessionCombine,'vgat',pSig);
    title('vgat');
    subplot(4,5,10);
    fScatterCmpEarlyLate(TAUCepochSessionCombine,'ALM terminal',pSig);
    title('ALM terminal');
    subplot(4,5,13);
    fScatterLineCmpEarlyLate(TAUCepochSiteCombine,'vglut2',pSig);
    title('vglut2');
    subplot(4,5,14);
    fScatterLineCmpEarlyLate(TAUCepochSiteCombine,'vgat',pSig);
    title('vgat');
    subplot(4,5,15);
    fScatterLineCmpEarlyLate(TAUCepochSiteCombine,'ALM terminal',pSig);
    title('ALM terminal');
    subplot(4,5,16);
    fBarShiftedAUC(TAUCepochSiteCombine,celltypestr,color_celltype_mean,n_celltype_toshow);
    subplot(4,5,17);
    fScatterCmpEarlyLateCellType(TAUCepochSiteCombine,celltypestr,color_celltype_mean,mkr);
    subplot(4,5,18);
    fScatterCmpEarlyLate(TAUCepochSiteCombine,'vglut2',pSig);
    title('vglut2');
    subplot(4,5,19);
    fScatterCmpEarlyLate(TAUCepochSiteCombine,'vgat',pSig);
    title('vgat');
    subplot(4,5,20);
    fScatterCmpEarlyLate(TAUCepochSiteCombine,'ALM terminal',pSig);
    title('ALM terminal');
    saveas(figAUCepoch,[savepath,filesep,regionstr{i_region},'-epoch AUC comparison.pdf'],'pdf');
    
%     %consistency of selectivity from during late stim to late delay
% 
%     figCI=figure;
%     set(gcf,'paperPosition',[0,0,4,2]);
%     subplot(1,2,1);%by session;    
%     fplotAUCconsitency(TAUCepochSessionCombine,'correct',pSig);
%     subplot(1,2,2);%by site;
%     fplotAUCconsitency(TAUCepochSiteCombine,'correct',pSig);
%     saveas(figCI,[savepath,filesep,regionstr{i_region},'-epoch AUC consistency.pdf'],'pdf');
%     
%     %scatter plot of different cell type onset time
%     figT=figure;
%     set(gcf, 'position', [0 0 600 400]);
%     tlabelstr='time of separation';
%     for irow=1:2%t_AUC/t_ttest
%         for icol=1:2%cor/err
%             if irow==1
%                 x_ttest=fCell2Mat(t_sep_AUC{1}(:,icol));
%                 y_ttest=fCell2Mat(t_sep_AUC{2}(:,icol));
%             else
%                 x_ttest=fCell2Mat(t_sep_ttest{1}(:,icol));
%                 y_ttest=fCell2Mat(t_sep_ttest{2}(:,icol));
%             end
%             subplot(2,2,icol+irow*2-2);
%             [~,p_t_sep]=ttest2(x_ttest,y_ttest);
%             scatter(ones(size(x_ttest)),ts(x_ttest),20,color_celltype_mean{1});hold on;
%             scatter(1+ones(size(y_ttest)),ts(y_ttest),20,color_celltype_mean{2});
%             fPlotMeanSem(ts(x_ttest),0.8,color_celltype_mean{1});
%             fPlotMeanSem(ts(y_ttest),2.2,color_celltype_mean{2});
%             y_lim=get(gca,'Ylim');
%             x_lim=[0,3];
%             plot([1,2],[1,1]*y_lim(end)*0.9,'k-');
%             plot(x_lim,[0,0],'k--');
%             plot(x_lim,[-0.5,-0.5],'k--');
%             text(1.5,y_lim(end),['p=',num2str(p_t_sep)]);
%             title(titlestr{icol});
%             ylabel(tlabelstr);            
%             set(gca,'XTick',[1,2],'XTickLabel',celltypestr,'Ylim',[-frameNumTime(1),frameNumTime(end)],'FontSize',12,'FontName','Arial');
%         end
%     end
%     
%     saveas(figT,[savepath,filesep,regionstr{i_region},'-algin to ',behEventAlign,'-Time_Separate.pdf'],'pdf');
end
%}
writetable(TAUCepochSessionCombine,[savepath,filesep,'EpochAUC.csv']);
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
function [matdata,matdataraw]=fDFFbyTrialType(dff,trialtypecell,method)
trialtype=trialtypecell{1};
matdata=cell(size(trialtype,1),size(trialtype,2));
matdataraw=cell(size(trialtype,1),size(trialtype,2));
for nResult=1:size(trialtype,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
    for nStim=1:size(trialtype,2) %for each stimulus%[1,6]%just 2 end trials
        matdataraw{nResult,nStim}=zeros(length(dff),size(dff{1},2));
        for idatapoint=1:length(dff)
            selectedTrialInd=trialtypecell{1,idatapoint}(nResult,nStim,:);
            selectedTrialInd=logical(squeeze(selectedTrialInd))';
            rawdata=dff{1,idatapoint}(selectedTrialInd,:);
            matdataraw{nResult,nStim}(idatapoint,:)=nanmean(rawdata,1);
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
function [Tout] = fScatterAUCCmpOneEpoch(Tin,pSig,epoch)
if strcmp(epoch,'early')
    ind_col=3:4;
    ylabelstr='stim AUC';
elseif strcmp(epoch,'late')
    ind_col=5:6;
    ylabelstr='late delay AUC';
end
TAUCVglut2=Tin(logical((Tin.celltype=='vglut2').*(Tin.Answer == 'correct')),ind_col);
TAUCVglut2=table2array(TAUCVglut2);
indSigVglut2=logical((TAUCVglut2(:,2)<pSig/2)+(TAUCVglut2(:,2)>1-pSig/2));
scatter(ones(length(TAUCVglut2(~indSigVglut2,1)),1),TAUCVglut2(~indSigVglut2,1),10,'k','filled')
hold on;
scatter(ones(length(TAUCVglut2(indSigVglut2,1)),1),TAUCVglut2(indSigVglut2,1),10,'r','filled')
TAUCVgat=Tin(logical((Tin.celltype=='vgat').*(Tin.Answer == 'correct')),ind_col);
TAUCVgat=table2array(TAUCVgat);
indSigVgat=logical((TAUCVgat(:,2)<pSig/2)+(TAUCVgat(:,2)>1-pSig/2));
scatter(ones(length(TAUCVgat(~indSigVgat,1)),1)*2,TAUCVgat(~indSigVgat,1),10,'k','filled')
scatter(ones(length(TAUCVgat(indSigVgat,1)),1)*2,TAUCVgat(indSigVgat,1),10,'r','filled')
% TAUCALM=Tin(logical((Tin.celltype=='ALM terminal').*(Tin.Answer == 'correct')),ind_col);
% TAUCALM=table2array(TAUCALM);
% indSigALM=logical((TAUCALM(:,2)<pSig/2)+(TAUCALM(:,2)>1-pSig/2));
% scatter(ones(length(TAUCALM(~indSigALM,1)),1)*3,TAUCALM(~indSigALM,1),10,'k','filled')
% scatter(ones(length(TAUCALM(indSigALM,1)),1)*3,TAUCALM(indSigALM,1),10,'r','filled')
% % [h,p1]=ttest2(TAUCVglut2(:,1),TAUCVgat(:,1));
% % [h,p2]=ttest2(TAUCVglut2(:,1),TAUCALM(:,1));
% % [h,p3]=ttest2(TAUCALM(:,1),TAUCVgat(:,1));
[h,varp1]=vartest2(TAUCVglut2(:,1),TAUCVgat(:,1));
% [h,varp2]=vartest2(TAUCVglut2(:,1),TAUCALM(:,1));
% [h,varp3]=vartest2(TAUCALM(:,1),TAUCVgat(:,1));
p1=ranksum(TAUCVglut2(:,1),TAUCVgat(:,1));
% p2=ranksum(TAUCVglut2(:,1),TAUCALM(:,1));
% p3=ranksum(TAUCALM(:,1),TAUCVgat(:,1));
xlim([0,3]);
ylim([0,1]);
plot([1,2],[0.8,0.8],'k-');
text(1,0.9,[plabelsymbol(p1),'var',plabelsymbol(varp1)]);
% plot([2.1,3],[0.8,0.8],'k-');
% text(2.5,0.9,plabelsymbol(p2));
% plot([1,3],[0.9,0.9],'k-');
% text(1.5,1,plabelsymbol(p3));
indSigContraVglut2=logical(TAUCVglut2(:,2)<pSig/2);
indSigContraVgat=logical(TAUCVgat(:,2)<pSig/2);
% indSigContraALM=logical(TAUCALM(:,2)<pSig/2);
text(0.5,0.1,strcat(num2str(sum(indSigContraVglut2)),'/',num2str(length(indSigContraVglut2)),'contra'));
text(1.5,0.1,strcat(num2str(sum(indSigContraVgat)),'/',num2str(length(indSigContraVgat)),'contra'));
% text(2.5,0.1,strcat(num2str(sum(indSigContraALM)),'/',num2str(length(indSigContraALM)),'contra'));
% set(gca,'XTick',[1,2,3],'XTickLabel',{'vglut2','vgat','ALM'});
set(gca,'XTick',[1,2],'XTickLabel',{'vglut2','vgat'});
ylabel(ylabelstr);
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
% text(0,1,['pearson correlation correct=',num2str(corr(TAUC.AUCearly,TAUC.AUClate))]);
xlabel('stim AUC');
ylabel('late delay AUC');
set(gca,'Xlim',[0,1],'Ylim',[0,1],'FontSize',12,'FontName','Arial');
plot([0,1],[0.5,0.5],'k--');
plot([0.5,0.5],[0,1],'k--');
%for error
TAUCerr=Tin(logical((Tin.celltype==celltype).*(Tin.Answer == 'error')),3:6);
indSigEarlyErr=logical((TAUCerr.pAUCearly<pSig/2)+(TAUCerr.pAUCearly>1-pSig/2));
indSigLateErr=logical((TAUCerr.pAUClate<pSig/2)+(TAUCerr.pAUClate>1-pSig/2));
indSig=logical(indSigEarlyErr.*indSigLateErr);
curve3=scatter(TAUCerr.AUCearly(indSig),TAUCerr.AUClate(indSig),10,'r');
hold on;
curve4=scatter(TAUCerr.AUCearly(~indSig),TAUCerr.AUClate(~indSig),10,'k');
% text(0,0.9,['pearson correlation error=',num2str(corr(TAUCerr.AUCearly,TAUCerr.AUClate))]);
xlabel('stim AUC');
ylabel('late delay AUC');
set(gca,'Xlim',[0,1],'Ylim',[0,1],'FontSize',12,'FontName','Arial');
plot([0,1],[0.5,0.5],'k--');
plot([0.5,0.5],[0,1],'k--');
% legend([curve1,curve2,curve3,curve4],'significant correct','n.s. correct','significant error','n.s. error');
end
%cmp AUC early vs late for different cell type
function [Tout]=fScatterCmpEarlyLateCellType(Tin,celltype,color,mkr)
for icelltype=length(celltype):-1:1
    TAUC=Tin(logical((Tin.celltype==celltype{icelltype}).*(Tin.Answer == 'correct')),3:6);
    curve_scatter(icelltype)=scatter(TAUC.AUCearly,TAUC.AUClate,10,color{icelltype},mkr{icelltype});
    hold on;
end
% hl=legend(curve_scatter(:),'vglut2' , 'vgat' , 'ALM terminal','AutoUpdate','off');
hl=legend(curve_scatter(:),'vglut2' , 'vgat','AutoUpdate','off');
set(hl,'Box','Off');
% text(0,1,['pearson correlation correct=',num2str(corr(TAUC.AUCearly,TAUC.AUClate))]);
xlabel('stim AUC');
ylabel('late delay AUC');
set(gca,'Xlim',[0.3,0.7],'Ylim',[0.05,0.95],'FontSize',12,'FontName','Arial');
plot([0,1],[0.5,0.5],'k--');
plot([0.5,0.5],[0,1],'k--');
end
%cmp AUC early vs late for different cell type, histogram of shifted
%session
function [Tout] = fBarShiftedAUC(Tin,celltype,color,n_celltype_toshow)
plot([0,n_celltype_toshow+1],[0.5,0.5],'k-');hold on;
for icelltype=1:length(celltype(1:n_celltype_toshow))
    TAUC=Tin(logical((Tin.celltype==celltype{icelltype}).*(Tin.Answer == 'correct')),3:6);
    temp=(TAUC.AUCearly-0.5).*(TAUC.AUClate-0.5);
    flag=(temp>0);
    bar(icelltype,sum(flag)/length(flag),'FaceColor',color{icelltype});
    text(icelltype,0.1,[num2str(sum(flag)),'/',num2str(length(flag))]);
end
set(gca,'Xlim',[0,n_celltype_toshow+1],'Ylim',[0,1],'XTick',1:n_celltype_toshow,'XTickLabel',celltype(1:n_celltype_toshow));
box off;
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
%calculate index of consistency
function [Tout] = fplotAUCconsitency(Tin,answer,pSig)
consistentIndexVglut2 =fAUCconsitency(Tin,answer,'vglut2',pSig);
consistentIndexVgat =fAUCconsitency(Tin,answer,'vgat',pSig);
scatter(ones(length(consistentIndexVglut2),1),consistentIndexVglut2,10,'r');hold on;
scatter(ones(length(consistentIndexVgat),1)*2,consistentIndexVgat,10,'b');
% [h,p]=ttest2(consistentIndexVglut2,consistentIndexVgat);

[h,p2]=vartest2(consistentIndexVglut2,consistentIndexVgat);
if h
    p=ranksum(consistentIndexVglut2,consistentIndexVgat);
else
    [~,p]=ttest2(consistentIndexVglut2,consistentIndexVgat);
end
plot([0,3],[0,0],'k--');
set(gca,'XTick',[1,2],'XTickLabel',{'vglut2','vgat'});
ylabel('consistency index');
xlim([0,3]);
% ylim([0,1]);
plot([1,2],[0.8,0.8],'k-');
text(0.5,0.9,strcat('p=',num2str(p)),'Unit','Normalized');
text(0.5,0.95,strcat('vartest2 p=',num2str(p2)),'Unit','Normalized');
set(gca,'FontSize',12,'FontName','Arial');
Tout=Tin;
end
function [consistentIndex] =fAUCconsitency(Tin,answer,celltype,pSig)
Tauc=Tin(logical((Tin.Answer==answer).*(Tin.celltype==celltype)),:);
% %nx method
% indSig=logical((Tauc.pAUCearly<pSig/2)+(Tauc.pAUCearly>1-pSig/2));
% matAUC=[Tauc.AUCearly(indSig),Tauc.AUClate(indSig)];
% consistentIndex=(1-abs(matAUC(:,2)-matAUC(:,1)))./(1+abs(matAUC(:,2)-matAUC(:,1)));
% consistentIndex=(1-(matAUC(:,2)-matAUC(:,1)))./(1+(matAUC(:,2)-matAUC(:,1)));
%pyx method
matAUC=[Tauc.AUCearly,Tauc.AUClate];
% indSig1=logical((Tauc.pAUCearly<pSig/2)+(Tauc.pAUCearly>1-pSig/2));
% indSig2=logical((Tauc.pAUClate<pSig/2)+(Tauc.pAUClate>1-pSig/2));
% indSig=logical(indSig1);
% matAUC=[Tauc.AUCearly(indSig),Tauc.AUClate(indSig)];
% consistentIndex=((matAUC(:,2)-0.5)./(matAUC(:,1)-0.5));
% directional consistentIndex
consistentIndex=(abs(matAUC(:,2)-0.5)./(matAUC(:,2)-0.5)).*(abs(matAUC(:,1)-0.5)./(matAUC(:,1)-0.5)).*abs(matAUC(:,2)-matAUC(:,1));
end