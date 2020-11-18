%This code need to be run at >2018b
% dbstop if error
clear;
dbstop if error;
close all;
savepath='H:\FP\summary';
[num,txt,raw] =xlsread('H:\FP\fiber photometry data summary.xlsx');%criteria to choose sessions come from this file
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
fiberstr={'Soma'};%fiberstr={'Soma','Terminal'};
dffstr={'dff baseline correction first','dff motion correction first','dff470','dff410'};
behEventAlign='go cue';%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
behEventSort='first lick';% string can be in{'first lick','reward','go cue'};
masklick='no';%'yes' if mask lickings or 'no' when want to see all activity
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
whichAnimal='';
celltype={'SC vglut2','SC vgat'};
n_siteCell=cell(1,2);%each cell store a cell type
n_datapointCell=cell(1,2);%each cell store a cell type
dataProcessStr=strcat('-grouped by ',selectivitystr{i_selectivity},'-shuffle',num2str(nshuffle),'-lick time1s');%string of data process,used to label file name
%
for i_celltype=1:2
    %select sessions for summary analysis
    experiment=celltype{i_celltype};
    if strcmp(selectivitystr{i_selectivity},'stimuli')%|| strcmp(selectivitystr{i_selectivity},'sensory difficulty')%here only include sessions with probe trials
        ind_session=(T.delay_max==1.5).*(T.delay_min==0.3).*strcmp(T.experiment,experiment).*strcmp(T.used_as_data,'yes').*(T.retract==0).*((contains(T.brain_region,region{i_region}))).*(T.probe_fraction>0);%
    else
        ind_session=(T.delay_max==1.5).*(T.delay_min==0.3).*strcmp(T.experiment,experiment).*strcmp(T.used_as_data,'yes').*(T.retract==0).*((contains(T.brain_region,region{i_region})));%.*(T.probe_fraction>0)
    end
    %ind_session=strcmp(T.animal,whichAnimal).*strcmp(T.used_as_data,'yes').*(T.probe_fraction>0).*(T.retract==0);
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
    %
    [dff_aligned_cat_animal,trialType_cat_animal,behEvent_aligned_cat_animal,licking_aligned_cat_animal,Tsignal_cat_animal]=deal(cell(2,n_site(1)));%two fibers so two cells to store,for multiple sessions, store 1D--1-soma,2-terminal data,2D--sites; note soma sites are more than terminal sites
%     [dff_aligned_by_session,trialType_by_session,behEvent_aligned_by_session,licking_aligned_by_session]=deal(cell(2,n_datapoint(1)));%2D-data points, each means one site one session
    for nfiber=1:length(fiberstr)%here means fiber bilaterally implanted,1-soma, 2-terminal
        %store AUC and p-value
        dataInfoStr=strcat(experiment,'-',regionstr{i_region},'-',fiberstr{nfiber});
        
        fileNameAUC=strcat(savepath,filesep,dataInfoStr,dataProcessStr,'-n',num2str(n_datapoint(nfiber)),'-sites',num2str(n_site(nfiber)),'-EpochAUC.mat');
        datapointNameCell=cell(n_datapoint(nfiber),1);
        varTypes = {'string','categorical','double','double','double','double','double'};
        T_AUC=table('Size',[n_datapoint(nfiber)*2,7],'VariableTypes',varTypes,...
            'VariableNames',{'SessionName','Answer','ITI','sound','delay','response','lick'});
        T_pvalues=table('Size',[n_datapoint(nfiber)*2,7],'VariableTypes',varTypes,...
            'VariableNames',{'SessionName','Answer','ITI','sound','delay','response','lick'});
        siteNameCell=[];%later will be struct with field 'nameSite','nameCases'  
        T_AUC_site=table('Size',[n_site(nfiber)*2,7],'VariableTypes',varTypes,...
            'VariableNames',{'SiteName','Answer','ITI','sound','delay','response','lick'});
        T_pvalues_site=table('Size',[n_site(nfiber)*2,7],'VariableTypes',varTypes,...
            'VariableNames',{'SiteName','Answer','ITI','sound','delay','response','lick'});
        %store ttest p-value
        fileNameTtest=strcat(savepath,filesep,dataInfoStr,dataProcessStr,'-n',num2str(n_datapoint(nfiber)),'-sites',num2str(n_site(nfiber)),'-EpochTtest.mat');
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
                [trialType,behrule,trialTypeStr] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix',fiberSide{i_fiber});%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
%                 %this block choose only easy trials
%                 tempTrialType=zeros(size(trialType,1),2,size(trialType,3));
%                 tempTrialType(:,1,:)=trialType(:,1,:);
%                 tempTrialType(:,2,:)=trialType(:,end,:);
%                 trialType=tempTrialType;
%                 trialTypeStr=strcat(trialTypeStr,'-easy tirals');
                disp([trialTypeStr,'rule of this session is ',behrule]);
                % }
                trialType_cat_animal{nfiber,ind_site}=cat(3,trialType_cat_animal{nfiber,ind_site},trialType);
                trialType_by_session{nfiber,ind_datapoint}=trialType;
                %% concatenate dff, trial type maxtrix and behavior event, etc.
                for ndff=[1] %plot for each dff, see which is better and test whether 410 signals have similar trend
                    T_SigbyEpoch = fGetSigBehEpoch(behEventFrameIndex,dff{1,ndff}(ind_fiber(i_fiber),:),frT);
                    if isempty(Tsignal_cat_animal{nfiber,ind_site})
                        Tsignal_cat_animal{nfiber,ind_site}=T_SigbyEpoch;
                    else
                        Tsignal_cat_animal{nfiber,ind_site}=vertcat(Tsignal_cat_animal{nfiber,ind_site},T_SigbyEpoch);
                    end
                end
                %save corresponding datapoint name- which session and which sites
                fileNameAUC_datapoint=strcat(savepath,filesep,dataInfoStr,'-',datapointNameCell{ind_datapoint,1},dataProcessStr,'-epochAUC.mat');%adding datapoint name to the AUC file name
                stranswer={'correct','error'};
                if exist(fileNameAUC_datapoint,'file')
                    load(fileNameAUC_datapoint);
                    for nResult=1:size(trialType,1)-2
                        T_AUC(2*ind_datapoint+nResult-2,:) = T_AUC_case(nResult,:);
                        T_pvalues(2*ind_datapoint+nResult-2,:) = T_pvalues_case(nResult,:);
                    end
                    if ~strcmp(datapointNameCell{ind_datapoint,1},T_AUC.SessionName(2*ind_datapoint+nResult-2))
                        warning('Name of data point mismatch, double check');
                    end
                end
                for nResult=1:size(trialType,1)
                    %calculate moving p to decide when selectivity becoming significant
                    label = fTrialType2Label(trialType,2);
                    indTrial=trialType(nResult,:,:);
                    indTrial=sum(squeeze(indTrial),1);
                    ind_trial=logical(squeeze(indTrial));
                    poslabel=2;
                    %plot moving AUC, method 1
                    if nResult<size(trialType,1)-1 && nfiber==1 %only calculate AUC from cor/err soma
                        if  ismissing(T_AUC.SessionName(2*ind_datapoint+nResult-2)) % initialized but not assigned value
                            %calculate AUC in different epoch
                            [T_AUC.SessionName(2*ind_datapoint+nResult-2), T_pvalues.SessionName(2*ind_datapoint+nResult-2)]= deal(datapointNameCell{ind_datapoint,1});
                            [T_AUC.Answer(2*ind_datapoint+nResult-2),T_pvalues.Answer(2*ind_datapoint+nResult-2)] = deal(stranswer{nResult});
                            [T_AUC.ITI(2*ind_datapoint+nResult-2),T_pvalues.ITI(2*ind_datapoint+nResult-2)]=fAUC(label(ind_trial),T_SigbyEpoch.ITI(ind_trial),poslabel,1000);
                            [T_AUC.sound(2*ind_datapoint+nResult-2),T_pvalues.sound(2*ind_datapoint+nResult-2)]=fAUC(label(ind_trial),T_SigbyEpoch.sound(ind_trial),poslabel,1000);
                            [T_AUC.delay(2*ind_datapoint+nResult-2),T_pvalues.delay(2*ind_datapoint+nResult-2)]=fAUC(label(ind_trial),T_SigbyEpoch.delay(ind_trial),poslabel,1000);
                            [T_AUC.response(2*ind_datapoint+nResult-2),T_pvalues.response(2*ind_datapoint+nResult-2)]=fAUC(label(ind_trial),T_SigbyEpoch.response(ind_trial),poslabel,1000);
                            [T_AUC.lick(2*ind_datapoint+nResult-2),T_pvalues.lick(2*ind_datapoint+nResult-2)]=fAUC(label(ind_trial),T_SigbyEpoch.lick(ind_trial),poslabel,1000);   
                            T_AUC_case(nResult,:)=T_AUC(2*ind_datapoint+nResult-2,:);
                            T_pvalues_case(nResult,:) = T_pvalues(2*ind_datapoint+nResult-2,:);
                        end
                    end
                end
                save(fileNameAUC_datapoint,'T_AUC_case','T_pvalues_case');   
                ind_datapoint=ind_datapoint+1;
            end
        end
        %calculate AUC for large sessions grouped by sites, only calculated soma, not terminal
        for ind_site=1:n_site(1)
            fileNameAUC_site=strcat(savepath,filesep,dataInfoStr,'-',siteNameCell(ind_site).nameSite,dataProcessStr,'-epochAUC.mat');%adding datapoint name to the AUC file name
            fileNameAUC_site=fileNameAUC_site{1};
            stranswer={'correct','error'};
            if exist(fileNameAUC_site,'file')
                load(fileNameAUC_site);
                %check whether site include same dataset
                if fEqual(nameCases,siteNameCell(ind_site).nameCases)
                    for nResult=1:size(trialType,1)-2
                        T_AUC_site(2*ind_site+nResult-2,:) = T_AUC_site_case(nResult,:);
                        T_pvalues_site(2*ind_site+nResult-2,:) = T_pvalues_site_case(nResult,:);
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
            for nResult=1:2 %only calculated cor/err trials
                label = fTrialType2Label(trialType_cat_animal{1,ind_site},2);
                indTrial=trialType_cat_animal{1,ind_site}(nResult,:,:);
                indTrial=sum(squeeze(indTrial),1);
                ind_trial=logical(squeeze(indTrial));
                poslabel=2;
                if  ismissing(T_AUC_site.SiteName(2*ind_site+nResult-2))
                    %calculate AUC in different epoch
                    [T_AUC_site.SiteName(2*ind_site+nResult-2), T_pvalues_site.SiteName(2*ind_site+nResult-2)]= deal(siteNameCell(ind_site).nameSite);
                    [T_AUC_site.Answer(2*ind_site+nResult-2),T_pvalues_site.Answer(2*ind_site+nResult-2)] = deal(stranswer{nResult});
                    [T_AUC_site.ITI(2*ind_site+nResult-2),T_pvalues_site.ITI(2*ind_site+nResult-2)]=fAUC(label(ind_trial),Tsignal_cat_animal{nfiber,ind_site}.ITI(ind_trial),poslabel,1000);
                    [T_AUC_site.sound(2*ind_site+nResult-2),T_pvalues_site.sound(2*ind_site+nResult-2)]=fAUC(label(ind_trial),Tsignal_cat_animal{nfiber,ind_site}.sound(ind_trial),poslabel,1000);
                    [T_AUC_site.delay(2*ind_site+nResult-2),T_pvalues_site.delay(2*ind_site+nResult-2)]=fAUC(label(ind_trial),Tsignal_cat_animal{nfiber,ind_site}.delay(ind_trial),poslabel,1000);
                    [T_AUC_site.response(2*ind_site+nResult-2),T_pvalues_site.response(2*ind_site+nResult-2)]=fAUC(label(ind_trial),Tsignal_cat_animal{nfiber,ind_site}.response(ind_trial),poslabel,1000);
                    [T_AUC_site.lick(2*ind_site+nResult-2),T_pvalues_site.lick(2*ind_site+nResult-2)]=fAUC(label(ind_trial),Tsignal_cat_animal{nfiber,ind_site}.lick(ind_trial),poslabel,1000);
                    T_AUC_site_case(nResult,:)=T_AUC_site(2*ind_site+nResult-2,:);
                    T_pvalues_site_case(nResult,:) = T_pvalues_site(2*ind_site+nResult-2,:);
                end
             end
            save(fileNameAUC_site,'T_AUC_site_case','T_pvalues_site_case','nameCases');   
        end
        %save AUC, because calculation is time consuming
        save(fileNameAUC,'T_AUC','T_pvalues','datapointNameCell','T_AUC_site','T_pvalues_site','siteNameCell');
        save(fileNameTtest,'pTtestcell','pTtestcellSite');
    end
    %}
    n_datapointCell{i_celltype}=n_datapoint;
    n_siteCell{i_celltype}=n_site;
end

% msgbox('EpochAUC Done');
color_celltype={[1 0.5 0.5],[0.5 0.5 1]};%blue-vgat,red-vglut2
color_celltype_mean={[1 0 0],[0 0 1]};%blue-vglut2,red-vgat
experiment={'SC vglut2','SC vgat'};
celltypestr={'vglut2','vgat'};
ylabelstr={'AUC','p of AUC','p of t-test'};
titlestr={'Correct','Error','Miss','Violation'};
ytext=[1,0.1];
figP=figure;%3*4,1st row-AUC,2nd row-pAUC,3rd row-pTtest
set(gcf, 'PaperPosition', [0 0 4 4]);
ts=double((-frameNum(1):binstep:frameNum(2))*frT/1000);
t_sep_AUC=cell(2,1);%
t_sep_ttest=cell(2,1);

for ncelltype=1:2
    dataInfoStr=strcat(experiment{ncelltype},'-',regionstr{i_region},'-',fiberstr{1});
    fileNameAUC=strcat(savepath,filesep,dataInfoStr,dataProcessStr,'-n',num2str(n_datapointCell{ncelltype}(1)),'-sites',num2str(n_siteCell{ncelltype}(1)),'-EpochAUC.mat');
    load(fileNameAUC);
    fileNameTtest=strcat(savepath,filesep,dataInfoStr,dataProcessStr,'-n',num2str(n_datapointCell{ncelltype}(1)),'-sites',num2str(n_siteCell{ncelltype}(1)),'-EpochTtest.mat');
    load(fileNameTtest);
    pSig=0.05;
    figure(figP);
    subplot(2,2,2*ncelltype-1);%by session
    fplotScatterConnected(T_AUC,T_pvalues,pSig,'correct');
    title([experiment{ncelltype},'-by session']);
    subplot(2,2,2*ncelltype);%by site
    fplotScatterConnected(T_AUC_site,T_pvalues_site,pSig,'correct');
    title([experiment{ncelltype},'-by site']);
%     suptitle(experiment{ncelltype});
end
saveas(figP,[savepath,filesep,'AUCepoch.pdf'],'pdf');
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
%plot AUCepoch change
function[] =fplotScatterConnected(T_AUC,T_pvalues,pSig,Answer)
indrow=(T_AUC.Answer==Answer);
AUCmat=table2array(T_AUC(indrow,3:end));
pAUCmat=table2array(T_pvalues(indrow,3:end));
hold on;
for i=1:size(AUCmat,2)
    indSig=logical((pAUCmat(:,i)<pSig/2)+(pAUCmat(:,i)>1-pSig/2));
    scatter(i*ones(sum(~indSig),1),AUCmat(~indSig,i),10,'k');
    scatter(i*ones(sum(indSig),1),AUCmat(indSig,i),10,'r');
end
plot(AUCmat','color',[0.5,0.5,0.5],'LineWidth',0.5);
ylabel('AUC');
xlim([0,6]);
ylim([0,1.3]);
plot([0,6],[0.5,0.5],'k--');
[h,p_stim_delay]=ttest(AUCmat(:,3)-AUCmat(:,2));
[h,p_delay_lick]=ttest(AUCmat(:,5)-AUCmat(:,3));
[h,p_stim_lick]=ttest(AUCmat(:,5)-AUCmat(:,2));
plot([2,2.9],[1,1],'k-');
text(1,1.1,['p=',num2str(round(p_stim_delay,3))]);
plot([3.1,5],[1,1],'k-');
text(4,1.1,['p=',num2str(round(p_delay_lick,3))]);
plot([2,5],[1.2,1.2],'k-');
text(3,1.3,['p=',num2str(round(p_stim_lick,3))]);

set(gca,'XTick',1:5,'XTickLabel',{'ITI','sound','delay','response','lick'});
%改成旋转后的label以节省空间
xtb = get(gca,'XTickLabel');% 获取横坐标轴标签句柄
xt = get(gca,'XTick');% 获取横坐标轴刻度句柄
yt = get(gca,'YTick'); % 获取纵坐标轴刻度句柄
xtextp=xt;%每个标签放置位置的横坐标，这个自然应该和原来的一样了。
ytextp=-yt(1)*ones(1,length(xt)); % 设置显示标签的位置，写法不唯一，这里其实是在为每个标签找放置位置的纵坐标
% rotation，正的旋转角度代表逆时针旋转，旋转轴可以由HorizontalAlignment属性来设定，
% 有3个属性值：left，right，center，这里可以改这三个值，以及rotation后的角度，这里写的是45
text(xtextp,ytextp,xtb,'HorizontalAlignment','right','rotation',45,'fontsize',12);
set(gca,'xticklabel',[]);% 将原有的标签隐去

set(gca,'FontSize',12);

end