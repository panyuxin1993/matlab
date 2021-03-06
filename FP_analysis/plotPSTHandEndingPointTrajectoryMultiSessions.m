%This code choose some sessions and group(but not concanetate) data together 
% to plot PSTH and Ending point trajectories; both show individual sessions
% and animals
% dbstop if error
clear;
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
selectivitystr={'stimuli','difficulty','sensory','choice'};%sensory means grouping difficulties;
trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
fiberstr={'Soma'};%fiberstr={'Soma','Terminal'};
dffstr={'dff baseline correction first','dff motion correction first','dff470','dff410'};
behEventAlign='stimOnset';%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
masklick='yes';%'yes' if mask lickings or 'no' when want to see all activity
stimDur=0.5;%0.5s click duration
frameNumTime=[0.5,2.2];%from 2s before align point to 5s after align point
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
for i_celltype=1:2
    %select sessions for summary analysis
    experiment=celltype{i_celltype};
    if strcmp(selectivitystr{i_selectivity},'stimuli')|| strcmp(selectivitystr{i_selectivity},'difficulty')%here only include sessions with probe trials
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
    [dff_aligned_cat_animal,trialType_cat_animal,behEvent_aligned_cat_animal,licking_aligned_cat_animal]=deal(cell(2,n_site(1)));%two fibers so two cells to store,for multiple sessions, store 1D--1-soma,2-terminal data,2D--sites; note soma sites are more than terminal sites
    [dff_aligned_by_session,trialType_by_session,behEvent_aligned_by_session,licking_aligned_by_session]=deal(cell(2,n_datapoint(1)));%2D-data points, each means one site one session
    for nfiber=1:length(fiberstr)%here means fiber bilaterally implanted,1-soma, 2-terminal
        %store AUC and p-value
        fileNameAUCrecent=strcat(savepath,filesep,experiment,'-recentAUCfile');
        fileNameAUC=strcat(savepath,filesep,experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-binCentered-size',num2str(binsize),'-n',num2str(n_datapoint(nfiber)),'-sites',num2str(n_site(nfiber)),'-shuffle',num2str(nshuffle),'-AUC.mat');
        if exist(fileNameAUC,'file')
            load(fileNameAUC);
        elseif exist(fileNameAUCrecent,'file')%if no such file, start from recent file and update.
            load(fileNameAUCrecent);
            if size(AUCcell,1)<n_datapoint(nfiber)%have new data, so need to update cell array
                tempcell=cell(n_datapoint(nfiber)-size(AUCcell,1),4);
                tempNameCell=cell(n_datapoint(nfiber)-size(AUCcell,1),1);
                AUCcell=cat(1,AUCcell,tempcell);
                pAUCcell=cat(1,pAUCcell,tempcell);
                datapointNameCell=cat(1,datapointNameCell,tempNameCell);
            end
            if size(AUCcellSite,1)<n_site(nfiber)
                tempcell=cell(n_site(nfiber)-size(AUCcellSite,1),4);
                tempNameCell=cell(n_site(nfiber)-size(AUCcellSite,1),1);
                AUCcellSite=cat(1,AUCcellSite,tempcell);
                pAUCcellSite=cat(1,pAUCcellSite,tempcell);
                siteNameCell=cat(1,siteNameCell,tempNameCell);
            end
        else
            AUCcell=cell(n_datapoint(nfiber),4);%each cell 1d-data points, 2d-one result(cor/err/miss/vio)
            pAUCcell=cell(n_datapoint(nfiber),4);
            datapointNameCell=cell(n_datapoint(nfiber),1);
            AUCcellSite=cell(n_site(nfiber),4);
            pAUCcellSite=cell(n_site(nfiber),4);
            siteNameCell=cell(n_site(nfiber),1);
        end
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
                %name this site
                tempSiteName=strcat(T.animal(ind_session(i)),'-',fiberSide{i_fiber});
                tempIndSite=cellfun(@(x) strcmp(x,tempSiteName),siteNameCell);
                tempNSite=cellfun(@isempty, siteNameCell);
                tempNSite=sum(~tempNSite);
                if sum(tempIndSite)>0
                    ind_site=find(tempIndSite);
                else
                    ind_site=tempNSite+1;%create a new name in a new empty cell
                    siteNameCell{ind_site}=tempSiteName;
                end
                %get trial type
                [trialType,behrule,trialTypeStr] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix',fiberSide{i_fiber});%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
                trialType_cat_animal{nfiber,ind_site}=cat(3,trialType_cat_animal{nfiber,ind_site},trialType);
                trialType_by_session{nfiber,ind_datapoint}=trialType;
                %% concatenate dff, trial type maxtrix and behavior event, etc.
                for ndff=[1] %plot for each dff, see which is better and test whether 410 signals have similar trend
                    if strcmp(behEventAlign,'stimOnset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
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
                for nResult=1:size(trialType,1)%plot raster plot of DLC body parts
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                        end              
                    end
                    %calculate moving t-test, method 2
                    pTtestcell{ind_datapoint,nResult}=fMovingTtest(label(indTrial),dff_aligned(indTrial,:),binsize,binstep);                  
                end
                %save corresponding datapoint name- which session and which sites
                tempdatapointNameCell=strsplit(rootpath,filesep);
                datapointNameCell{ind_datapoint,1}=strcat(tempdatapointNameCell{end},'_',fiberSide{i_fiber});
                ind_datapoint=ind_datapoint+1;
            end
        end
        %calculate AUC for large sessions grouped by sites, only calculated soma, not terminal
        for ind_site=1:n_site(1)
            for nResult=1:4 %only calculated cor/err trials
                label = fTrialType2Label(trialType_cat_animal{1,ind_site},2);
                indTrial=trialType_cat_animal{1,ind_site}(nResult,:,:);
                indTrial=sum(squeeze(indTrial),1);
                indTrial=logical(squeeze(indTrial));
                if isempty(AUCcellSite{ind_site,nResult}) && nResult<=2
                    [AUCcellSite{ind_site,nResult},pAUCcellSite{ind_site,nResult}] = fMovingAUC(label(indTrial),dff_aligned_cat_animal{1,ind_site}(indTrial,:),2,nshuffle,binsize,binstep);
                end
                %calculate moving t-test, method 2
                pTtestcellSite{ind_site,nResult}=fMovingTtest(label(indTrial),dff_aligned_cat_animal{1,ind_site}(indTrial,:),binsize,binstep);
            end

        end
        neuralActivity_by_site=fDFFbyTrialType(dff_aligned_cat_animal,trialType_cat_animal);
        neuralActivity_by_datapoint=fDFFbyTrialType(dff_aligned_by_session,trialType_by_session);
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
%             figMovingAUC=figure;%plot moving AUC
%             set(gcf, 'position', [0 200 1400 600]);
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
                    else
                        color_mean_trace=fMeanTraceColor(behrule,size(trialType,2));
                    end

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
                    
                    %plot mean trace
                    figure(figMeanTrace);%save mean trace
                    subplot(2,size(neuralActivity2Plot{figrow},1),nResult+(figrow-1)*size(neuralActivity2Plot{figrow},1));
                    neuralActivity=neuralActivity2Plot{figrow}{nResult,nStim};
                    curve_meanTrace(nStim)=fPSTHBinDelayEndPoint(neuralActivity,color_mean_trace{nStim},'EndSmooth','smooth','EndSmoothMethod','bin','ParaEndSmoothMethod',0.3*1000/frT);
                    if ~exist('yrange','var')
                        y_lim=get(gca,'Ylim');
                    else
                        y_lim=yrange;
                    end
%                     set(gca,'Ylim',yrange);%for comparing different cell type
                    plot(frameNumTime(1)*round(1000/frT)*ones(1,2)+1,[y_lim(1),y_lim(2)],'k-');%align to a behavior event
                    if strcmp(behEventAlign,'stimOnset')
                        plot((frameNumTime(1)+stimDur)*round(1000/frT)*ones(1,2)+1,[y_lim(1),y_lim(2)],'k-');%if align to stimOnset, also plot delay onset
                    end
                    xlabel(['time(s) from ',behEventAlign]);
                    title(strcat(titlestr{nResult},titlestr_row{figrow}));
                    set(gca,'xtick',[round(1000/frT*(frameNumTime(1)-ceil(frameNumTime(1))))+1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
                    x_lim=get(gca,'Xlim');
                    xlim([x_lim(1),size(neuralActivity,2)]);
                    %   ylim([nanmin(nanmin(neuralActivityGrouped)),nanmax(nanmax(neuralActivityGrouped))]);
                    set(gca,'FontName','Arial','FontSize',14);
                    box off;
                    if nResult==1
                        ylabel('normalized \it\DeltaF/F');
                    end
                end
%                 %calculate moving p to decide when separation happens
%                 label = fTrialType2Label(trialType_cat,2);
%                 indTrial=trialType_cat(nResult,:,:);
%                 indTrial=sum(squeeze(indTrial),1);
%                 indTrial=logical(squeeze(indTrial));
%                 ts=double((-frameNum(1)+binsize/2:binstep:frameNum(2)-binsize/2+1)*frT/1000);
%                 datapointNameCell{end,1}='large session';
%                 %plot moving AUC, method 1
%                 if nResult<size(trialType_cat,1)-1 && nfiber==1 %only calculate AUC from cor/err soma
%                     if isempty(AUCcell{end,nResult})
%                         [AUCcell{end,nResult},pAUCcell{end,nResult}] = fMovingAUC(label(indTrial),dff_aligned_cat_animal{nfiber}(indTrial,:),2,nshuffle,binsize,binstep);
%                     end
%                     figure(figMovingAUC);%save mean trace
%                     subplot(1,size(trialType_cat,1)-2,nResult);
%                     plot(ts,AUCcell{end,nResult},'-k');
%                     %text(0,1,['Time of separation is',num2str(tp)]);
%                 end
%                 %calculate moving t-test, method 2
%                 pTtestcell{end,nResult}=fMovingTtest(label(indTrial),dff_aligned_cat_animal{nfiber}(indTrial,:),binsize,binstep);
%                 %label time of separation on the PSTH figure
%                 t_sep=fOnsetPhaseChange(pTtestcell{end,nResult}<p_sepOnset);
%                 t_sepAUC=fOnsetPhaseChange(pAUCcell{end,nResult}<p_sepOnset);
%                 figure(figMeanTrace);%save mean trace
%                 subplot(1,size(trialType_cat,1),nResult);
%                 text(x_lim(1),y_lim(1),['time of separation=',num2str(ts(t_sep)),'s(t-test); ',num2str(ts(t_sepAUC)),'s(AUC)']);
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
            text(x_lim(1),y_lim(end),strcat('n=',num2str(n_datapoint),' datapoints',num2str(n_session),' sessions,',num2str(n_site(nfiber)),' site,',num2str(n_animal),' animal'),'FontName','Arial','FontSize',14);
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
                %             saveas(fig,[savepath,filesep,experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-sort ',behEventSort,'-masklick-',masklick,'-raster.pdf'],'pdf');
            end
%             saveas(figMovingAUC,[savepath,filesep,experiment,'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-algin to ',behEventAlign,'-n',num2str(n_datapoint(nfiber)),'-AUC.pdf'],'pdf');
%             close all;
        end
        n_datapointCell{i_celltype}=n_datapoint;
        n_siteCell{i_celltype}=n_site;
        %save AUC, because calculation is time consuming
        save(fileNameAUC,'AUCcell','pAUCcell','datapointNameCell','AUCcellSite','pAUCcellSite','siteNameCell');
        save(fileNameAUCrecent,'AUCcell','pAUCcell','datapointNameCell','AUCcellSite','pAUCcellSite','siteNameCell');
        save(fileNameTtest,'pTtestcell','pTtestcellSite');
    end
end
%}
%% compare different celltype p value(AUC, ttest)
if i_selectivity>2%only when 2 group exist
    color_celltype={[1 0.5 0.5],[0.5 0.5 1]};%blue-vgat,red-vglut2
    color_celltype_mean={[1 0 0],[0 0 1]};%blue-vglut2,red-vgat
    experiment={'SC vglut2','SC vgat'};
    celltypestr={'vglut2','vgat'};
    ylabelstr={'AUC','p of AUC','p of t-test'};
    titlestr={'Correct','Error','Miss','Violation'};
    ytext=[1,0.1];
    figP=figure;%3*4,1st row-AUC,2nd row-pAUC,3rd row-pTtest
    set(gcf, 'position', [0 0 1000 600]);
    figAUC=figure;
    set(gcf,'position',[0 0 600 600]);
    ts=double((-frameNum(1)+binsize/2:binstep:frameNum(2)-binsize/2+1)*frT/1000);
    t_sep_AUC=cell(2,1);%
    t_sep_ttest=cell(2,1);
    for ncelltype=1:2
        fileNameAUC=strcat(savepath,filesep,experiment{ncelltype},'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-binCentered-size',num2str(binsize),'-n',num2str(n_datapointCell{ncelltype}(nfiber)),'-sites',num2str(n_siteCell{ncelltype}(nfiber)),'-shuffle',num2str(nshuffle),'-AUC.mat');
        load(fileNameAUC);
        fileNameTtest=strcat(savepath,filesep,experiment{ncelltype},'-',regionstr{i_region},'-',fiberstr{nfiber},'-grouped by ',selectivitystr{i_selectivity},'-binCentered-size',num2str(binsize),'-n',num2str(n_datapointCell{ncelltype}(nfiber)),'-sites',num2str(n_siteCell{ncelltype}(nfiber)),'-Ttest.mat');
        load(fileNameTtest);    
%         pAUCcell=pAUCcellSite;%%%%%%%%%%%%here decide to use pAUCcellSite/pAUCcell
%         AUCcell=AUCcellSite;%%%%%%%%%%%%here decide to use AUCcellSite/AUCcell
%         pTtestcell=pTtestcellSite;%%%%%%%%%%%%here decide to use pTtestcell/pTtestcellSite
        tempstr=cell(size(AUCcell));
        tempts=cell(size(AUCcell));
        for i=1:length(titlestr)
            for j=1:size(AUCcell,1)%for each session
                tempstr{j,i}=[celltypestr{ncelltype},'-',num2str(j),'-',titlestr{i}];
            end
        end
        [tempts{:}]=deal(ts);
        threshold=zeros(size(AUCcell))+(frameNum(1)-binsize/2)/binstep;%rule out t_sep before stimuli onset
        threshold=num2cell(threshold);
        p_sig_AUC=cellfun(@(x) (x<p_sepOnset/2)+((1-x)<p_sepOnset/2),pAUCcell,'UniformOutput',0);
        p_sig_ttest=cellfun(@(x) x<p_sepOnset,pTtestcell,'UniformOutput',0);
        t_sep_AUC{ncelltype}=cellfun(@fOnsetPhaseChangeExamAUCCase,p_sig_AUC,AUCcell,tempts,threshold,tempstr,'UniformOutput',0);
%         t_sep_ttest{ncelltype}=cellfun(@fOnsetPhaseChangeExamAUCCase,p_sig_ttest,pTtestcell,tempts,threshold,tempstr,'UniformOutput',0);
%         t_sep_AUC{ncelltype}=cellfun(@fOnsetPhaseChange,p_sig_AUC,threshold,'UniformOutput',0);
        t_sep_ttest{ncelltype}=cellfun(@fOnsetPhaseChange,p_sig_ttest,threshold,'UniformOutput',0);
        for i_data=1:size(AUCcell,1)
            for icol=1:2
                figure(figP);
                subplot(3,4,icol);%first row
                curve(ncelltype)=plot(ts,AUCcell{i_data,icol},'Color',color_celltype{ncelltype},'LineWidth',0.5);
                hold on;
                plot(ts,ones(length(ts),1)*0.5,'Color',[0 0 0]);%show 0.5
                subplot(3,4,icol+4);%2nd row
                curve(ncelltype)=semilogy(ts,pAUCcell{i_data,icol},'Color',color_celltype{ncelltype},'LineWidth',0.5);%here the true p is 1-pvalue
                hold on;
                semilogy(ts,ones(length(ts),1)*0.01,'Color',[0 0 0]);%show significant level
                %text(0,ytext(ncelltype),[experiment{ncelltype},'t_s_e_p=',num2str(ts(t_sep_AUC{ncelltype}{i_data,icol})),'s']);
                subplot(3,4,icol+8);%2nd row
                curve(ncelltype)=semilogy(ts,pTtestcell{i_data,icol},'Color',color_celltype{ncelltype},'LineWidth',0.5);
                hold on;
                semilogy(ts,ones(length(ts),1)*0.01,'Color',[0 0 0]);%show significant level
                %text(0,ytext(ncelltype),[experiment{ncelltype},'t_s_e_p=',num2str(ts(t_sep_ttest{ncelltype}{i_data,icol})),'s']);
                for irow=1:3%label
                    subplot(3,4,irow*4-3);
                    ylabel(ylabelstr{irow});
                end
            end
        end
        iRes=1;
        AUCbySession{ncelltype}=cellfun(@(x) reshape(x(:),1,[]),AUCcell(:,iRes),'UniformOutput' ,false);
        AUCbySession{ncelltype}=cell2mat(AUCbySession{ncelltype});
        AUCbySite{ncelltype}=cellfun(@(x) reshape(x(:),1,[]),AUCcellSite(:,iRes),'UniformOutput' ,false);
        AUCbySite{ncelltype}=cell2mat(AUCbySite{ncelltype});
    end
    %plot AUC separately for vglut2/vgat
    for icol=1:2%each col one cell type
        figure(figAUC);
        subplot(2,2,icol);%1st row
        plot(ts,AUCbySession{icol},'Color',color_celltype{icol},'LineWidth',0.5); hold on;
        curve(icol)=fPlotMean_CI(ts,AUCbySession{icol},0.05,color_celltype_mean{icol});        
        titlestrAUC=strcat(celltypestr{icol},'-by-session');
        title(titlestrAUC);
        subplot(2,2,icol+2);%2nd row
        plot(ts,AUCbySite{icol},'Color',color_celltype{icol},'LineWidth',0.5);   hold on;
        curve(icol)=fPlotMean_CI(ts,AUCbySite{icol},0.05,color_celltype_mean{icol});
        titlestrAUC=strcat(celltypestr{icol},'-by-site');
        title(titlestrAUC);
        for irow=1:2
            subplot(2,2,icol+irow*2-2);%2nd row
            plot(ts,ones(length(ts),1)*0.5,'Color',[0 0 0]);%show 0.5
            plot([0,0],[0,1],'k-');
            plot([0.5,0.5],[0,1],'k-');
            ylabel('AUC');
            xlabel('time(s) from stimuli onset');
            set(gca,'FontName','Arial','FontSize',12);
        end  
    end
    saveas(figAUC,[savepath,filesep,'AUC-algin to-',behEventAlign,'.pdf'],'pdf');
    for irow=1:3
        for icol=1:4
            figure(figP);
            subplot(3,4,irow*4-4+icol);
            set(gca,'FontSize',12,'FontName','Arial','Xlim',[-frameNumTime(1),frameNumTime(2)],'XTick',0:2,'XTickLabel',0:2);
            plot([0,0],[0,1],'k-');
            plot([0.5,0.5],[0,1],'k-');
            if irow==2
                x_ttest=fCell2Mat(t_sep_AUC{1}(:,icol));
                y_ttest=fCell2Mat(t_sep_AUC{2}(:,icol));
                [~,p_t_sep]=ttest2(x_ttest,y_ttest);
                text(0,0.01,['p=',num2str(p_t_sep)]);
            end
            if irow==3
                x_ttest=fCell2Mat(t_sep_ttest{1}(:,icol));
                y_ttest=fCell2Mat(t_sep_ttest{2}(:,icol));
                [~,p_t_sep]=ttest2(x_ttest,y_ttest);
                text(0,0.01,['p=',num2str(p_t_sep)]);
            end
        end
    end
%     subplot(3,4,4);
%     legend(curve_mean,experiment);
    saveas(figP,[savepath,filesep,regionstr{i_region},'-algin to ',behEventAlign,'-pAUC_Ttest.pdf'],'pdf');
    %scatter plot of different cell type onset time
    figT=figure;
    set(gcf, 'position', [0 0 600 400]);
    tlabelstr='time of separation';
    for irow=1:2%t_AUC/t_ttest
        for icol=1:2%cor/err
            if irow==1
                x_ttest=fCell2Mat(t_sep_AUC{1}(:,icol));
                y_ttest=fCell2Mat(t_sep_AUC{2}(:,icol));
            else
                x_ttest=fCell2Mat(t_sep_ttest{1}(:,icol));
                y_ttest=fCell2Mat(t_sep_ttest{2}(:,icol));
            end
            subplot(2,2,icol+irow*2-2);
            [~,p_t_sep]=ttest2(x_ttest,y_ttest);
            scatter(ones(size(x_ttest)),ts(x_ttest),20,color_celltype_mean{1});hold on;
            scatter(1+ones(size(y_ttest)),ts(y_ttest),20,color_celltype_mean{2});
            fPlotMeanSem(ts(x_ttest),0.8,color_celltype_mean{1});
            fPlotMeanSem(ts(y_ttest),2.2,color_celltype_mean{2});
            y_lim=get(gca,'Ylim');
            x_lim=[0,3];
            plot([1,2],[1,1]*y_lim(end)*0.9,'k-');
            plot(x_lim,[0,0],'k--');
            plot(x_lim,[0.5,0.5],'k--');
            text(1.5,y_lim(end),['p=',num2str(p_t_sep)]);
            title(titlestr{icol});
            ylabel(tlabelstr);            
            set(gca,'XTick',[1,2],'XTickLabel',celltypestr,'Ylim',[-frameNumTime(1),frameNumTime(end)],'FontSize',12,'FontName','Arial');
        end
    end
    
    saveas(figT,[savepath,filesep,regionstr{i_region},'-algin to ',behEventAlign,'-Time_Separate.pdf'],'pdf');
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
function [matdata]=fDFFbyTrialType(dff,trialtypecell)
trialtype=trialtypecell{1};
matdata=cell(size(trialtype,1),size(trialtype,2));
for nResult=1:size(trialtype,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
    for nStim=1:size(trialtype,2) %for each stimulus%[1,6]%just 2 end trials
        matdata{nResult,nStim}=zeros(length(dff),size(dff{1},2));
        for idatapoint=1:length(dff)
            selectedTrialInd=trialtypecell{1,idatapoint}(nResult,nStim,:);
            selectedTrialInd=logical(squeeze(selectedTrialInd))';
            rawdata=dff{1,idatapoint}(selectedTrialInd,:);
            matdata{nResult,nStim}(idatapoint,:)=nanmean(rawdata,1);
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

%shadow as ci
xpatch=[ts(1:size(neuralActivity,2)), fliplr(ts(1:size(neuralActivity,2)))];
ypatch=[neuralActivityCI(1,:),fliplr(neuralActivityCI(2,:))];
p=patch(xpatch,ypatch,color);%plot confidence interval
p.FaceAlpha=0.05;
p.EdgeColor='none';%color;%'none';
hold on;

outputcurve=plot(ts(1:size(neuralActivity,2)),neuralActivityMean,'Color',color,'linewidth',2);
hold on;
end