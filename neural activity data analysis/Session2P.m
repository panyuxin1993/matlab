classdef Session2P
    %SESSION2P store and process 2P data from a session
    %   Detailed explanation goes here
    
    properties
        name        %session name
        metadata    %metadata, such as frame rate
        path        %important path of files
        filename    %file names of important files
        dff         %dff file
        zscored_dff %z-scored dff
        spkr        %spiking rate file based on deconvolution from dff
        zscored_spkr %z-scored spiking rate file based on deconvolution from dff
        Data_extract %behavior data file
        behaviorEvent%preprocessed behavioral data
        SavedCaTrials
        DLC         %struct read from a .xls file, including dcnum, dctxt
        behavior_performance%behavior performance stored in a table
        TSVM        %table storing the cross epoch and time SVM
        TAUC        %table storing the cross epoch and time AUC
        Tmean       %table storing the cross epoch and time mean activities
    end
    
    methods
        function obj = Session2P(sessionname,rootpath,trial2include,trial2exclude,varargin)
            %SESSION2P Construct an instance of this class
            %   input is some important path, and some metadata, such as
            %   session name, trials to include/exclude
            obj.name=sessionname;
            temp=strsplit(sessionname,'_');
            obj.metadata.animal=temp{1};
            obj.metadata.date=temp{2};
            obj.metadata.datetime=datetime(str2double(obj.metadata.date(1:4)),str2double(obj.metadata.date(5:6)),str2double(obj.metadata.date(7:8)));
            if length(temp)>2
                obj.metadata.field=temp{3};
            else
                obj.metadata.field=[];
            end
            obj.metadata.trial2include=trial2include;
            obj.metadata.trial2exclude=trial2exclude;
            
            if ~isempty(varargin)
                obj.metadata.ind_tr_1=varargin{1};
            else
                obj.metadata.ind_tr_1=1;
            end
            obj.path.root=rootpath;
            dirmat=strcat(rootpath,filesep,'*.mat');
            dirs=dir(dirmat);
            dircell=struct2cell(dirs);
            filenames=dircell(1,:);
            file_imaging=cellfun(@(x) contains(x,'CaTrialsSIM_'), filenames);
            file_beh=cellfun(@(x) contains(x,'Virables'), filenames);
            if sum(file_beh)==0 %if still no Data_extract variable
                obj.filename.behdata=fDataExtract(1,obj.path.root,'*imaging.mat');%my SC data names
%                 obj.filename.behdata =fDataExtract(1,obj.path.root,'*2AFC.mat'); %CD M2 data names
            else
                obj.filename.behdata = filenames{file_beh};
            end
            obj.path.behdata=strcat(rootpath,filesep,obj.filename.behdata);
            load(obj.path.behdata);
            obj.Data_extract=Data_extract;

            %calcium activity data and extract related metadata
            obj.filename.Ca = filenames{file_imaging};
            obj.path.Ca=strcat(rootpath,filesep,obj.filename.Ca);
            load(obj.path.Ca);
            obj.SavedCaTrials=SavedCaTrials;
            obj.metadata.nROI=size(SavedCaTrials.f_raw{1},1);%ROI number
            obj.metadata.frT = SavedCaTrials.FrameTime;%frame time
            obj.metadata.ntr=length(SavedCaTrials.f_raw);%trial number of f_raw
            obj.metadata.nFrameEachTrial=cellfun(@(x) size(x,2),SavedCaTrials.f_raw);%frame number each trials
            ind_1stFrame=zeros(1,length(obj.metadata.nFrameEachTrial));
            ind_1stFrame(1)=1;
            ind_1stFrame(2:end)=cumsum(obj.metadata.nFrameEachTrial(1:end-1))+1;
            obj.metadata.ind_1stFrame=ind_1stFrame(obj.metadata.ind_tr_1:obj.metadata.ind_tr_1+obj.metadata.ntr-1);%if number of trials unsed for analysis is not whole but part of trials
            
            %get index of trial used for analyses
            obj.metadata.indTrial2include=fIncludeTrials(trial2include,ind_1stFrame,'logical');
            obj.metadata.indTrialAfterExcluded=fExcludeTrials(trial2exclude,ind_1stFrame,'logical');
            obj.metadata.indTrial2use=logical(obj.metadata.indTrial2include.*obj.metadata.indTrialAfterExcluded);
            dirs_struct=dir(rootpath);
            ind_movement_file=arrayfun(@(x) strcmp(x.name,'movement_trials.xlsx'),dirs_struct);
            if sum(ind_movement_file)>0 %exist the file with movement
                T_ind_movement=readtable([rootpath,filesep,'movement_trials.xlsx']);
                obj.metadata.indTrialMovementExcluded=fExcludeTrials(T_ind_movement.Ind_trial_movement,ind_1stFrame,'logical');
                obj.metadata.indTrial2use=logical(obj.metadata.indTrial2use.*obj.metadata.indTrialMovementExcluded);
            end

            %dff
            obj.filename.dff = 'dff.mat';
            if exist(strcat(obj.path.root,filesep,obj.filename.dff),'file')
                load(strcat(obj.path.root,filesep,obj.filename.dff));
                disp(['already have dff data and load ',obj.name]);
            else
                temp=strsplit(obj.path.root,filesep);
                savefolder=temp{end-2};
                dff=obj.mGetDff(savefolder,'save figure');
            end
            obj.dff=dff;
            obj.filename.dff_field = 'dff_NPseg.mat';
            obj.filename.zscored_dff = 'zscored_dff.mat';
            if exist(strcat(obj.path.root,filesep,obj.filename.zscored_dff),'file')
                load(strcat(obj.path.root,filesep,obj.filename.zscored_dff));
            else
                zscored_dff=obj.mGetZscoredDff();    
            end
            obj.zscored_dff=zscored_dff;
            %deconvolved firing rate
            obj.filename.spkr='deconvolution.mat';
            if exist(strcat(obj.path.root,filesep,obj.filename.spkr),'file')
                load(strcat(obj.path.root,filesep,obj.filename.spkr));
                disp(['already have decovolution data and load ',obj.name]);
            else
                flag_plot='plot_result';
                [spiking_rate, denoised_trace, zscored_spkr] = obj.mGetDeconvolvedFR(flag_plot);
            end
            obj.spkr=spiking_rate;
            obj.zscored_spkr=zscored_spkr;
            
            %preprocess of behavioral data, e.g. aligned to calcium signals
            % align to behavior event
            [obj.behaviorEvent.behEventFrameIndex,obj.behaviorEvent.lickingFrameIndex] = fGetBehEventTime( obj.Data_extract, ind_1stFrame, obj.SavedCaTrials.FrameTime ,obj.metadata.ind_tr_1);%get behavior event time

            %video related information
            if exist([rootpath,filesep,'video'],'dir')
                dirmat=strcat(rootpath,filesep,'video',filesep,'*.avi');
                dirs=dir(dirmat);
                dircell=struct2cell(dirs);
                filenames=dircell(1,:);
                file_video=cellfun(@(x) ~contains(x,'crop'), filenames);
                obj.filename.video = filenames{file_video};
                obj.filename.cropvideo=strrep(obj.filename.video,'.','-crop-mjpeg.');
                obj.filename.OLED = strrep(obj.filename.video,'.avi','-OLED.mat');
                path_OLED=strcat(rootpath,filesep,'video',filesep,obj.filename.OLED);
                v_reader=VideoReader(strcat(rootpath,filesep,'video',filesep,obj.filename.cropvideo));%video may have wrong frame rate, but cropped video will not,since it was transformed using ffmpeg
                video_frameRate=v_reader.FrameRate;
                obj.metadata.video_frameRate=video_frameRate;
                if ~exist(path_OLED,'file')
                    fileOLEDcsv=strrep(obj.filename.video,'.avi','-OLED.csv');
                    T=readtable(strcat(rootpath,filesep,'video',filesep,fileOLEDcsv));%different sessions may need mannually set parameters and save result separately
                    timeTrialStartBeh=cellfun(@fTrialStartStr2double,Data_extract.Time_trialStart);% time(s) relative to start of session day
                    obj.metadata.frameTrialStartVideo=fTrialStartByOLEDrefArduino( T, video_frameRate ,timeTrialStartBeh);
                    save(path_OLED,'obj.metadata.frameTrialStartVideo');
                else
                    load(path_OLED);
                    obj.metadata.frameTrialStartVideo=frameTrialStartVideo;
                end
                dirmat=strcat(rootpath,filesep,'video',filesep,'*.h5');
                dirs=dir(dirmat);
                dircell=struct2cell(dirs);
                filenames=dircell(1,:);
                obj.filename.DLCcoordinate = strrep(filenames{1},'h5','mat');
                obj.path.DLCcoordinate=strcat(rootpath,filesep,'video',filesep,obj.filename.DLCcoordinate );
                file_trace=strrep(obj.filename.DLCcoordinate,'mat','csv');
                if ~exist(obj.path.DLCcoordinate,'file')
                    [dcnum,dctxt,~] = xlsread(file_trace,1);%this process is time consuming, so not do very often
                    save(obj.path.DLCcoordinate,'dcnum','dctxt');
                else
                    load(obj.path.DLCcoordinate)
                end
                obj.DLC.dcnum=dcnum;
                obj.DLC.dctxt=dctxt;
            else % no corresponding video file
                obj.filename.video = '';
                obj.filename.cropvideo='';
                obj.filename.OLED = '';
                obj.metadata.video_frameRate=[];
                obj.metadata.frameTrialStartVideo=[];
                obj.filename.DLCcoordinate = '';
                obj.path.DLCcoordinate='';
                obj.DLC.dcnum=[];
                obj.DLC.dctxt=[];
            end
        end
        
        function [Tperformance,obj]=mSummaryBehavior(obj,f)
            %summary behavior
            %f- ref fGetTrialType, way to summary behavior
            [ trialType ,obj.metadata.sensory_motor_rule] = fGetTrialType( obj.Data_extract,[] ,f ,'matrix');
            switch f
                case 1
                    if strcmp(obj.metadata.sensory_motor_rule,'low click rate-right')
                        trialType=flip(trialType,2);%so now, the order is ipsi to contra
                    end
                    stimuli=num2cell(sort(unique(obj.Data_extract.Stimulus))');
                    rowName=cellfun(@num2str,stimuli,'UniformOutput',false);
                case 3 
                    rowName={'ipsi','contra'};
                case 4
                    rowName={'ipsi','contra'};
            end
            rowName{end+1}='total';
            [var_miss,var_vio,var_cor]=deal(zeros(size(trialType,2)+1,1));
            for i_stim=1:size(trialType,2)
                var_miss(i_stim)=sum(trialType(3,i_stim,:))/sum(trialType(:,i_stim,:),'all');
                var_vio(i_stim)=sum(trialType(4,i_stim,:))/(sum(trialType(:,i_stim,:),'all')-sum(trialType(3,i_stim,:)));
                var_cor(i_stim)=sum(trialType(1,i_stim,:))/(sum(trialType(1,i_stim,:))+sum(trialType(2,i_stim,:)));
            end
            var_miss(end)=sum(trialType(3,:,:),'all')/size(trialType,3);
            var_vio(end)=sum(trialType(4,:,:),'all')/(size(trialType,3)-sum(trialType(3,:,:),'all'));
            var_cor(end)=sum(trialType(1,:,:),'all')/(sum(trialType(1,:,:),'all')+sum(trialType(2,:,:),'all'));
            obj.behavior_performance=table(var_cor,var_vio,var_miss,...
                'VariableNames',{'correct','violation','miss'},...
                'RowNames',rowName);
            Tperformance=obj.behavior_performance;
        end
        
        function [dff,obj] = mGetDff(obj,savefolder,varargin)
            savestr='dff';
            figHow='save';
            if ~isempty(varargin)
                if strcmp(varargin{1},'save figure')
                    figHow='save';
                else
                    figHow='see';
                end
                if length(varargin)>1 && strcmp(varargin{2},'field_NPseg')
                    savestr='dff_NPseg';
                end
            end

            if ~exist([obj.path.root,filesep,savefolder])
                mkdir([obj.path.root,filesep,savefolder]);
            end
            
            f_cat  = [];
            ind_1stFr(1) = 1;
            ind_tr_1=1;
            ntr = length(obj.SavedCaTrials.f_raw); %cancate all trials together, or  just set the number of trials to use
            if strcmp(savestr,'dff')
                for i = ind_tr_1:ntr
                    f_cat = [f_cat obj.SavedCaTrials.f_raw{i}];
                    ind_1stFr(i+1) = size(f_cat,2) + 1;
                end
                ind_1stFr(i+1) = [];
            elseif strcmp(savestr,'dff_NPseg')
                if ~isfield(obj.SavedCaTrials,'SegNPdataAll')
                    warning([savefolder,'no field named SegNPdataAll, so skip it and set dff empty']);
                    dff=[];
                    obj.dff=[];
                    return;
                end
                for i = ind_tr_1:ntr
                    f_cat = [f_cat nanmean(obj.SavedCaTrials.SegNPdataAll{i}([6,7,10,11],:),1)];
                    ind_1stFr(i+1) = size(f_cat,2) + 1;
                end
                ind_1stFr(i+1) = [];
            end
            dff=zeros(size(f_cat));
            for roiNo = 1:size(f_cat,1) %SavedCaTrials.nROIs may be not true
                %%
                %substract changing baseline
                totalFr = size(f_cat,2);
                frT = obj.SavedCaTrials.FrameTime;
                span = round(20*1000/frT/2); %every 40s, get a mean and substrate to compensate for baseline change
                ind_x = 0;
                x = [];
                for i = 1:totalFr
                    %     ind_x = ind_x + 1;
                    ind1 = i- span;
                    ind2 = i+span;
                    if ind1 <= 0
                        ind1 = 1;
                    end
                    if ind2 > totalFr
                        ind2 = totalFr;
                    end
                    x(i) = prctile(f_cat(roiNo,ind1:ind2),5);
                end
                f = f_cat(roiNo,:) - x + mean(x);
                % figure; histogram(f_cat_sub(roiNo,:),100)
                %%
                %get f mode as f0
                [N,edges,bin] = histcounts(f,100);
                f0 = edges(N == max(N));
                if length(f0)==2
                    f0=mean(f0);
                end
                %get dff
                dff(roiNo,:) = (f - f0)/f0;
                if exist('figHow','var')
                    B=figure;
                    set(gcf, 'position', [0 0 1500 500]);
                    plot(f_cat(roiNo,:),'k-');
                    hold on;
                    plot(x,'c-');
                    plot(f,'g-');
                    plot(dff(roiNo,:),'b-');
                    legend('f cat','moving f mean','f baseline correction','dff');
                    set(gca,'FontName','Arial','FontSize',14);
                    if strcmp(figHow,'save')
                        saveas(B,[obj.path.root,filesep,savefolder,filesep,'ROI-',num2str(roiNo),'-cat_f.jpg'],'jpg');
                    end
                    close all;
                end
            end
            obj.dff=dff;
            save([obj.path.root,filesep,savestr],'dff');
        end
        
        function zscored_dff=mGetZscoredDff(obj)
            %smooth the dff and calculated a z-scored dff
            nROI=size(obj.dff,1);
            zscored_dff=zeros(size(obj.dff));
            for roiNo = 1:nROI
                smoothed_sr=smooth(obj.dff(roiNo,:),1);%smooth the data
                sr_mean=nanmean(smoothed_sr);%mean for each roi
                sr_std=nanstd(smoothed_sr);%std for each roi
                zscored_dff(roiNo,:)=(smoothed_sr-sr_mean)/sr_std;
            end
            save([obj.path.root, filesep,'zscored_dff.mat'],  'zscored_dff');
            
        end
        
        function [spiking_rate, denoised_trace, zscored_spkr] = ...
                mGetDeconvolvedFR(obj,flag_plot)
            Paras={[],[],1000/obj.metadata.frT,[]};

            [spiking_rate, denoised_trace, zscored_spkr,ROIFitCoefs] = fxy_GetDeconvolvedFR(obj.dff,Paras,obj.name,obj.path.root,flag_plot);
        end
        
        function fr_smoothed=mSmoothFR(obj,activity_type,binsize)
            %smooth the firing rate data, binsize (in ms)
            switch activity_type
                case 'dff'
                    activities=obj.dff;
                case 'spkr'
                    activities=obj.spkr;
                case 'zscored_dff'
                    activities=obj.zscored_dff;
                case 'zscored_spkr'
                    activities=obj.zscored_spkr;
            end
            span=ceil((binsize/obj.metadata.frT-1)/2);
            fr_smoothed=zeros(size(activities));
            for t=1:size(activities,2)
                t1=max(1,t-span);
                t2=min(size(activities,2),t+span);
                fr_smoothed(:,t)=nanmean(activities(:,t1:t2),2);
            end
        end
        
        function population_stability=mGetPopulationStability(obj,activity_type,tau)
           %using dff/spkr etc. to calcultate a population trajectory and
           %calcualte the angles of vectors pointing to different points at the
           %trajectory
           %activity_type- dff/spkr etc, decide which activities to be used
           %tau- time between 2 ypetime points, which defined the vector and
           %used to calculate the angle
           %population_stability- 1-by-n vector, n is the total time points
           %of input activities.
           switch activity_type
               case 'dff' 
                   activities=obj.dff;
               case 'spkr'
                   activities=obj.spkr;
               case 'zscored_dff'
                   activities=obj.zscored_dff;
               case 'zscored_spkr'
                   activities=obj.zscored_spkr;
           end
           population_stability=zeros(1,size(activities,2));
           for t=1:size(activities,2)
               t1=min(1,t-tau);
               t2=max(size(activities,2),t+tau);
               v1=activities(:,t1);
               v2=activities(:,t2);
               population_stability(t)=dot(v1,v2)/(norm(v1)*norm(v2));
           end
        end
        
        function bodyco = mGetDLCcoordinate(obj,bodyparts,coordinates,...
                treshold4likelihood, varargin)
            %get coordinate of body parts, method analogy to dff

            indcol1=cellfun(@(x) strcmp(bodyparts,x),obj.DLC.dctxt(2,:));
            indcol2=cellfun(@(x) strcmp(coordinates,x),obj.DLC.dctxt(3,:));
            indcol_likelihood=cellfun(@(x) strcmp('likelihood',x),obj.DLC.dctxt(3,:));
            indcol_co=logical(indcol1.*indcol2);
            indcol_like=logical(indcol1.*indcol_likelihood);
            bodyco=obj.DLC.dcnum(:,indcol_co);
            bodycoli=obj.DLC.dcnum(:,indcol_like);
            bodyco(bodycoli<treshold4likelihood)=nan;%rule out those low likelihood data
            figtemp=figure;
            subplot(2,1,1);
            plot(obj.DLC.dcnum(:,indcol_co),'k');hold on;
            plot(bodyco,'r'); 
            if ~strcmp(coordinates,'likelihood')
                if contains(bodyparts,'LickPort')
                    bodyco=fBaselineCorrection(bodyco,5*obj.metadata.video_frameRate);%40s as span
                end
                bodyco=bodyco-nanmean(bodyco);%calculate pixel shift, if it is likelihood, no need for normalization
                bodyco(isnan(bodyco))=0;%transform nan to 0
%                 % rule out those body parts with low likelihood
%                 if mean(bodycoli)<treshold4likelihood
%                     bodyco(:)=nan;
%                 end
            end
            plot(bodyco,'b');
            legend('raw','high likelihood','normalized');
            subplot(2,1,2);
            histogram(bodycoli,'BinWidth',0.1);
            xlabel('likelihood');
            close(figtemp);
        end
        
        function curve_meanTrace = mPlotDffPSTH(obj,frameNumTime,...
                behEventAlign,masklick,i_selectivity,...
                savename_figdff,title_fig)
            %mPlotDffPSTH plot PSTH of dff this session
            %   using other existed function
            obj.metadata.behEventAlign=behEventAlign;
            obj.metadata.masklick=masklick;
            obj.metadata.i_selectivity=i_selectivity;
            curve_meanTrace = fPlotDffPSTH_ASession(obj.dff,obj.metadata.ind_tr_1,...
                obj.Data_extract,obj.SavedCaTrials,frameNumTime,behEventAlign,...
                masklick,i_selectivity,obj.metadata.trial2include,obj.metadata.trial2exclude,...
                savename_figdff,title_fig);
        end
        
        function fig_rasterPSTH=mPlotActivityRasterPSTH(obj,activity_type,ind_ROI,...
                behEventAlignPool,masklickPool,i_selectivity,behEventSortPool)
            %   using other existed function
            obj.metadata.behEventAlign=behEventAlignPool;
            obj.metadata.masklick=masklickPool;
            obj.metadata.i_selectivity=i_selectivity;
            if isempty(ind_ROI) || min(ind_ROI)>obj.metadata.nROI || max(ind_ROI)<1
                iStart=1;
                iEnd=obj.metadata.nROI;
            else
                iStart=max(1,min(ind_ROI));
                iEnd=min(obj.metadata.nROI,max(ind_ROI));
            end
            for i_align=1:length(behEventAlignPool)
                subplot(1,length(behEventAlignPool),i_align);
                hold on;
                %decide some global variable
                behEventAlign=behEventAlignPool{i_align};%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
                behEventSort=behEventSortPool{i_align};% string can be in{'first lick','reward','go cue'};
                masklick=masklickPool{i_align};
                
                selectivitystr={'stimuli','sensory difficulty','sensory','choice'};%sensory means grouping difficulties;
                trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
                trialTypeStr=selectivitystr{i_selectivity};
                %plotting settings
                if strcmp(behEventAlign,'stim onset')
                    frameNumTime=[0.5,2];%from 5s before align point to 5s after align point
                elseif strcmp(behEventAlign,'delay onset')
                    frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
                else
                    frameNumTime=[2,3.5];%from 5s before align point to 5s after align point
                end
                mkdir([obj.path.root,filesep,obj.name]);
                if strcmp(activity_type,'dff')
                    savename_fig_activity=[obj.path.root,filesep,obj.name,filesep,obj.name,'-cat_f_inclusion'];
                    [figDff] = fLabelROIsDffIncludedTrialRange(obj.dff,savename_fig_activity,obj.SavedCaTrials,obj.metadata.trial2include);
                    for iROI=iStart:iEnd
                        savename_fig=[obj.path.root,filesep,obj.name,filesep,obj.name,'ROI-',num2str(iROI),'-alignTo',behEventAlign,'-sort',behEventSort,'-dff-rasterPSTH'];
                        disp(savename_fig);
                        title_fig=[obj.name,'ROI-',num2str(iROI)];
                        title_fig=strrep(title_fig,'_','\_');%下标变为转义字符的下划线
                        fig_rasterPSTH = fPlotDffRasterPSTH_ASession(obj.dff(iROI,:),obj.metadata.ind_tr_1,...
                            obj.Data_extract,obj.SavedCaTrials,frameNumTime,behEventAlign,...
                            masklick,i_selectivity,behEventSort,obj.metadata.indTrial2use,...
                            savename_fig,title_fig,trialTypeStr);
                    end
                elseif strcmp(activity_type,'spkr')
                    savename_fig_activity=[obj.path.root,filesep,obj.name,filesep,obj.name,'-spkr_inclusion'];
                    [figDff] = fLabelROIsDffIncludedTrialRange(obj.spkr,savename_fig_activity,obj.SavedCaTrials,obj.metadata.trial2include);
                    for iROI=iStart:iEnd
                        savename_fig=[obj.path.root,filesep,obj.name,filesep,obj.name,'ROI-',num2str(iROI),'-alignTo',behEventAlign,'-sort',behEventSort,'-spkr-rasterPSTH'];
                        disp(savename_fig);
                        title_fig=[obj.name,'ROI-',num2str(iROI)];
                        title_fig=strrep(title_fig,'_','\_');%下标变为转义字符的下划线
                        fig_rasterPSTH = fPlotRasterPSTH_ASession(obj.spkr(iROI,:),'spkr',obj.metadata.ind_tr_1,...
                            obj.Data_extract,obj.SavedCaTrials,frameNumTime,behEventAlign,...
                            masklick,i_selectivity,behEventSort,obj.metadata.indTrial2use,...
                            savename_fig,title_fig,trialTypeStr);
                    end
                elseif  strcmp(activity_type,'zscored_spkr')
                    savename_fig_activity=[obj.path.root,filesep,obj.name,filesep,obj.name,'-zscored_spkr_inclusion'];
                    [figDff] = fLabelROIsDffIncludedTrialRange(obj.dff,savename_fig_activity,obj.SavedCaTrials,obj.metadata.trial2include);
                    for iROI=iStart:iEnd
                        savename_fig=[obj.path.root,filesep,obj.name,filesep,obj.name,'ROI-',num2str(iROI),'-alignTo',behEventAlign,'-sort',behEventSort,'-zscored_spkr-rasterPSTH'];
                        disp(savename_fig);
                        title_fig=[obj.name,'ROI-',num2str(iROI)];
                        title_fig=strrep(title_fig,'_','\_');%下标变为转义字符的下划线
                        fig_rasterPSTH = fPlotRasterPSTH_ASession(obj.zscored_spkr(iROI,:),'zscored_spkr',obj.metadata.ind_tr_1,...
                            obj.Data_extract,obj.SavedCaTrials,frameNumTime,behEventAlign,...
                            masklick,i_selectivity,behEventSort,obj.metadata.indTrial2use,...
                            savename_fig,title_fig,trialTypeStr);
                    end
                end
                    
            end
        end
        
        function [fig_A,fig_B]=mInspectActivity(obj,iROI,activity_form,behEventAlign,masklick,behEventSort,chosen_result,chosen_stim,n_trial_show)
            selectivitystr={'stimuli','sensory difficulty','choice','sensory'};%sensory means grouping difficulties;
            trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
            i_selectivity=3;
            trialTypeStr=selectivitystr{i_selectivity};
            %plotting settings
            if strcmp(behEventAlign,'stim onset') && strcmp(masklick,'yes')
                frameNumTime=[0.5,2];%from 5s before align point to 5s after align point
            elseif strcmp(behEventAlign,'delay onset') && strcmp(masklick,'yes')
                frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
            else
                frameNumTime=[2,3.5];%from 5s before align point to 5s after align point
            end
            if strcmp(activity_form,'dff')
                title_fig=[obj.name,'ROI-',num2str(iROI),'-dff'];
                [fig_A,fig_B] = fActivityInspect(obj.dff(iROI,:),'dff',obj.metadata.ind_tr_1,...
                    obj.Data_extract,obj.SavedCaTrials,frameNumTime,behEventAlign,...
                    masklick,i_selectivity,behEventSort,obj.metadata.trial2include,obj.metadata.trial2exclude,...
                    title_fig,trialTypeStr,chosen_result,chosen_stim,n_trial_show);
            elseif strcmp(activity_form,'spkr')
                title_fig=[obj.name,'ROI-',num2str(iROI),'-spkr'];
                [fig_A,fig_B] = fActivityInspect(obj.spkr(iROI,:),'spkr',obj.metadata.ind_tr_1,...
                    obj.Data_extract,obj.SavedCaTrials,frameNumTime,behEventAlign,...
                    masklick,i_selectivity,behEventSort,obj.metadata.trial2include,obj.metadata.trial2exclude,...
                    title_fig,trialTypeStr,chosen_result,chosen_stim,n_trial_show);
            end
        end
        
        function [meanActivityByTrialType,cellActivityByTrialType, Tout]=mMeanPSTHbyROI(obj,activity_type,ind_ROI,behEventAlign,masklick,i_selectivity)
            %meanActivityByTrialType-mean activities of each ROI,each
            %condition one vector
            %cellActivityByTrialType-cell of activities of each ROI, each
            %condition one cell
            %Tout- table combine previous mean and cell output with
            %additional information about ROI identity
            switch activity_type
                case 'dff'
                    activities_data=obj.dff;
                case 'zscored_dff'
                    activities_data=obj.zscored_dff;
                case 'spkr'
                    activities_data=obj.spkr;
                case 'zscored_spkr'
                    activities_data=obj.zscored_spkr;
            end

            %plotting settings
            if strcmp(behEventAlign,'stim onset')
                frameNumTime=[0.5,2];%from 5s before align point to 5s after align point
            elseif strcmp(behEventAlign,'delay onset')
                frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
            else
                frameNumTime=[1,2];%from 1s before align point to 2s after align point
            end
            if isempty(ind_ROI) || min(ind_ROI)>obj.metadata.nROI || max(ind_ROI)<1
                iStart=1;
                iEnd=obj.metadata.nROI;
            else
                iStart=max(1,min(ind_ROI));
                iEnd=min(obj.metadata.nROI,max(ind_ROI));
            end

            varROI=iStart:iEnd;
            varROI=reshape(varROI,[],1);
            nROI=iEnd-iStart+1;
            [varMean,varCell]=deal(cell(iEnd-iStart+1,1));
            [varanimal{1:nROI}]=deal(obj.metadata.animal);
            varanimal=reshape(varanimal,[],1);
            [vardate{1:nROI}]=deal(obj.metadata.date);
            vardate=reshape(vardate,[],1);
            for iROI=iStart:iEnd
                [meanActivityByTrialType_ROI,cellActivityByTrialType_ROI] = fMeanPSTHbyROI(activities_data(iROI,:),...
                    obj.metadata.ind_tr_1,obj.Data_extract,obj.SavedCaTrials,...
                    frameNumTime,behEventAlign,masklick,i_selectivity,...
                    obj.metadata.trial2include,obj.metadata.trial2exclude);
                varMean{iROI,1}=meanActivityByTrialType_ROI;
                varCell{iROI,1}=cellActivityByTrialType_ROI;
                if exist('meanActivityByTrialType','var')
                    meanActivityByTrialType=cellfun(@(x,y) vertcat(x,y),meanActivityByTrialType,meanActivityByTrialType_ROI,'UniformOutput',false);
                    cellActivityByTrialType=cellfun(@(x,y) vertcat(x,y),cellActivityByTrialType,cellActivityByTrialType_ROI,'UniformOutput',false);
                else
                    meanActivityByTrialType=meanActivityByTrialType_ROI;
                    cellActivityByTrialType=cellActivityByTrialType_ROI;
                end
            end
            Tout=table(varanimal,vardate,varROI,varMean,varCell,...
                'VariableNames',{'animal','date','nROI','meanActivityByTrialType','cellActivityByTrialType'});
        end
        
        function [indBaseline] = mIndexBaselineAroundMode(obj,signals,varargin)
            file_name=[signals,'_IndBaselineAroundMode.mat'];
            if exist(strcat(obj.path.root,filesep,file_name),'file')
                load(strcat(obj.path.root,filesep,file_name));
            else
                if isempty(varargin)
                    sigThreshSD=2;
                else
                    if mod(length(varargin),2)==0
                        for i=1:2:length(varargin)
                            switch varargin{i}
                                case 'sigThreshSTD'
                                    sigThreshSD=varargin{i+1};
                            end
                        end
                    else
                        warning('odd varargin input argument');
                    end
                end
                % calculate a index of baseline around signal mode
                if strcmp(signals,'spkr')
                    %indBaseline = fIndexBaselineAroundMode(obj.spkr,'sigThreshSTD',sigThreshSD);% 2 times std below as baseline
                    indBaseline = true(size(obj.spkr)); %seems that this is currently the best practise
                elseif strcmp(signals,'dff')
                    indBaseline = fIndexBaselineAroundMode(obj.dff,'sigThreshSTD',sigThreshSD);% 2 times std below as baseline
                end
                save(strcat(obj.path.root,filesep,file_name),'indBaseline');
            end
        end
        
        function [frac_event,frac_event_start]=mSignificant(obj,signals,criteria,threshold_dur,sample_method,flag_plot)
            % signals, usually obj.spkr
            %   signals - m-by-n matrix, m ROIs and n time points
            %   criteria- char, eg. 2STD,
            %   threshold_dur- number of time (unit) needed to be included
            %   as significant, eg. 3
            
            indBaseline=obj.mIndexBaselineAroundMode(signals);
            % calculate a significant matrix
            if strcmp(signals,'spkr')     
                [significant_matrix,start_sig_mat] = fSignificant(obj.spkr,criteria,threshold_dur,'BaselineIndex',indBaseline);%continues 3 frames above 2 times std
            elseif strcmp(signals,'dff')
                [significant_matrix,start_sig_mat] = fSignificant(obj.dff,criteria,threshold_dur,'BaselineIndex',indBaseline);%continues 3 frames above 2 times std
            end
            % compare duration/start probability of significant signals happen in different epochs
            indITI=[];
            indSound=[];
            indDelay=[];
            indLick=[];
            behEventFrameIndex=obj.behaviorEvent.behEventFrameIndex;%for simplicity
            if strcmp(sample_method,'random')% method 1, randomly sample time points
                for indTrial=1:obj.metadata.ntr-1
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
            elseif strcmp(sample_method,'all') % method 2, sample all the points and get a mean
                for indTrial=1:obj.metadata.ntr-1
                    indITI_current=behEventFrameIndex.start(indTrial):behEventFrameIndex.stimOnset(indTrial);
                    indITI=[indITI,indITI_current];
                    indSound_current=behEventFrameIndex.stimOnset(indTrial):behEventFrameIndex.stimOffset(indTrial);
                    indSound=[indSound,indSound_current];
                    indDelay_current=behEventFrameIndex.stimOffset(indTrial):behEventFrameIndex.go(indTrial);
                    indDelay=[indDelay,indDelay_current];
                    indLick_current=behEventFrameIndex.go(indTrial):30+behEventFrameIndex.go(indTrial);
                    indLick=[indLick,indLick_current];
                end
            end
            
            % calculate the probability of significant event or their start
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
            % summarize across sessions and compare the proportion of significant activities happened
            if strcmp(flag_plot,'show_case')
                fig_SigProb=figure;
                set(gcf,'Position',[100,100,600,300]);
                ax1=subplot(1,2,1);
                ax1=fScatterStat(ax1, frac_event,'Probability of significant activity');
                ax2=subplot(1,2,2);
                ax2=fScatterStat(ax2, frac_event_start,'Probability of significant activity start');
                saveas(fig_SigProb,[obj.path.root,filesep,obj.name,filesep,'probability_significance.jpg'],'jpg');
            end
        end
        
        function [ind_NSDelayMovement,DelayMovement_criteria]=mGetIndNSDdelayMovement(obj,baseline,varargin)
            bodyparts='Tongue';
            coordinates='x';
            treshold4likelihood=0.1;
            str_nFrames='500ms';
            trialTypeStr='cor and err';
%             baseline='same';%same as bodyparts and coordinates
            if ~isempty(varargin)
                nvar=length(varargin);
                if mod(nvar/2)~=0
                    warning('invalid input number for member function: mGetIndNSDdelayMovement, use default settings');
                else
                    for i=1:2:nvar
                        switch varargin{i}
                            case 'bodyparts'
                                bodyparts=varargin{i+1};
                            case 'coordinates'
                                coordinates=varargin{i+1};
                            case 'treshold4likelihood'
                                treshold4likelihood=varargin{i+1};
                            case 'str_nFrames'
                                str_nFrames=varargin{i+1};
                            case 'trialTypeStr'
                                trialTypeStr=varargin{i+1};
                        end
                    end
                end
            end
            switch trialTypeStr
                case 'cor'
                    combineCorErr='divideCorErr';%and only use correct trial
                case 'cor and err'
                    combineCorErr='combineCorErr';%{'combineCorErr','divideCorErr'}
            end
            [trialType,~,~] = fGetTrialType( obj.Data_extract,[],3,'matrix','left',combineCorErr); %default, trial type by choice
            ind_trial=logical(reshape(sum(trialType(1:end-2,:,:),2),[],1));%only correct/error trials
            if isempty(obj.DLC.dcnum) || isempty(baseline)%no video and DLC related data
                DelayMovement_criteria='All_trial_include';
                ind_NSDelayMovement=true(obj.metadata.ntr,1);
            else
                [behEventFrameIndex_dlc,~] = fGetBehEventTime( obj.Data_extract,...
                obj.metadata.frameTrialStartVideo, obj.metadata.video_frameRate );%get behavior event time
                bodyco = obj.mGetDLCcoordinate(bodyparts,coordinates,treshold4likelihood);
                [T_DLCbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex_dlc,...
                    bodyco,1000/obj.metadata.video_frameRate,str_nFrames);%DLC
                if strcmp(baseline,'same')
                    bsdata=T_DLCbyEpoch.sound(ind_trial);
                    DelayMovement_criteria='Sound_as_baseline';
                elseif strcmp(baseline,'LickPort')
                    bodycoLickport = obj.mGetDLCcoordinate('LeftLickPort',coordinates,treshold4likelihood);
                    [T_DLCbyEpoch_LP,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex_dlc,...
                        bodycoLickport,1000/obj.metadata.video_frameRate,str_nFrames);%DLC
                    bsdata=T_DLCbyEpoch_LP.delay(ind_trial);
                    DelayMovement_criteria='LickPort_as_baseline';
                end
                baseline=[mean(bsdata)-2*std(bsdata),mean(bsdata)+2*std(bsdata)];
                xdlc=T_DLCbyEpoch.delay;
                ind_sig= logical((xdlc<baseline(1))+(xdlc>baseline(2)));
                ind_NSDelayMovement=(~ind_sig);
            end
        end
        
        function [nfig,struct_rho,struct_p] = mPlotDLCvsDff(obj,epoch,...
                str_nFrames,trialTypeStr,ind_NSDelayMovement,...
                bodyparts,coordinates,treshold4likelihood,savepath,savename)
            %plot correlation between DLC movement and dff amplitude
            switch trialTypeStr
                case 'cor'
                    combineCorErr='divideCorErr';%and only use correct trial
                case 'cor and err'
                    combineCorErr='combineCorErr';%{'combineCorErr','divideCorErr'}
            end
            [trialType,~,~] = fGetTrialType( obj.Data_extract,[],3,'matrix','left',combineCorErr); %default, trial type by choice
            ind_trial=logical(reshape(sum(trialType(1,:,:),2),[],1));%only correct trials
            [behEventFrameIndex_dff,~] = fGetBehEventTime( obj.Data_extract, ...
                obj.metadata.ind_1stFrame, obj.SavedCaTrials.FrameTime );%get behavior event time
            [behEventFrameIndex_dlc,~] = fGetBehEventTime( obj.Data_extract,...
                obj.metadata.frameTrialStartVideo, obj.metadata.video_frameRate );%get behavior event time
            bodyco = obj.mGetDLCcoordinate(bodyparts,coordinates,treshold4likelihood);
            [T_DLCbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex_dlc,...
                bodyco,1000/obj.metadata.video_frameRate,str_nFrames);%DLC
            baseline=[mean(T_DLCbyEpoch.sound(ind_trial))-2*std(T_DLCbyEpoch.sound(ind_trial)),mean(T_DLCbyEpoch.sound(ind_trial))+2*std(T_DLCbyEpoch.sound(ind_trial))];
            ncol=6;
            nrow=10;
            nfig=ceil(ceil(obj.metadata.nROI/ncol)/nrow);
            ifig=0;
            for roiNo = 1:obj.metadata.nROI
                if mod(roiNo,nrow*ncol)==1
                    ifig=ifig+1;
                    figure;
                    set(gcf,'Position',[10,10,600,1000]);
                end
                [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex_dff,...
                    obj.dff(roiNo,:),obj.metadata.frT,str_nFrames);%dff
                
                figure(ifig);
                ind_subplot=mod(roiNo,nrow*ncol);
                if ind_subplot==0
                    ind_subplot=nrow*ncol;
                end
                subplot(nrow,ncol,ind_subplot);
                
                switch epoch
                    case 'stim'
                        xdlc=T_DLCbyEpoch.sound(ind_trial);
                        ydff=T_SigbyEpoch.sound(ind_trial);
                    case 'delay'
                        xdlc=T_DLCbyEpoch.delay(ind_trial);
                        ydff=T_SigbyEpoch.delay(ind_trial);
                    case 'response'
                        xdlc=T_DLCbyEpoch.response(ind_trial);
                        ydff=T_SigbyEpoch.response(ind_trial);
                    case 'lick'
                        xdlc=T_DLCbyEpoch.lick(ind_trial);
                        ydff=T_SigbyEpoch.lick(ind_trial);
                end
                ind_sig= (~ind_NSDelayMovement);
                scatter(xdlc(~ind_sig),ydff(~ind_sig),10,'k');hold on;
                [r_ns,p_ns]=corrcoef(xdlc(~ind_sig),ydff(~ind_sig),'Rows','pairwise');%both x,y not nan
                scatter(xdlc(ind_sig),ydff(ind_sig),10,'r');
                [r_sig,p_sig]=corrcoef(xdlc(ind_sig),ydff(ind_sig),'Rows','pairwise');%both x,y not nan
                if size(r_ns,1) ==1 %only one datapoint or no data
                    struct_rho(roiNo).ns=nan;
                    struct_p(roiNo).ns=nan;
                else
                    struct_rho(roiNo).ns=r_ns(2,1);
                    struct_p(roiNo).ns=p_ns(2,1);
                end
                if size(r_sig,1)==1
                    struct_rho(roiNo).sig=nan;
                    struct_p(roiNo).sig=nan;
                else
                    struct_rho(roiNo).sig=r_sig(2,1);
                    struct_p(roiNo).sig=p_sig(2,1);
                end
%                 text(0.1,0.9,['rho=',num2str(struct_rho(roiNo).ns),',p=',num2str(struct_p(roiNo).ns)],'Units','Normalized');
%                 text(0.1,1,['rho=',num2str(struct_rho(roiNo).sig),',p=',num2str(struct_p(roiNo).sig)],'Units','Normalized','Color','r');
                text(0.1,0.9,plabelsymbol(struct_p(roiNo).ns),'Units','Normalized');
                text(0.1,1,plabelsymbol(struct_p(roiNo).sig),'Units','Normalized','Color','r');
                
                if mod(roiNo,nrow*ncol)==1
                    ylabel(['\DeltaF/F during ',epoch]);
                end
                if mod(roiNo,nrow*ncol)==(nrow*ncol-ncol+1) || roiNo==obj.metadata.nROI
                    xlabel([bodyparts,'-',coordinates,' (pixel to ipsi)']);
                end
            end
            for ifig=1:nfig
                saveas(figure(ifig),[savepath,filesep,savename,'-',num2str(ifig),'.png'],'png');
                saveas(figure(ifig),[savepath,filesep,savename,'-',num2str(ifig),'.fig'],'fig');
            end
            close all;
        end
        
        function figDLC_PSTH=mPlotDLC_PSTH(obj,bodyparts,coordinates,ind_NSDelayMovement,...
                selectivity,pSigTtest,treshold4likelihood,behEventAlign,...
                yrange,titlestr,nonOverlappingBin,savepath,savename)
            %ref fCorrelationBehaviorDistributedSession
            selectivitystr={'stimuli','difficulty','sensory','choice'};%sensory means grouping difficulties;
            trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
            temp=cellfun(@(x) strcmp(x,selectivity),  selectivitystr);
            i_selectivity=find(temp);%*********variable**************
            outcomeType='combineCorErr';%'combineCorErr','divideCorErr',
            binsize=1;%bin size for moving ttest
            binstep=1;%bin step for moving ttest
            % pSigTtest=0.01;%significance level for moving ttest
            fr=obj.metadata.video_frameRate/nonOverlappingBin;
            nbodyparts=size(bodyparts,1)*size(bodyparts,2);
            frameTrialStartVideo=floor(obj.metadata.frameTrialStartVideo/nonOverlappingBin);
            % get behavior event from .beh
            ind_tr_1=1;%using data from trial 1
            [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( obj.Data_extract, frameTrialStartVideo, 1000/fr ,ind_tr_1);%get behavior event time
            [trialType,rule,trialTypeStr] = fGetTrialType( obj.Data_extract,[],trialTypeVar(i_selectivity),'matrix','left',outcomeType);%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
            bodycoSelectedTrial=cell(1,nbodyparts);%each cell represent one session mean
            bodycoDiffSelectedTrial=cell(1,nbodyparts);%each cell represent one session mean
            for iBodyPart=1:nbodyparts
                %align to behavior event
                behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
                if strcmp(behEventAlign{iBodyPart},'stim onset') || strcmp(behEventAlign{iBodyPart},'delay onset')
                    masklick='yes';
                else
                    masklick='no';
                end
                if strcmp(behEventAlign{iBodyPart},'stim onset')
                    frameNumTime=[0.5,2];%plotting settings, from 0.5s before align point to 2s after align point
                    tempXticklabel=[0,1,2];
                elseif strcmp(behEventAlign{iBodyPart},'delay onset')
                    frameNumTime=[1,1.4];
                    tempXticklabel=[-.5,0,.5,1];
                else
                    frameNumTime=[1,2];
                    tempXticklabel=[-1,0,1,2];
                end
                frameNum=double(floor(frameNumTime*fr));
                bodyco = obj.mGetDLCcoordinate(bodyparts{iBodyPart},coordinates{iBodyPart},treshold4likelihood);
                % aligned the coordinates change of body parts to behavior events
                if strcmp(behEventAlign{iBodyPart},'delay onset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
                    disp(bodyparts{iBodyPart});
                    [ bodyco_aligned, ~,~ ] = fAlignDelaySigal( bodyco, behEventFrameIndex,  frameNum ,'raw');
                else
                    [ bodyco_aligned, ~, ~ ] = fAlignSigBehEvent( bodyco, behEventFrameIndex,lickingFrameIndex,behEventAlign{iBodyPart},frameNum );%decide which behavior event to align to
                end
                nfigcol=1;%4 column(correct/error/miss/violation),here only use correct trials
                bodycoSelectedTrial{iBodyPart}=cell(size(trialType,2),nfigcol);
                bodycoDiffSelectedTrial{iBodyPart}=cell(1,nfigcol);%diff between two groups mean
                for nStim=1:size(trialType,2) %for each stimulus
                    for  nResult=1:nfigcol
                        selectedTrialInd=trialType(nResult,nStim,:);
                        selectedTrialInd=logical(squeeze(selectedTrialInd))';
                        selectedTrialInd=reshape(selectedTrialInd,[],1);
                        ind_NSDelayMovement=reshape(ind_NSDelayMovement,[],1);
                        selectedTrialInd_end=logical(selectedTrialInd.*ind_NSDelayMovement);
                        bodycoSelectedTrial{iBodyPart}{nStim,nResult}=bodyco_aligned(selectedTrialInd_end,:);
                    end
                end
                for  nResult=1:nfigcol
                    tempbodycoSelectedTrial1=cell2mat(bodycoSelectedTrial{iBodyPart}(1:size(trialType,2)/2,nResult));
                    tempbodycoSelectedTrial2=cell2mat(bodycoSelectedTrial{iBodyPart}(size(trialType,2)/2+1:end,nResult));
                    bodycoDiffSelectedTrial{iBodyPart}{1,nResult}=nanmean(tempbodycoSelectedTrial1,1)-nanmean(tempbodycoSelectedTrial2,1);
                end
            end
            %plot
            figDLC_PSTH=figure;%plot mean trace
            % set(gcf, 'position', [0 0 250*nfigcol 200*nbodyparts]);
            set(gcf, 'PaperPosition', [0 0 2*size(bodyparts,1) 2*size(bodyparts,2)]);
            for iBodyPart=1:nbodyparts
                % plot coordinates change, ref plotDffPSTH
                
                if size(trialType,2)==6
                    color_mean_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};
                elseif size(trialType,2)==2
                    color_mean_trace={[0 0 1],[1 0 0]};
                    color_cases={[0.5 0.5 1],[1 0.5 0.5]};
                end
                if strcmp(rule,'low click rate-right')
                    color_mean_trace=fliplr(color_mean_trace);
                end
                matCasesStimResult=cell(1,size(trialType,2));
                for  nStim=1:size(trialType,2) %for each stimulus
                    bodycoSelectedTrialCase=bodycoSelectedTrial{iBodyPart}{nStim,nResult};
                    matCasesStimResult{nStim}=bodycoSelectedTrialCase;
                end
                n_shuffle=1000;    %%%%%%%%%%%%%%%%%%%%%%%%%%%change for shuffle times
                %due to sample rate difference, matrix in cell may have different
                %length, which need fix by inter1
                temp1=fInter1(matCasesStimResult(:,1));
                temp2=fInter1(matCasesStimResult(:,2));
                category1=cell2mat(temp1);%note here that only useful for 2 stimuli, if multiple stimuli, then this cause error
                category2=cell2mat(temp2);
                activity=cat(1,category1,category2);
                label = cat(1,ones(size(category1,1),1),2*ones(size(category2,1),1));
                poslabel=2;
                fileAUC=strcat(obj.name,'-',bodyparts{iBodyPart},'-',coordinates{iBodyPart},'-',savename,'-align',behEventAlign{iBodyPart},'-binsize',num2str(binsize),'-nshuffle',num2str(n_shuffle),'-MovingAUC.mat');
                if exist([obj.path.root,filesep,'AUC',filesep,fileAUC],'file')
                    load([obj.path.root,filesep,'AUC',filesep,fileAUC]);
                else
                    [auc,pAUC,activitySmoothed] = fMovingAUC(label,activity,poslabel,n_shuffle,binsize,binstep);
                end
                if ~exist([obj.path.root,filesep,'AUC'],'dir')
                    mkdir([obj.path.root,filesep,'AUC']);
                end
                save([obj.path.root,filesep,'AUC',filesep,fileAUC],'auc','pAUC','activitySmoothed');
                
                %plot mean, ci, individual stimuli etc after AUC analysis.
                neuralActivityMean=cell(1,size(trialType,2));
                for  nStim=1:size(trialType,2) %for each stimulus; in fact only support 2 category
                    figure(figDLC_PSTH);%save mean trace
                    subplot(size(bodyparts,2),size(bodyparts,1),iBodyPart);
                    plot(1:size(matCasesStimResult{nStim},2),matCasesStimResult{nStim},'Color',color_cases{nStim},'linewidth',0.5);
                    hold on;
                    [neuralActivityMean{nStim},neuralActivityCI]=fMean_SE(activitySmoothed{nStim});
                    activity_Std=nanstd(activitySmoothed{nStim},1,1);
                    neuralActivityCI(1,:)=nanmean(activitySmoothed{nStim},1)+activity_Std;
                    neuralActivityCI(2,:)=nanmean(activitySmoothed{nStim},1)-activity_Std;
                    %plot CI
                    if size(trialType,2)==2%if only two stimuli, plot CI, otherwise, do not plot for the reason of simplicity
                        xpatch=[1:size(neuralActivityCI,2), fliplr(1:size(neuralActivityCI,2))];
                        ypatch=[neuralActivityCI(1,:),fliplr(neuralActivityCI(2,:))];
                        p=patch(xpatch,ypatch,color_mean_trace{nStim});%plot confidence interval
                        p.FaceAlpha=0.1;
                        p.EdgeColor=color_mean_trace{nStim};%'none';
                        hold on;
                    end
                    title(titlestr{iBodyPart});
                end
                %plot mean trace
                for  nStim=1:size(trialType,2) %for each stimulus; in fact only support 2 category          
                    curve_meanTrace(nStim)=plot(1:size(neuralActivityMean{nStim},2),neuralActivityMean{nStim},'Color',color_mean_trace{nStim},'linewidth',2);
                end
                
                if ~exist('yrange','var') || isempty(yrange)
                    y_lim=get(gca,'Ylim');
                else
                    y_lim=yrange{iBodyPart};
                end
                %       %plot significant bar (indicating time point of significant)
                indSig=logical((pAUC<pSigTtest/2)+(pAUC>1-pSigTtest/2));
                %         indSig=(pTtest<pSigTtest);
                ySig=ones(size(pAUC))*y_lim(2)*0.9;
                xSig=(1:length(pAUC));
                xSig(~indSig)=nan;
                ySig(~indSig)=nan;
                xSig=fRuleOutOccasional(xSig,1);
                ySig=fRuleOutOccasional(ySig,1);
                plot(xSig,ySig,'k-','LineWidth',1);
                %label the figure
                if strcmp(coordinates{iBodyPart},'likelihood')
                    yrange={[ 0 , 1.1 ]};
                end
                if ~exist('yrange','var')|| isempty(yrange)
                    y_lim=get(gca,'Ylim');
                else
                    y_lim=yrange{iBodyPart};
                end
                xlim([1,sum(frameNum)+1]);
                plot([frameNum(1)+1,frameNum(1)+1],[y_lim(1),y_lim(2)],'k-');%align to a behavior event
                hold on;
                if strcmp(behEventAlign{iBodyPart},'stim onset')
                    plot((frameNum(1)+1+0.5*fr)*ones(1,2),[y_lim(1),y_lim(2)],'k-');%plot delay onset
                elseif strcmp(behEventAlign{iBodyPart},'delay onset')
                    plot((frameNum(1)+1-0.5*fr)*ones(1,2),[y_lim(1),y_lim(2)],'k-');%plot stim onset
                end
                if ceil(iBodyPart/size(bodyparts,1))==size(bodyparts,2)
                    xlabel(['Time (s) from ',behEventAlign{iBodyPart}]);
                end
                if mod(iBodyPart,size(bodyparts,1))==1 || size(bodyparts,1)==1
                    if strcmp(coordinates{iBodyPart},'likelihood')
                        ylabel('likelihood');
                    else
                        ylabel('{\it\Delta}pixel');
                    end
                end
                set(gca,'xtick',frameNum(1)+fr*tempXticklabel+1,'xticklabel',tempXticklabel);
                %         set(gca,'xtick',[floor(fr*(frameNumTime(1)-floor(frameNumTime(1)))):floor(fr):size(neuralActivity,2)],'xticklabel',tempXticklabel);
                set(gca,'FontName','Arial','FontSize',12);
                set(gca,'Ylim',y_lim);
                box off;
            end
            subplot(size(bodyparts,2),size(bodyparts,1),1);
            if contains(trialTypeStr,'stimuli')
                h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','','contra hard','','contra easy'},'Location','best');
                %legend(curve_meanTrace(:),num2str(Data_extract.Stimuli(:)),'Location','best');
            elseif contains(trialTypeStr,'difficulty')
                h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','contra hard','contra easy'},'Location','best');
            elseif contains(trialTypeStr,'first lick')
                h=legend(curve_meanTrace(:),{'ipsi lick first','contra lick first'},'Location','best');
            elseif contains(trialTypeStr,'sensory')
                h=legend(curve_meanTrace(:),{'low click rate','high click rate'},'Location','best');
            end
            set(h,'box','off');
            saveas(figDLC_PSTH,[obj.path.root,filesep,'AUC',filesep,obj.name,'-DLC_PSTH',savename,'-algin',behEventAlign{iBodyPart},'-sort',behEventSort,'mean_trace.pdf'],'pdf');
            saveas(figDLC_PSTH,[obj.path.root,filesep,'AUC',filesep,obj.name,'-DLC_PSTH',savename,'-algin',behEventAlign{iBodyPart},'-sort',behEventSort,'mean_trace.fig'],'fig');
            saveas(figDLC_PSTH,[savepath,filesep,obj.name,'-DLC_PSTH',savename,'-algin',behEventAlign{iBodyPart},'-sort',behEventSort,'mean_trace.png'],'png');
            close all;
        end
        
        function [TAUC, Tmean,obj] = mGetEpochAUC(obj,celltype,ROItype,trialTypeStr,...
                AUCtype,AUCCorrectedMethod)
            date_str=datestr(obj.metadata.datetime,'yyyy-mm-dd');
            date_str=strrep(date_str,'-0','-');
            selectivitystr={'stimuli','sensory difficulty','sensory','choice'};%sensory means grouping difficulties;
            trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
            i_selectivity=cellfun(@(x) strcmp(x,AUCtype),selectivitystr);%*********variable**************
            i_selectivity=find(i_selectivity);
            str_nFrames='1s';%'500ms';%'1s'
            %naming the temp file by method of AUC calculation correction
            filepath_inSummary='E:\2P\summary\AUC\cases';%backup another copy in the summary folder
            fileName_datapars=[obj.metadata.animal,'-',date_str,'trialType',trialTypeStr];
            if strcmp(trialTypeStr,'cor') || strcmp(trialTypeStr,'err')
                combineCorErr='divideCorErr';%and only use correct trial
                corErrTrialNumber='raw';
                fileName_processpars=[AUCtype,'-trialNum',corErrTrialNumber,'-timeBin',str_nFrames];
            elseif strcmp(trialTypeStr,'cor and err')%used to test whether it is sensory or choice AUC,in order to be more fair, keep trial numbers balenced for correct and error trials
                combineCorErr='combineCorErr';%{'combineCorErr','divideCorErr'}
                switch AUCCorrectedMethod
                    case 'balencedCorErrTrialNum'
                        corErrTrialNumber='balence';
                        fileName_processpars=[AUCtype,'-trialNum',corErrTrialNumber,'-timeBin',str_nFrames];
                    case 'SensoryChoiceOrthogonalSubtraction'
                        corErrTrialNumber='raw';
                        if strcmp(AUCtype,'choice')|| strcmp(AUCtype,'sensory')
                            fileName_processpars=[AUCtype,'-AUCCorrectedMethod',AUCCorrectedMethod,'-timeBin',str_nFrames];
                        else
                            warning('error input combination of fGetEpochAUCtableASession function, SensoryChoiceOrthogonalSubtraction method');
                        end
                    otherwise
                        corErrTrialNumber='raw';
                        fileName_processpars=[AUCtype,'-trialNum',corErrTrialNumber,'-timeBin',str_nFrames];
                end
            end
            fileNameT=[obj.path.root,filesep,fileName_datapars,'-',fileName_processpars,'-EpochAUC.mat'];
            fileNameT_inSummary=[filepath_inSummary,filesep,fileName_datapars,'-',fileName_processpars,'-EpochAUC.mat'];
            
            if exist(fileNameT,'file')%CD058-2018-1-27trialTypecor and err-choice-AUCCorrectedMethodSensoryChoiceOrthogonalSubtraction-timeBin1s-EpochAUC
                load(fileNameT);
                n_ROI=size(TAUC,1);
                disp(['Table exist, use ',obj.metadata.animal,date_str,';nROI=',num2str(n_ROI)]);
                %delete duplicated variables
                standard_varNames={'animal','date','field','celltype','ROItype','nROI',...
                    'ipsi_performance', 'contra_performance','overall_performance',...
                    'ITI','sound','delay','response','lick','mid_delay','late_delay',...
                    'pITI','psound','pdelay','presponse','plick','pmid_delay','plate_delay','delayMovingAUC','pdelayMovingAUC'};
                T_names=TAUC.Properties.VariableNames;
                for i=1:length(T_names)
                    temp=cellfun(@(x) strcmp(x,T_names(i)), standard_varNames);
                    if sum(temp)==0
                        TAUC(:,T_names{i})=[];
                    end
                end
            else
                nROI=size(obj.SavedCaTrials.f_raw{1},1);
                [delayMovingAUC,pdelayMovingAUC]=deal(cell(nROI,1));
                [varanimal{1:nROI}]=deal(obj.metadata.animal);
                varanimal=reshape(varanimal,[],1);
                [vardate{1:nROI}]=deal(obj.metadata.date);
                vardate=reshape(vardate,[],1);
                [varfield{1:nROI}]=deal(obj.metadata.field);
                varfield=reshape(varfield,[],1);
                [varcelltype{1:nROI}]=deal(celltype);
                varcelltype=reshape(varcelltype,[],1);
                [varROItype{1:nROI}]=deal(ROItype);
                varROItype=reshape(varROItype,[],1);
                varROI=(1:size(obj.SavedCaTrials.f_raw{1},1));
                varROI=reshape(varROI,[],1);
                ind_tr_1=1;
                ntr=length(obj.SavedCaTrials.f_raw);
                frT = obj.SavedCaTrials.FrameTime;
                % align to behavior event
                nFrameEachTrial=cellfun(@(x) size(x,2),obj.SavedCaTrials.f_raw);
                ind_1stFrame=zeros(1,length(nFrameEachTrial));
                ind_1stFrame(1)=1;
                ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
                ind_1stFrame=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials
                frameNumTime=[1,1.5];%from 1s before delay onset to 1.5s after delay onset
                frameNum=double(round(frameNumTime*1000/frT));
                [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( obj.Data_extract, ind_1stFrame, obj.SavedCaTrials.FrameTime );%get behavior event time
                
                if strcmp(AUCtype,'choice')|| strcmp(AUCtype,'sensory')
                    [ITI, sound,delay, response, lick, pITI, psound, pdelay, presponse, plick]= deal(nan(nROI,1));
                    switch corErrTrialNumber
                        case 'balence'
                            [trialType,~,~] = fGetTrialType( obj.Data_extract,[],trialTypeVar(i_selectivity),'matrix','left','divideCorErr');
                            trialTypeIndCell=cell(2,2);%{ipsi cor, contra cor;ipsi err, contra err};
                            for i=1:2
                                for j=1:2
                                    trialTypeIndCell{i,j}=find(logical(reshape(trialType(i,j,:),[],1)));
                                end
                            end
                            temp=cellfun(@length,trialTypeIndCell);
                            trialNumEachType=floor(sum(temp,'all')/4);%keep total trial number stable and redistributed to each trial types
                            s = RandStream('mlfg6331_64');%for reproducibility
                            trialTypeIndCellBalenced=cellfun(@(x) reshape(datasample(s,x,trialNumEachType,'Replace',true),[],1), trialTypeIndCell,'UniformOutput',false);
                            trialTypeIndFinal=cell2mat(trialTypeIndCellBalenced);%1st col-cor; 2nd col-err
                            ind_trial=reshape(trialTypeIndFinal,[],1);
                            label_AUC=[ones(trialNumEachType*2,1);ones(trialNumEachType*2,1)*2];
                        case 'raw'
                            [trialType,~,~] = fGetTrialType( obj.Data_extract,[],trialTypeVar(i_selectivity),'matrix','left',combineCorErr);
                            label_choice = fTrialType2Label(trialType,2);
                            if contains(trialTypeStr,'cor') % 'cor'| 'cor and err'
                                ind_trial=logical(reshape(sum(trialType(1,:,:),2),[],1));%only correct trials
                            end
                            if strcmp(trialTypeStr,'err')
                                ind_trial=logical(reshape(sum(trialType(2,:,:),2),[],1));%only error trials
                            end
                            label_AUC=label_choice(ind_trial);
                            %used for orthogonal subtraction
                            [trialType_orthogonal,~,~] = fGetTrialType( obj.Data_extract,[],7-trialTypeVar(i_selectivity),'matrix','left',combineCorErr);
                            label_choice_orthogonal = fTrialType2Label(trialType_orthogonal,2);
                            label_AUC_orthogonal=label_choice_orthogonal(ind_trial);
                    end
                    for roiNo = 1:size(obj.SavedCaTrials.f_raw{1},1) %SavedCaTrials.nROIs may be not true
                        disp([obj.metadata.animal,obj.metadata.date,'rioNo',num2str(roiNo)]);
                        [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,obj.dff(roiNo,:),frT,str_nFrames);
                        if strcmp(AUCCorrectedMethod,'SensoryChoiceOrthogonalSubtraction')
                            T_SigbyEpoch=fOrthogonalSubtraction(T_SigbyEpoch,ind_trial,label_AUC,label_AUC_orthogonal);
                        else
                            T_SigbyEpoch=T_SigbyEpoch(ind_trial,:);
                        end
                        poslabel=2;
                        nshuffle=1000;%%%%check
                        [ITI(roiNo),pITI(roiNo)]=fAUC(label_AUC,T_SigbyEpoch.ITI,poslabel,nshuffle);
                        [sound(roiNo),psound(roiNo)]=fAUC(label_AUC,T_SigbyEpoch.sound,poslabel,nshuffle);
                        [delay(roiNo),pdelay(roiNo)]=fAUC(label_AUC,T_SigbyEpoch.delay,poslabel,nshuffle);
                        [response(roiNo),presponse(roiNo)]=fAUC(label_AUC,T_SigbyEpoch.response,poslabel,nshuffle);
                        [lick(roiNo),plick(roiNo)]=fAUC(label_AUC,T_SigbyEpoch.lick,poslabel,nshuffle);
                        %calculate moving AUC
                        [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( obj.dff(roiNo,:), behEventFrameIndex,  frameNum );
                        dff_late_delay=nanmean(dff_aligned(:,frameNum(1)+round(1*1000/frT):frameNum(1)+round(1.5*1000/frT)),2);
                        dff_mid_delay=nanmean(dff_aligned(:,frameNum(1)+round(0.3*1000/frT):frameNum(1)+round(1*1000/frT)),2);
                        if strcmp(AUCCorrectedMethod,'SensoryChoiceOrthogonalSubtraction')
                            dff_aligned_ortho_corrected=fOrthogonalSubtraction(dff_aligned,ind_trial,label_AUC,label_AUC_orthogonal);
                            dff_mid_delay4auc=fOrthogonalSubtraction(dff_mid_delay,ind_trial,label_AUC,label_AUC_orthogonal);
                            dff_late_delay4auc=fOrthogonalSubtraction(dff_late_delay,ind_trial,label_AUC,label_AUC_orthogonal);
                        else
                            dff_aligned=dff_aligned(ind_trial,:);
                            dff_mid_delay4auc=dff_mid_delay(ind_trial,:);
                            dff_late_delay4auc=dff_late_delay(ind_trial,:);
                        end
                        [mid_delay(roiNo),pmid_delay(roiNo)]=fAUC(label_AUC,dff_mid_delay4auc,poslabel,nshuffle);
                        [late_delay(roiNo),plate_delay(roiNo)]=fAUC(label_AUC,dff_late_delay4auc,poslabel,nshuffle);
                        
                        binsize=1;
                        binstep=1;
                        for nResult=1:size(trialType,1)-2
                            if strcmp(trialTypeStr,'cor and err')%combine cor and err together
                                indTrial=ind_trial;
                            else %calculate AUC for cor and err respectively
                                indTrial=trialType(nResult,:,:);
                                indTrial=sum(squeeze(indTrial),1);
                                indTrial=logical(squeeze(indTrial));
                            end
                            if nResult==1
                                if strcmp(trialTypeStr,'cor and err')
                                    if strcmp(AUCCorrectedMethod,'SensoryChoiceOrthogonalSubtraction')
                                        [delayMovingAUC{roiNo}.do,pdelayMovingAUC{roiNo}.do] = fMovingAUC(label_AUC,dff_aligned_ortho_corrected,2,nshuffle,binsize,binstep);
                                    else
                                        [delayMovingAUC{roiNo}.do,pdelayMovingAUC{roiNo}.do] = fMovingAUC(label_AUC,dff_aligned,2,nshuffle,binsize,binstep);
                                    end
                                elseif strcmp(trialTypeStr,'cor')
                                    [delayMovingAUC{roiNo}.cor,pdelayMovingAUC{roiNo}.cor] = fMovingAUC(label_AUC,dff_aligned,2,nshuffle,binsize,binstep);
                                end
                            elseif nResult==2 && strcmp(trialTypeStr,'err')
                                [delayMovingAUC{roiNo}.err,pdelayMovingAUC{roiNo}.err] = fMovingAUC(label_AUC,dff_aligned,2,nshuffle,binsize,binstep);
                            end
                        end
                    end
                    
                elseif strcmp(AUCtype,'stimuli') %here, compare auc of cor/err for each stimuli
                    [trialType,~,~] = fGetTrialType( obj.Data_extract,[],trialTypeVar(i_selectivity),'matrix','left','divideCorErr');
                    [ITI, sound,delay, response, lick, pITI, psound, pdelay, presponse, plick]= deal(nan(nROI,size(trialType,2)));
                    for roiNo = 1:size(obj.SavedCaTrials.f_raw{1},1)
                        disp([animal,date,'rioNo',num2str(roiNo)]);
                        [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,obj.dff(roiNo,:),frT,str_nFrames);
                        label_choice = fTrialType2Label(trialType(1:2,:,:),1);%only include cor and err trials
                        poslabel=1;
                        nshuffle=1000;
                        for nStim=1:size(trialType,2)
                            ind_trial=logical(reshape(sum(trialType(1:2,nStim,:),1),[],1));
                            [ITI(roiNo,nStim),pITI(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.ITI(ind_trial),poslabel,nshuffle);
                            [sound(roiNo,nStim),psound(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.sound(ind_trial),poslabel,nshuffle);
                            [delay(roiNo,nStim),pdelay(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.delay(ind_trial),poslabel,nshuffle);
                            [response(roiNo,nStim),presponse(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.response(ind_trial),poslabel,nshuffle);
                            [lick(roiNo,nStim),plick(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.lick(ind_trial),poslabel,nshuffle);
                        end
                        %calculate moving AUC
                        [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( obj.dff(roiNo,:), behEventFrameIndex,  frameNum );
                        binsize=1;
                        binstep=1;
                        label = label_choice;
                        for nStim=1:size(trialType,2)
                            indTrial=sum(trialType(1:2,nStim,:),1);%only for cor and err trials(1st d)
                            indTrial=logical(squeeze(indTrial));
                            [delayMovingAUC{roiNo}.stim{nStim},pdelayMovingAUC{roiNo}.stim{nStim}] = fMovingAUC(label(indTrial),dff_aligned(indTrial,:),poslabel,nshuffle,binsize,binstep);
                        end
                    end
                end
                %add task performance
                nROI=size(obj.SavedCaTrials.f_raw{1},1);
                session=[animal,'_',date,'_',field];
                trial2include='all';
                trial2exclude=[];
                objsession=Session2P(session,filepath,trial2include,trial2exclude);
                Tperformance=objsession.mSummaryBehavior(3);
                [ipsi_performance{1:nROI}]=deal(Tperformance{'ipsi','correct'});
                ipsi_performance=reshape(ipsi_performance,[],1);
                [contra_performance{1:nROI}]=deal(Tperformance{'contra','correct'});
                contra_performance=reshape(contra_performance,[],1);
                [overall_performance{1:nROI}]=deal(Tperformance{'total','correct'});
                overall_performance=reshape(overall_performance,[],1);
                varROI=(1:size(obj.SavedCaTrials.f_raw{1},1));
                varROI=reshape(varROI,[],1);
                TAUC=table(varanimal,vardate,varfield,varcelltype,varROItype,varROI,...
                    ipsi_performance, contra_performance,overall_performance,...
                    ITI, sound,delay, response, lick, mid_delay,late_delay,...
                    pITI, psound, pdelay, presponse, plick,pmid_delay,plate_delay,delayMovingAUC,pdelayMovingAUC,...
                    'VariableNames',{'animal','date','field','celltype','ROItype','nROI',...
                    'ipsi_performance', 'contra_performance','overall_performance',...
                    'ITI','sound','delay','response','lick','mid_delay','late_delay',...
                    'pITI','psound','pdelay','presponse','plick','pmid_delay','plate_delay','delayMovingAUC','pdelayMovingAUC'});
            end
            %calculate mean activities for each time epoch and whether they
            %different from ITI as baseline
            if ~exist('Tmean','var')
                nROI=size(obj.SavedCaTrials.f_raw{1},1);
                %     [ITI_m, sound_m,delay_m, response_m, lick_m, pITI_m, psound_m, pdelay_m, presponse_m, plick_m]= deal(nan(nROI,1));
                [trialType,~,~] = fGetTrialType( obj.Data_extract,[],trialTypeVar(i_selectivity),'matrix','left',combineCorErr);
                [trialType_byChoice,~,~] = fGetTrialType( obj.Data_extract,[],3,'matrix','left','combineCorErr');
                if contains(trialTypeStr,'cor') % 'cor'| 'cor and err'
                    ind_trial=logical(reshape(sum(trialType(1,:,:),2),[],1));%only correct trials
                end
                ind_ipsi=logical(reshape(sum(trialType_byChoice(1,1,:),2),[],1));
                ind_contra=logical(reshape(sum(trialType_byChoice(1,2,:),2),[],1));
                ind_tr_1=1;
                ntr=length(obj.SavedCaTrials.f_raw);
                frT = obj.SavedCaTrials.FrameTime;
                % align to behavior event
                nFrameEachTrial=cellfun(@(x) size(x,2),obj.SavedCaTrials.f_raw);
                ind_1stFrame=zeros(1,length(nFrameEachTrial));
                ind_1stFrame(1)=1;
                ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
                ind_1stFrame=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials
                frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
                frameNum=double(round(frameNumTime*1000/frT));
                [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( obj.Data_extract, ind_1stFrame, obj.SavedCaTrials.FrameTime );%get behavior event time
                [sound_pselectivity,delay_pselectivity,response_pselectivity,lick_pselectivity]=deal(ones(size(obj.SavedCaTrials.f_raw{1},1),1));
                for roiNo = 1:size(obj.SavedCaTrials.f_raw{1},1) %SavedCaTrials.nROIs may be not true
                    [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,obj.dff(roiNo,:),frT,str_nFrames);
                    T_SigbyEpoch_mean=T_SigbyEpoch(ind_trial,:);
                    [ITI_m(roiNo).mean,pITI_m(roiNo).mean]=fMeanPdiff(T_SigbyEpoch_mean.ITI,T_SigbyEpoch_mean.ITI);
                    [sound_m(roiNo).mean, psound_m(roiNo).mean]=fMeanPdiff(T_SigbyEpoch_mean.ITI,T_SigbyEpoch_mean.sound);
                    [delay_m(roiNo).mean,pdelay_m(roiNo).mean]=fMeanPdiff(T_SigbyEpoch_mean.ITI,T_SigbyEpoch_mean.delay);
                    [response_m(roiNo).mean,presponse_m(roiNo).mean]=fMeanPdiff(T_SigbyEpoch_mean.ITI,T_SigbyEpoch_mean.response);
                    [lick_m(roiNo).mean,plick_m(roiNo).mean]=fMeanPdiff(T_SigbyEpoch_mean.ITI,T_SigbyEpoch_mean.lick);
                    T_SigbyEpoch_ipsi=T_SigbyEpoch(logical(ind_trial.*ind_ipsi),:);
                    [ITI_m(roiNo).ipsi,pITI_m(roiNo).ipsi]=fMeanPdiff(T_SigbyEpoch_ipsi.ITI,T_SigbyEpoch_ipsi.ITI);
                    [sound_m(roiNo).ipsi, psound_m(roiNo).ipsi]=fMeanPdiff(T_SigbyEpoch_ipsi.ITI,T_SigbyEpoch_ipsi.sound);
                    [delay_m(roiNo).ipsi,pdelay_m(roiNo).ipsi]=fMeanPdiff(T_SigbyEpoch_ipsi.ITI,T_SigbyEpoch_ipsi.delay);
                    [response_m(roiNo).ipsi,presponse_m(roiNo).ipsi]=fMeanPdiff(T_SigbyEpoch_ipsi.ITI,T_SigbyEpoch_ipsi.response);
                    [lick_m(roiNo).ipsi,plick_m(roiNo).ipsi]=fMeanPdiff(T_SigbyEpoch_ipsi.ITI,T_SigbyEpoch_ipsi.lick);
                    T_SigbyEpoch_contra=T_SigbyEpoch(logical(ind_trial.*ind_contra),:);
                    [ITI_m(roiNo).contra,pITI_m(roiNo).contra]=fMeanPdiff(T_SigbyEpoch_contra.ITI,T_SigbyEpoch_contra.ITI);
                    [sound_m(roiNo).contra, psound_m(roiNo).contra]=fMeanPdiff(T_SigbyEpoch_contra.ITI,T_SigbyEpoch_contra.sound);
                    [delay_m(roiNo).contra,pdelay_m(roiNo).contra]=fMeanPdiff(T_SigbyEpoch_contra.ITI,T_SigbyEpoch_contra.delay);
                    [response_m(roiNo).contra,presponse_m(roiNo).contra]=fMeanPdiff(T_SigbyEpoch_contra.ITI,T_SigbyEpoch_contra.response);
                    [lick_m(roiNo).contra,plick_m(roiNo).contra]=fMeanPdiff(T_SigbyEpoch_contra.ITI,T_SigbyEpoch_contra.lick);
                    sound_pselectivity(roiNo)=ranksum(T_SigbyEpoch_contra.sound,T_SigbyEpoch_ipsi.sound);
                    delay_pselectivity(roiNo)=ranksum(T_SigbyEpoch_contra.delay,T_SigbyEpoch_ipsi.delay);
                    response_pselectivity(roiNo)=ranksum(T_SigbyEpoch_contra.response,T_SigbyEpoch_ipsi.response);
                    lick_pselectivity(roiNo)=ranksum(T_SigbyEpoch_contra.lick,T_SigbyEpoch_ipsi.lick);
                end
                [varanimal{1:nROI}]=deal(obj.metadata.animal);
                varanimal=reshape(varanimal,[],1);
                [vardate{1:nROI}]=deal(obj.metadata.date);
                vardate=reshape(vardate,[],1);
                [varfield{1:nROI}]=deal(obj.metadata.field);
                varfield=reshape(varfield,[],1);
                [varcelltype{1:nROI}]=deal(celltype);
                varcelltype=reshape(varcelltype,[],1);
                [varROItype{1:nROI}]=deal(ROItype);
                varROItype=reshape(varROItype,[],1);
                varROI=(1:size(obj.SavedCaTrials.f_raw{1},1));
                varROI=reshape(varROI,[],1);
                
                Tmean=table(varanimal,vardate,varfield,varcelltype,varROItype,varROI,...
                    ITI_m', sound_m',delay_m', response_m', lick_m', ...
                    pITI_m', psound_m', pdelay_m', presponse_m', plick_m',...
                    sound_pselectivity,delay_pselectivity,response_pselectivity,...
                    lick_pselectivity,'VariableNames',...
                    {'animal','date','field','celltype','ROItype','nROI',...
                    'ITI','sound','delay','response','lick',...
                    'pITI','psound','pdelay','presponse','plick',...
                    'sound_pselectivity','delay_pselectivity','response_pselectivity',...
                    'lick_pselectivity'});
            end
            save(fileNameT,'TAUC','Tmean');
            save(fileNameT_inSummary,'TAUC','Tmean');
            obj.TAUC=TAUC;
            obj.Tmean=Tmean;
        end

        function [TSVM,obj]=mSVMscore(obj,celltype,dataForm,indROI,strROI,trialTypeStr,...
                SVMtype,nRepeat, pTraining,CorrectedMethod,smoothBinsize,CTflag)
            %ref fGetEpochSVMscoreASession V22.5.26, to calculate a SVM score
            activity_data_smooth=obj.mSmoothFR(dataForm,smoothBinsize);%{objsession.zscored_spkr,objsession.dff,objsession.spkr}%choose among several acitivities form
            activity_data_raw=obj.mSmoothFR(dataForm,1);%when extracting epoch activities, do not smooth
            if ~isempty(indROI)
                activity_data=activity_data_smooth(indROI,:);%use part of the ROI to perform SVM
            else
                activity_data=activity_data_smooth;
            end
            indTrial2use=reshape(obj.metadata.indTrial2use,[],1);
            datestr=strrep(obj.metadata.date,'/','-');
            if strcmp(trialTypeStr,'cor')
                combineCorErr='divideCorErr';%and only use correct trial
            elseif strcmp(trialTypeStr,'cor and err')
                combineCorErr='combineCorErr';%{'combineCorErr','divideCorErr'}
            end
            selectivitystr={'stimuli','sensory difficulty','sensory','choice'};%sensory means grouping difficulties;
            trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
            i_selectivity=cellfun(@(x) strcmp(x,SVMtype),selectivitystr);%*********variable**************
            i_selectivity=find(i_selectivity);
            str_nFrames='1s';%'500ms';%'1s'
            frameNumTime_delay=[1,1.2];%from 1s before delay onset to 1s after that, since 1-1.5
            frameNumTime_go=[1,1.5];%from 5s before align point to 5s after align point
            %calculate moving SVM
            binsize=5;
            binstep=5;
            
            %naming the temp file by method of AUC calculation correction
            filepath_inSummary='E:\2P\summary\SVM\cases';%backup another copy in the summary folder
            fileName_datapars=[obj.metadata.animal,'-',datestr,'-',strROI,'_',dataForm,'-trialType',trialTypeStr,...
                '-epochBin',str_nFrames,'-Tdelay',num2str(frameNumTime_delay(1)),'-',num2str(frameNumTime_delay(2)),...
                '_go',num2str(frameNumTime_go(1)),'-',num2str(frameNumTime_go(2))];
            if strcmp(trialTypeStr,'cor') || strcmp(trialTypeStr,'err')
                combineCorErr='divideCorErr';%and only use correct trial
                corErrTrialNumber='raw';
                fileName_processpars=[SVMtype,'-smooth_binsize',num2str(smoothBinsize),'-movingTimeBin',num2str(binsize),'-',CTflag,'-pTraining',num2str(pTraining)];
            elseif strcmp(trialTypeStr,'cor and err')%used to test whether it is sensory or choice AUC,in order to be more fair, keep trial numbers balenced for correct and error trials
                combineCorErr='combineCorErr';%{'combineCorErr','divideCorErr'}
                fileName_processpars=[SVMtype,'-',CorrectedMethod,'-smooth_binsize',num2str(smoothBinsize),'-movingTimeBin',num2str(binsize),'-',CTflag,'-pTraining',num2str(pTraining)];
                switch CorrectedMethod
                    case 'balencedCorErrTrialNum'%probably problematic, since the trial are resampled and may let training and testing dataset be same
                        corErrTrialNumber='balence';
                    case 'SensoryChoiceOrthogonalSubtraction'
                        corErrTrialNumber='raw';
                        if (~strcmp(SVMtype,'choice'))&& (~strcmp(SVMtype,'sensory'))
                            warning('error input combination of fGetEpochAUCtableASession function, SensoryChoiceOrthogonalSubtraction method');
                        end
                    otherwise
                        corErrTrialNumber='raw';
                end
            end
            fileNameT=[obj.path.root,filesep,fileName_datapars,'-',fileName_processpars,'-EpochSVM.mat'];
            fileNameT_inSummary=[filepath_inSummary,filesep,fileName_datapars,'-',fileName_processpars,'-EpochSVM.mat'];
            
            if exist(fileNameT,'file')
                load(fileNameT);
                disp(['Table exist, use ',obj.metadata.animal,obj.metadata.date]);
            elseif exist(fileNameT_inSummary,'file')
                load(fileNameT_inSummary);
                disp(['Table exist, use ',obj.metadata.animal,obj.metadata.date]);
            else
                %
                disp(['analyzing',fileNameT]);
                nROI=1;%each session was viewed as one ROI
                varNROI=size(activity_data,1);
                [varindROI{1:nROI}]=deal(indROI);
                varindROI=reshape(varindROI,[],1);
                [varanimal{1:nROI}]=deal(obj.metadata.animal);
                varanimal=reshape(varanimal,[],1);
                [vardate{1:nROI}]=deal(obj.metadata.date);
                vardate=reshape(vardate,[],1);
                [varfield{1:nROI}]=deal(obj.metadata.field);
                varfield=reshape(varfield,[],1);
                [varcelltype{1:nROI}]=deal(celltype);
                varcelltype=reshape(varcelltype,[],1);
                
                ind_tr_1=1;
                ntr=length(obj.SavedCaTrials.f_raw);
                frT = obj.SavedCaTrials.FrameTime;
                % align to behavior event
                nFrameEachTrial=cellfun(@(x) size(x,2),obj.SavedCaTrials.f_raw);
                ind_1stFrame=zeros(1,length(nFrameEachTrial));
                ind_1stFrame(1)=1;
                ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
                ind_1stFrame=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials
                %time range for delay and go algnment
                frameNum_delay=double(round(frameNumTime_delay*1000/frT));
                ts_delay=-frameNumTime_delay(1):frT/1000:frameNumTime_delay(2);
                frameNum_go=double(round(frameNumTime_go*1000/frT));
                ts_go=-frameNumTime_go(1):frT/1000:frameNumTime_go(2);
                [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( obj.Data_extract, ind_1stFrame, obj.SavedCaTrials.FrameTime );%get behavior event time
                
                if strcmp(SVMtype,'choice')|| strcmp(SVMtype,'sensory')
                    switch corErrTrialNumber
                        case 'balence'
                            [trialType,~,~] = fGetTrialType( obj.Data_extract,[],trialTypeVar(i_selectivity),'matrix','left','divideCorErr');
                            trialTypeIndCell=cell(2,2);%{ipsi cor, contra cor;ipsi err, contra err};
                            for i=1:2
                                for j=1:2
                                    trialTypeIndCell{i,j}=find(logical(reshape(trialType(i,j,:),[],1)));
                                end
                            end
                            temp=cellfun(@length,trialTypeIndCell);
                            trialNumEachType=floor(sum(temp,'all')/4);%keep total trial number stable and redistributed to each trial types
                            s = RandStream('mlfg6331_64');%for reproducibility
                            trialTypeIndCellBalenced=cellfun(@(x) reshape(datasample(s,x,trialNumEachType,'Replace',true),[],1), trialTypeIndCell,'UniformOutput',false);
                            trialTypeIndFinal=cell2mat(trialTypeIndCellBalenced);%1st col-cor; 2nd col-err
                            ind_trial=reshape(trialTypeIndFinal,[],1);
                            ind_trial=logical(ind_trial.*indTrial2use);
                            label_SVM=[ones(trialNumEachType*2,1);ones(trialNumEachType*2,1)*2];
                        case 'raw'
                            [trialType,~,~] = fGetTrialType( obj.Data_extract,[],trialTypeVar(i_selectivity),'matrix','left',combineCorErr);
                            label_choice = fTrialType2Label(trialType,2);
                            if contains(trialTypeStr,'cor') % 'cor'| 'cor and err'
                                ind_trial=logical(reshape(sum(trialType(1,:,:),2),[],1));%only correct trials
                            elseif strcmp(trialTypeStr,'err')
                                ind_trial=logical(reshape(sum(trialType(2,:,:),2),[],1));%only correct trials
                            end
                            ind_trial=logical(ind_trial.*indTrial2use);
                            label_SVM=label_choice(ind_trial);
                            %used for orthogonal subtraction
                            [trialType_orthogonal,~,~] = fGetTrialType( obj.Data_extract,[],7-trialTypeVar(i_selectivity),'matrix','left',combineCorErr);
                            label_choice_orthogonal = fTrialType2Label(trialType_orthogonal,2);
                            label_SVM_orthogonal=label_choice_orthogonal(ind_trial);
                    end
                    
                    sigbyEpoch=[];
                    sigMoving=[];
                    sigMovingGo=[];
                    for roiNo = 1:size(activity_data,1)
                        %data for epoch activity
                        [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,activity_data_raw(roiNo,:),frT,str_nFrames);%when extracting epoch data, not smooth
                        if strcmp(CorrectedMethod,'SensoryChoiceOrthogonalSubtraction')
                            T_SigbyEpoch=fOrthogonalSubtraction(T_SigbyEpoch,ind_trial,label_SVM,label_SVM_orthogonal);
                        else
                            T_SigbyEpoch=T_SigbyEpoch(ind_trial,:);
                        end
                        sigbyEpoch=cat(3,sigbyEpoch,table2array(T_SigbyEpoch));
                        %data for moving activity
                        [ activity_data_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( activity_data(roiNo,:), behEventFrameIndex,  frameNum_delay );
                        if strcmp(CorrectedMethod,'SensoryChoiceOrthogonalSubtraction')
                            activity_data_aligned_ortho_corrected=fOrthogonalSubtraction(activity_data_aligned,ind_trial,label_SVM,label_SVM_orthogonal);
                        else
                            activity_data_aligned_ortho_corrected=activity_data_aligned(ind_trial,:);
                        end
                        sigMoving=cat(3,sigMoving,activity_data_aligned_ortho_corrected);
                        %data for moving activity around go cue
                        behEventAlign='go cue';
                        [ activity_data_aligned_go, ~, ~ ] = fAlignSigBehEvent( activity_data(roiNo,:), behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum_go );
                        if strcmp(CorrectedMethod,'SensoryChoiceOrthogonalSubtraction')
                            activity_data_aligned_ortho_corrected=fOrthogonalSubtraction(activity_data_aligned_go,ind_trial,label_SVM,label_SVM_orthogonal);
                        else
                            activity_data_aligned_ortho_corrected=activity_data_aligned_go(ind_trial,:);
                        end
                        sigMovingGo=cat(3,sigMovingGo,activity_data_aligned_ortho_corrected);
                    end
                    sigbyEpoch=permute(sigbyEpoch,[1,3,2]);%now, 1d trial number, 2d roi, 3d frames
                    sigMoving=permute(sigMoving,[1,3,2]);
                    sigMovingGo=permute(sigMovingGo,[1,3,2]);
                    
                    [CESVM.score, CESVM.score_shuffle,CESVM.p,CESVM.p_shuffled,CESVM.ts] = fCrossTimeSVM(sigbyEpoch,label_SVM, nRepeat, pTraining,1,1,1:5,CTflag);
                    %calculate moving SVM
                    if strcmp(CTflag,'CT') ||strcmp(CTflag,'TT')
                        [delayCTSVM.score, delayCTSVM.score_shuffle, delayCTSVM.p, delayCTSVM.p_shuffled, delayCTSVM.ts]=fCrossTimeSVM(sigMoving,label_SVM, nRepeat, pTraining,binsize,binstep,ts_delay,CTflag);
                        [goCTSVM.score, goCTSVM.score_shuffle, goCTSVM.p, goCTSVM.p_shuffled, goCTSVM.ts]=fCrossTimeSVM(sigMovingGo,label_SVM, nRepeat, pTraining,binsize,binstep,ts_go,CTflag);
                    else
                        [delayCTSVM.score, delayCTSVM.score_shuffle, delayCTSVM.p, delayCTSVM.p_shuffled, delayCTSVM.ts]=fCrossTimeSVM(sigMoving,label_SVM, nRepeat, pTraining,binsize,binstep,ts_delay,'nan');
                        [goCTSVM.score, goCTSVM.score_shuffle, goCTSVM.p, goCTSVM.p_shuffled, goCTSVM.ts]=fCrossTimeSVM(sigMovingGo,label_SVM, nRepeat, pTraining,binsize,binstep,ts_go,'nan');
                    end
                elseif strcmp(SVMtype,'stimuli') %here, compare auc of cor/err for each stimuli
                    [ITI, sound,delay, response, lick, pITI, psound, pdelay, presponse, plick]= deal(nan(nROI,size(trialType,2)));
                end
                
                TSVM=table(varanimal,vardate,varfield,varcelltype,varNROI,varindROI,CESVM,delayCTSVM,goCTSVM,...
                    'VariableNames',{'animal','date','field','celltype','nROI','indROI','CESVM','delayCTSVM','goCTSVM'});
                %}
            end
            save(fileNameT,'TSVM');
            save(fileNameT_inSummary,'TSVM');%save another copy in the summary folder
            obj.TSVM=TSVM;
        end
        
        function [figMeanSig] = mPlotCrossTimeSVM(obj,CTSVM_label,pSig,SVMtype,varargin)
            %plot accuracy of each session
            switch CTSVM_label
                case 'delayCTSVM'
                    CTSVM=obj.TSVM.delayCTSVM;
                case 'CESVM'
                    CTSVM=obj.TSVM.CESVM;
                case 'goCTSVM'
                    CTSVM=obj.TSVM.goCTSVM;
            end
            score=CTSVM.score;
            score_shuffle=CTSVM.score_shuffle;
            pCT=CTSVM.p;
            pCT_shuffle=CTSVM.p_shuffled;
            [score_tt,score_shuffle_tt]=deal(zeros(size(score,1),size(score,2)));
            for i=1:size(score,1)
                score_tt(i,:)=diag(squeeze(score(i,:,:)));
                score_shuffle_tt(i,:)=diag(squeeze(score_shuffle(i,:,:)));
            end
            pTT=diag(pCT);
            pTT_shuffle=diag(pCT_shuffle);
            
            figMeanSig=figure;
            set(gcf,'position',[200,200,700,200]);
            colorScore={[1,0,0],[0,0,1],[0.5,0.5,0.5]};
            
            subplot(1,3,1);
            if strcmp(SVMtype,'choice')
                fPlotMean_CI(CTSVM.ts,score_tt,colorScore{1},pSig);
            elseif strcmp(SVMtype,'sensory')
                fPlotMean_CI(CTSVM.ts,score_tt,colorScore{2},pSig);
            end
            hold on;
            fPlotMean_CI(CTSVM.ts,score_shuffle_tt,colorScore{3},pSig);
            %label significance
            indSig=pTT<0.05;
            markersize=fSigMarkerSize(pTT);
            y_lim=get(gca,'Ylim');
            ySig=ones(size(pTT))*y_lim(2)*0.9;
            xSig=CTSVM.ts(1:length(pTT));
            xSig(~indSig)=[];
            ySig(~indSig)=[];
            markersize(~indSig)=[];
            if strcmp(SVMtype,'choice')
                scatter(xSig,ySig,'Marker','*','SizeData', markersize,'MarkerEdgeColor','r');
            elseif strcmp(SVMtype,'sensory')
                scatter(xSig,ySig,'Marker','*','SizeData', markersize,'MarkerEdgeColor','b');
            end
            hold on;
            %color plot of the cross temporal decoding
            ax2=subplot(1,3,2);
            fPlotMeanSig2D_SVM(ax2, score,pCT,CTSVM.ts,pSig,'no bar');
            ylabel('Training time');
            xlabel('Testing time');
            title('Real data');
            ax3=subplot(1,3,3);
            fPlotMeanSig2D_SVM(ax3, score_shuffle,pCT_shuffle, CTSVM.ts,pSig,'colorbar');
            ylabel('Training time');
            xlabel('Testing time');
            title('Shuffled data');
            %for epoch SVM, label epochs
            if ~isempty(varargin)
                ts_label=varargin{1};
                subplot(1,3,1);
                set(gca,'XTick',1:length(ts_label),'XTickLabel',ts_label);
                for i=2:3
                    subplot(1,3,i);
                    set(gca,'XTick',1:length(ts_label),'XTickLabel',fliplr(ts_label));
                    set(gca,'YTick',1:length(ts_label),'YTickLabel',fliplr(ts_label));
                end
            end
        end
    end
    
end

%% auxiliary fuction
function [timeTrialStart]=fTrialStartStr2double(strTrialStart)
%transform the trial start time from string form to double form
t=strsplit(strTrialStart,'_');
t=str2double(t);
timeTrialStart=t(4)*3600+t(5)*60+t(6)+t(7)/1000; %time(s) relative to trial1
end

function [out]=fBaselineCorrection(in,span)
%span = round(20*1000/frT/2); %every 40s, get a mean and substrate to compensate for baseline change
x = [];
for i = 1:length(in)
    %     ind_x = ind_x + 1;
    ind1 = max(i- span,1);
    ind2 = min(i+ span,length(in));
    x(i) = prctile(in(ind1:ind2),8);
end
x=reshape(x,[],1);
out = in - x + nanmean(x);
end

function [meanout,pout] = fBootstrpMeanP(activity,nboot,baseline)
bootmean=bootstrp(nboot,@nanmean,activity);
meanout=nanmean(bootmean);
pout=sum(bootmean<baseline)./sum(~isnan(bootmean));
pout(pout<0.5)=pout(pout<0.5)*2;
pout(pout>=0.5)=(1-pout(pout>=0.5))*2;
end

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

function markersize=fSigMarkerSize(p)
markersize(p>0.05)=0;
markersize(logical((p<=0.05).*(p>0.01)))=4;
markersize(logical((p<=0.01).*(p>0.001)))=8;
markersize(logical(p<=0.001))=12;
end

function [mean, p] = fMeanPdiff(baseline,data,varargin)
%Input- matrix of baseline/data, n-by-m (n repeats, m time points)
%   method- 'ranksum'(default)|'ttest', etc
%Output- vector of mean, p, 1-by-m
if isempty(varargin)
    method ='ranksum';
else
    method =varargin{1};
end
mean=nanmean(data,1);
switch method
    case 'ranksum'
        if size(data,2)~=size(baseline,2)
            warnning('input data not same column number');
            return;
        end
        for i=1:size(data,2)
            p(i)=ranksum(baseline(:,i), data(:,i));
        end
    case 'ttest'
        [~,p]=ttest(baseline,data);
end
end