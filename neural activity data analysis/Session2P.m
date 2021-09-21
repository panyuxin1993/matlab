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
        SavedCaTrials
        DLC         %struct read from a .xls file, including dcnum, dctxt
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
            if length(temp)>2
                obj.metadata.field=temp{3};
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
                            masklick,i_selectivity,behEventSort,obj.metadata.trial2include,obj.metadata.trial2exclude,...
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
                            masklick,i_selectivity,behEventSort,obj.metadata.trial2include,obj.metadata.trial2exclude,...
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
                            masklick,i_selectivity,behEventSort,obj.metadata.trial2include,obj.metadata.trial2exclude,...
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
        
        function meanActivityByTrialType=mMeanPSTHbyROI(obj,activity_type,ind_ROI,behEventAlign,masklick,i_selectivity)
            if strcmp(activity_type,'dff')
                activities_data=obj.zscored_dff;
            elseif strcmp(activity_type,'spkr')
                activities_data=obj.zscored_spkr;
            end
            %plotting settings
            if strcmp(behEventAlign,'stim onset')
                frameNumTime=[0.5,2];%from 5s before align point to 5s after align point
            elseif strcmp(behEventAlign,'delay onset')
                frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
            else
                frameNumTime=[2,3.5];%from 5s before align point to 5s after align point
            end
            if isempty(ind_ROI) || min(ind_ROI)>obj.metadata.nROI || max(ind_ROI)<1
                iStart=1;
                iEnd=obj.metadata.nROI;
            else
                iStart=max(1,min(ind_ROI));
                iEnd=min(obj.metadata.nROI,max(ind_ROI));
            end
            for iROI=iStart:iEnd
                meanActivityByTrialType_ROI = fMeanPSTHbyROI(activities_data(iROI,:),...
                    obj.metadata.ind_tr_1,obj.Data_extract,obj.SavedCaTrials,...
                    frameNumTime,behEventAlign,masklick,i_selectivity,...
                    obj.metadata.trial2include,obj.metadata.trial2exclude);
                if exist('meanActivityByTrialType','var')
                    meanActivityByTrialType=cellfun(@(x,y) vertcat(x,y),meanActivityByTrialType,meanActivityByTrialType_ROI,'UniformOutput',false);
                else
                    meanActivityByTrialType=meanActivityByTrialType_ROI;
                end
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

        function TSVM=mSVMscore(obj,celltype,trialTypeStr,ind_NSDelayMovement,SVMtype,nRepeat, pTraining)
            %ref fGetEpochSVMscoreASession, to calculate a SVM score
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
            str_nFrames='500ms';%'500ms';%'1s'
            fileNameT=[obj.path.root,filesep,obj.metadata.animal,'-',datestr,'trialType',trialTypeStr,'-',SVMtype,'-timeBin',str_nFrames,'-pTraining',num2str(pTraining),'-EpochSVM.mat'];
            if exist(fileNameT,'file')
                load(fileNameT);
                disp(['Table exist, use ',obj.name]);
            else
                disp(['analyzing',fileNameT]);
                nROI=1;%each session was viewed as one ROI
                [delayMovingSVM,pdelayMovingSVM]=deal([]);
                [varanimal{1:nROI}]=deal(obj.metadata.animal);
                varanimal=reshape(varanimal,[],1);
                [vardate{1:nROI}]=deal(obj.metadata.date);
                vardate=reshape(vardate,[],1);
                [varfield{1:nROI}]=deal(obj.name);
                varfield=reshape(varfield,[],1);
                [varcelltype{1:nROI}]=deal(celltype);
                varcelltype=reshape(varcelltype,[],1);
                frT = obj.SavedCaTrials.FrameTime;
                % align to behavior event
                frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
                frameNum=double(round(frameNumTime*1000/frT));
                [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( obj.Data_extract, obj.metadata.ind_1stFrame, obj.SavedCaTrials.FrameTime );%get behavior event time
                [trialType,~,~] = fGetTrialType( obj.Data_extract,[],trialTypeVar(i_selectivity),'matrix','left',combineCorErr);
                f
                if strcmp(SVMtype,'choice')|| strcmp(SVMtype,'sensory')
                    [ITI, sound,delay, response,lick,pITI, psound, pdelay, presponse, plick] = deal(cell(ones(1,1)));
                    label_choice = fTrialType2Label(trialType,2);
                    if contains(trialTypeStr,'cor') % 'cor'| 'cor and err'
                        ind_trial=logical(reshape(sum(trialType(1,:,:),2),[],1));%only correct trials
                    end
                    if ~empty(ind_NSDelayMovement)%for sessions without DLC data and no index of trials selected, just use all data
                        ind_trial=logical(ind_trial.*reshape(ind_NSDelayMovement,[],1));%choose only the trials without significant movement during the delay
                    end
                    sigbyEpoch=[];
                    sigMoving=[];
                    for roiNo = 1:size(obj.SavedCaTrials.f_raw{1},1)
                        %data for epoch activity
                        [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,obj.dff(roiNo,:),frT,str_nFrames);
                        sigbyEpoch=cat(3,sigbyEpoch,table2array(T_SigbyEpoch));
                        %data for moving activity
                        [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( obj.dff(roiNo,:), behEventFrameIndex,  frameNum );
                        sigMoving=cat(3,sigMoving,dff_aligned);
                    end
                    sigbyEpoch=permute(sigbyEpoch,[1,3,2]);
                    sigMoving=permute(sigMoving,[1,3,2]);
                    
                    score_shuffle=cell(1,1);
                    [score,score_shuffle{1}] = fSVM(sigbyEpoch(ind_trial,:,:),label_choice(ind_trial), nRepeat, pTraining,1,1);
                    [meanout,pout] = fBootstrpMeanP(score,1000,0.5);
                    [ITI{1}, sound{1},delay{1}, response{1},lick{1}] = deal(score(:,1),score(:,2),score(:,3),score(:,4),score(:,5));
                    [pITI{1}, psound{1},pdelay{1}, presponse{1},plick{1}] = deal(pout(:,1),pout(:,2),pout(:,3),pout(:,4),pout(:,5));
                    %calculate moving SVM
                    binsize=3;
                    binstep=1;
                    label = label_choice;
                    for nResult=1:size(trialType,1)-2
                        indTrial=trialType(nResult,:,:);
                        indTrial=sum(squeeze(indTrial),1);
                        indTrial=logical(squeeze(indTrial));
                        if nResult==1
                            if strcmp(trialTypeStr,'cor and err')
                                [delayMovingSVM.do,delayMovingSVM.shuffle_do]=fSVM(sigMoving(indTrial,:,:),label_choice(indTrial), nRepeat, pTraining,binsize,binstep);
                                [~,pdelayMovingSVM.do]=fBootstrpMeanP(delayMovingSVM.do,1000,0.5);
                            else
                                [delayMovingSVM.cor,delayMovingSVM.shuffle_cor]=fSVM(sigMoving(indTrial,:,:),label_choice(indTrial), nRepeat, pTraining,binsize,binstep);
                                [~,pdelayMovingSVM.cor]=fBootstrpMeanP(delayMovingSVM.cor,1000,0.5);
                            end
                        elseif nResult==2
                            [delayMovingSVM.err,delayMovingSVM.shuffle_err]=fSVM(sigMoving(indTrial,:,:),label_choice(indTrial), nRepeat, pTraining,binsize,binstep);
                            [~,pdelayMovingSVM.err]=fBootstrpMeanP(delayMovingSVM.err,1000,0.5);
                        end
                    end
                    
                elseif strcmp(SVMtype,'stimuli') %here, compare auc of cor/err for each stimuli
                    [ITI, sound,delay, response, lick, pITI, psound, pdelay, presponse, plick]= deal(nan(nROI,size(trialType,2)));
                    %         for roiNo = 1
                    %             [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,dff(roiNo,:),frT,str_nFrames);
                    %             label_choice = fTrialType2Label(trialType(1:2,:,:),1);%only include cor and err trials
                    %             poslabel=1;
                    %             nshuffle=1000;
                    %             for nStim=1:size(trialType,2)
                    %                 ind_trial=logical(reshape(sum(trialType(1:2,nStim,:),1),[],1));
                    %                 [ITI(roiNo,nStim),pITI(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.ITI(ind_trial),poslabel,nshuffle);
                    %                 [sound(roiNo,nStim),psound(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.sound(ind_trial),poslabel,nshuffle);
                    %                 [delay(roiNo,nStim),pdelay(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.delay(ind_trial),poslabel,nshuffle);
                    %                 [response(roiNo,nStim),presponse(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.response(ind_trial),poslabel,nshuffle);
                    %                 [lick(roiNo,nStim),plick(roiNo,nStim)]=fAUC(label_choice(ind_trial),T_SigbyEpoch.lick(ind_trial),poslabel,nshuffle);
                    %             end
                    %             %calculate moving AUC
                    %             [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff(roiNo,:), behEventFrameIndex,  frameNum );
                    %             binsize=1;
                    %             binstep=1;
                    %             label = label_choice;
                    %             for nStim=1:size(trialType,2)
                    %                 indTrial=sum(trialType(1:2,nStim,:),1);%only for cor and err trials(1st d)
                    %                 indTrial=logical(squeeze(indTrial));
                    %                 [delayMovingSVM.stim{nStim},pdelayMovingSVM.stim{nStim}] = fMovingAUC(label(indTrial),dff_aligned(indTrial,:),poslabel,nshuffle,binsize,binstep);
                    %             end
                    %         end
                end
                
                TSVM=table(varanimal,vardate,varfield,varcelltype,ITI, sound,delay, response, ...
                    lick, pITI, psound, pdelay, presponse, plick,delayMovingSVM,pdelayMovingSVM,score_shuffle,'VariableNames',...
                    {'animal','date','field','celltype','ITI','sound','delay','response',...
                    'lick','pITI','psound','pdelay','presponse','plick','delayMovingSVM','pdelayMovingSVM','shuffleEpochSVMscore'});
            end
            save(fileNameT,'TSVM');
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