%continuing sampling, so link each trials one by one to get a 'long trial'
% load('F:\pyx083-20180522\im_data_reg_cpu\result_save\CaTrialsSIM_pyx083_20180522_920nm_power50_zoom4x_dftReg_.mat');
% load('D:\xulab\behavior\pyx083\Data_Virables\2018_05_22_pyx083-imaging_Virables.mat');
clear;
close all;
%plot data from summary xlsx file, batch analysis
%{
[num,txt,raw] =xlsread('D:\xulab\project\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
T=cell2table(raw(2:end,1:14));
T.Properties.VariableNames=strrep(raw(1,1:14),' ','_');%table variable name can't have ' ',so replace them
ind_session=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control');
ind_session=find(ind_session);
n_session=length(ind_session);
animal_unique=unique(T.animal(ind_session));
n_animal=length(animal_unique);
sessionNamePool=cell(n_session,1);
pathPool=cell(n_session,1);
for i=1:n_session
    sessionNamePool{i}=T.session{ind_session(i)};
    pathPool{i}=T.file_path{ind_session(i)};
    trial2includecell{i}=T.used_trial{ind_session(i)};
end
% trial2include=[1,64;69,82;93,215];
trial2exclude=[];
%}
%plot single session
%
sessionNamePool={'pyx340_20210202'};
pathPool={'D:\pyx340_20210202\im_data_reg\result_save'};
trial2includecell={[1,193]};
%}

behEventAlignPool={'delay onset','go cue','first lick'};
masklickPool={'yes','no','no'};
% behEventAlignPool={'go cue','first lick'};
% masklickPool={'no','no'};

for i_align=1:length(behEventAlignPool)
    %decide some global variable
    behEventAlign=behEventAlignPool{i_align};%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
    behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
    masklick=masklickPool{i_align};
    
    i_selectivity=1;%*********variable**************
    selectivitystr={'stimuli','sensory difficulty','sensory','choice'};%sensory means grouping difficulties;
    trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
    trialTypeStr=selectivitystr{i_selectivity};
    
    for i_session=1:length(sessionNamePool)
        sessionName=sessionNamePool{i_session};
        datapath=['D:\',sessionName];
        
        
        %plotting settings
        if strcmp(behEventAlign,'stim onset')
            frameNumTime=[0.5,2];%from 5s before align point to 5s after align point
        elseif strcmp(behEventAlign,'delay onset')
            frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
        else
            frameNumTime=[1,2];%from 5s before align point to 5s after align point
        end
        % yrange=[ ];
        clear CurrFolder;
        cd([datapath,'\im_data_reg\result_save']);
%         cd(pathPool{i_session});
        CurrFolder=pwd;
        savefolder=sessionName;%'1-200trials';%'segNP';
        % savefolder='F:\2P\example\';
        
        % load([CurrFolder,'\',sessionName, '_imaging_Virables.mat']);%load behavior data
        dirmat=strcat(CurrFolder,'\*.mat');
        dirs=dir(dirmat);
        dircell=struct2cell(dirs);
        filenames=dircell(1,:);
        file_imaging=cellfun(@(x) contains(x,'Ca'), filenames);
        i_file_imaging=find(file_imaging);
        load(filenames{1,i_file_imaging});%load imaging data
        file_beh=cellfun(@(x) contains(x,'imaging_Virables'), filenames);
        i_file_beh=find(file_beh);
        load(filenames{1,i_file_beh});%load behavior data
        
        if ~exist(savefolder)
            mkdir(savefolder);
        end
        
        if exist('dff.mat','file')
            load([CurrFolder,'\','dff.mat']);%load dff
        else
            dff=fnx_getDff(CurrFolder,sessionName,'save figure');
        end
        
        ind_tr_1=1;%using data from trial 1
        ntr = length(SavedCaTrials.f_raw); %use all data, or  just set the number of trials to use
        frT = SavedCaTrials.FrameTime;
        % align to behavior event
        nFrameEachTrial=cellfun(@(x) size(x,2),SavedCaTrials.f_raw);
        ind_1stFrame=zeros(1,length(nFrameEachTrial));
        ind_1stFrame(1)=1;
        ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
        ind_1stFrame=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials
        if isnan(trial2includecell{i_session})            
            indTrial2include=fExcludeTrials(trial2exclude,ind_1stFrame,'logical');
        else
            indTrial2include=fIncludeTrials(trial2includecell{i_session},ind_1stFrame,'logical');
        end
        frameNum=double(round(frameNumTime*1000/frT));
        [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( Data_extract, ind_1stFrame, SavedCaTrials.FrameTime ,ind_tr_1);%get behavior event time
        
        
        for roiNo = 1:size(SavedCaTrials.f_raw{1},1) %SavedCaTrials.nROIs may be not true
            %replot f_raw picture, labeled with used trial range
            trial2include=[1,length(ind_1stFrame)];
            fLabelDffIncludedTrialRange(dff(roiNo,:),[CurrFolder,'\',savefolder,'\',savefolder,'ROI-',num2str(roiNo),'-cat_f_label_used_range.jpg'],ind_1stFrame,trial2include)
            %   align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
            if strcmp(behEventAlign,'delay onset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
                [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff(roiNo,:), behEventFrameIndex,  frameNum );
            else
                [ dff_aligned, behEvent_aligned, licking_aligned ] = fAlignSigBehEvent( dff(roiNo,:), behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
            end
            
            [trialType,rule] = fGetTrialType( Data_extract,[],i_selectivity,'matrix','left','divideCorErr');%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
            disp(strcat('processing...',sessionName,'-ROI',num2str(roiNo)));
            %%
            %plot color plot,4 column(correct/error/miss/violation), 7 row(for each
            %stimlus)+ mean trace
            A=figure;
            set(gcf, 'position', [0 0 1400 200*(size(trialType,2)+1)]);%last row, mean trace
            ColLimit = prctile(dff(roiNo,:),90);
            titlestr={'Correct','Error','Miss','Violation'};
            titlestr=strcat('ROI-',num2str(roiNo),titlestr);
            if size(trialType,2)==6
                color_mean_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};
            elseif size(trialType,2)==2
                color_mean_trace={[0 0 1],[1 0 0]};
            end
            %     if strcmp(rule,'low click rate-right')
            %         color_mean_trace=fliplr(color_mean_trace);
            %     end
            for nStim=1:size(trialType,2) %for each stimulus
                for  nResult=1:size(trialType,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
                    selectedTrialInd=trialType(nResult,nStim,:);
                    selectedTrialInd=logical(squeeze(selectedTrialInd))';
                    selectedTrialInd=logical(selectedTrialInd.*indTrial2include');
                    sot=behEvent_aligned.stimOnset(selectedTrialInd);% stim onset, white
                    flt=behEvent_aligned.lickFirst(selectedTrialInd);% first lick, black
                    flt_l=behEvent_aligned.lickFirst_left(selectedTrialInd);%left lick, magenta dots
                    flt_r=behEvent_aligned.lickFirst_right(selectedTrialInd);%right lick, red dots
                    rwt=behEvent_aligned.rewTime(selectedTrialInd);%reward, cyan dots
                    go=behEvent_aligned.go(selectedTrialInd);%go cue,white
                    neuralActivity=dff_aligned(selectedTrialInd,:);
                    temp=find(selectedTrialInd);
                    clear lickTime;
                    for i=1:length(temp)
                        lickTime(i).lickleft= licking_aligned.leftLick{temp(i)};%unit, frame number
                        lickTime(i).lickright=licking_aligned.rightLick{temp(i)};
                    end
                    main_licking_dir=Data_extract.Action_choice(temp(i));%this trial type, the mainly licking direction
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
                    figure(A);
                    %             if nStim == 1 %plot neural activity; here end trials have more trials so wider
                    %                 subplot(size(trialType,2)+3,2*size(trialType,1),[2*nResult-1,2*nResult-1+2*size(trialType,1)]);
                    %             elseif nStim  == size(trialType,2)
                    %                 subplot(size(trialType,2)+3,2*size(trialType,1),[2*nResult-1+2*size(trialType,1)*size(trialType,2),2*nResult-1+2*size(trialType,1)*size(trialType,2)+2*size(trialType,1)]);
                    %             else
                    %                 subplot(size(trialType,2)+3,2*size(trialType,1),2*nResult-1+2*size(trialType,1)*nStim);
                    %             end
                    subplot(size(trialType,2)+1,2*size(trialType,1),2*nResult-1+2*size(trialType,1)*(nStim-1));
                    hold on;
                    imagesc(neuralActivity(I,:));
                    x_lim=[0,size(neuralActivity,2)];%get(gca,'Xlim');
                    plot(sot,[1:sum(selectedTrialInd)],'w.','markersize',8);
                    hold on;
                    plot(go,[1:sum(selectedTrialInd)],'w.','markersize',8);
                    hold on;
                    plot(flt,[1:sum(selectedTrialInd)],'k.','markersize',8);
                    hold on;
                    plot(flt_l,[1:sum(selectedTrialInd)],'m.','markersize',8);
                    hold on;
                    plot(flt_r,[1:sum(selectedTrialInd)],'r.','markersize',8);
                    hold on;
                    plot(rwt,[1:sum(selectedTrialInd)],'c.','markersize',8);
                    hold on;
                    set(gca,'clim',[0 ColLimit]);
                    set(gca,'ytick',sum(selectedTrialInd),'yticklabel',sum(selectedTrialInd),'ydir','normal');
                    if nStim==size(trialType,2)
                        set(gca,'xtick',[round(1000/frT*(frameNumTime(1)-floor(frameNumTime(1))))+1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
                        %                 xlabel(['time(s) from ',behEventAlign],'FontName','Arial','FontSize',14);
                    else
                        set(gca,'xtick',[]);
                    end
                    ColBar = colorbar;
                    set(ColBar,'position',[0.91 0.1100 0.02 0.1445]);
                    if sum(selectedTrialInd)>0
                        ylim([0.5 sum(selectedTrialInd)+0.5]);
                    end
                    if nStim == 1
                        title(titlestr{nResult});
                    end
                    if nResult==1
                        if strcmp(Data_extract.rule,'low click rate-right')
                            ylabel([num2str(Data_extract.Stimuli(size(trialType,2)+1-nStim))]);
                        else
                            ylabel([num2str(Data_extract.Stimuli(nStim))]); %' clicks/s'
                        end
                    end
                    xlim(x_lim);
                    
                    %             if nStim == 1 %plot licking raster
                    %                 subplot(size(trialType,2)+3,2*size(trialType,1),[2*nResult,2*nResult+2*size(trialType,1)]);
                    %             elseif nStim  == size(trialType,2)
                    %                 subplot(size(trialType,2)+3,2*size(trialType,1),[2*nResult+2*size(trialType,1)*size(trialType,2),2*nResult+2*size(trialType,1)*size(trialType,2)+2*size(trialType,1)]);
                    %             else
                    %                 subplot(size(trialType,2)+3,2*size(trialType,1),2*nResult+2*size(trialType,1)*nStim);
                    %             end
                    subplot(size(trialType,2)+1,2*size(trialType,1),2*nResult+2*size(trialType,1)*(nStim-1));
                    if exist('lickTime')
                        for i=1:length(lickTime)%each trial as a row
                            for jl=1:length(lickTime(I(i)).lickleft)
                                line([lickTime(I(i)).lickleft(jl),lickTime(I(i)).lickleft(jl)],[i-0.5,i+0.5],'color','b','linewidth',1);%left lick
                                hold on;
                            end
                            for jr=1:length(lickTime(I(i)).lickright)
                                line([lickTime(I(i)).lickright(jr),lickTime(I(i)).lickright(jr)],[i-0.5,i+0.5],'color','r','linewidth',1);%right lick
                                hold on;
                            end
                        end
                        set(gca,'ytick',[]);%licking raster have corresponding color plot label, so not need
                        if nStim==size(trialType,2)
                            %set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
                            %set(gca,'xtick',[-round(1000/frT*frameNumTime(1))+1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
                            set(gca,'xtick',[round(1000/frT*(frameNumTime(1)-floor(frameNumTime(1))))+1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
                        else
                            set(gca,'xtick',[]);
                        end
                        xlim(x_lim);
                        if sum(selectedTrialInd)>0
                            ylim([0.5 sum(selectedTrialInd)+0.5]);
                        end
                    end
                    
                    %plot mean trace
                    subplot(size(trialType,2)+1,2*size(trialType,1),2*nResult-1+2*size(trialType,1)*size(trialType,2));
                    %[neuralActivityMean, neuralActivityCI]=fMean_CI(neuralActivity,0.05);
                    ts=double((-frameNum(1):1:frameNum(2))*frT/1000);
                    curve_meanTrace(nStim)=fPlotMean_SE(ts,neuralActivity,color_mean_trace{nStim});
                    
                    %plot mean trace of licking raster
                    if exist('lickTime')
                        subplot(size(trialType,2)+1,2*size(trialType,1),2*nResult+2*size(trialType,1)*size(trialType,2));
                        binSize=5;%frame number
                        binStep=1;
                        [ lickMat,lickRate ] = fLickRate( lickTime, binSize,sum(frameNum)+1,0, ones(length(lickTime),1), binStep);%r=unit ms, so need transform to frames
                        switch main_licking_dir
                            case 0
                                i_lickrate=2;
                            case 1
                                i_lickrate=3;
                            otherwise
                                i_lickrate=1;
                        end
                        ts=double((-frameNum(1):binStep:frameNum(2))*frT/1000);
                        indx=min(length(ts),size(lickRate,1));
                        curve_lick_meanTrace(nStim)=plot(ts(1:indx),lickRate(1:indx,i_lickrate)/frT,'color',color_mean_trace{nStim});hold on;
                    end
                end
            end
            
            %% label and general plots
            % find common ylim
            for  nResult=1:size(trialType,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
                subplot(size(trialType,2)+1,2*size(trialType,1),2*nResult-1+2*size(trialType,1)*size(trialType,2));
                flagclearyrange=0;
                if ~exist('yrange','var')
                    flagclearyrange=1;
                    if ~exist('y_lim','var')
                        y_lim=get(gca,'Ylim');
                    else
                        y_lim2=get(gca,'Ylim');
                        y_lim=[min(y_lim(1),y_lim2(1)),max(y_lim(2),y_lim2(2))];
                    end
                else
                    y_lim=yrange;
                end
            end
            for  nResult=1:size(trialType,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
                subplot(size(trialType,2)+1,2*size(trialType,1),2*nResult-1+2*size(trialType,1)*size(trialType,2));
                xlim([ts(1),ts(end)]);
                plot([0,0],[y_lim(1),y_lim(2)],'k-');%align to a behavior event
                if strcmp(behEventAlign,'stim onset')
                    plot(0.5*ones(1,2),[y_lim(1),y_lim(2)],'k-');%align to a behavior event
                elseif strcmp(behEventAlign,'delay onset')
                    plot(-0.5*ones(1,2),[y_lim(1),y_lim(2)],'k-');%align to a behavior event
                end
                xlabel(['time(s) from ',behEventAlign],'FontSize',12);
                if nResult==1
                    ylabel('\it\DeltaF/F');
                end
                %set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
                set(gca,'Ylim',y_lim,'xtick',[-floor(frameNumTime(1)):1:frameNumTime(2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
                %         set(gca,'FontName','Arial','FontSize',14);
            end
            
            % label and general plots
            % find common ylim
            clear y_lim;
            for  nResult=1:size(trialType,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
                subplot(size(trialType,2)+1,2*size(trialType,1),2*nResult+2*size(trialType,1)*size(trialType,2));
                flagclearyrange=0;
                if ~exist('yrange','var')
                    flagclearyrange=1;
                    if ~exist('y_lim','var')
                        y_lim=get(gca,'Ylim');
                    else
                        y_lim2=get(gca,'Ylim');
                        y_lim=[min(y_lim(1),y_lim2(1)),max(y_lim(2),y_lim2(2))];
                    end
                else
                    y_lim=yrange;
                end
            end
            for  nResult=1:size(trialType,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
                subplot(size(trialType,2)+1,2*size(trialType,1),2*nResult+2*size(trialType,1)*size(trialType,2));
                xlim([ts(1),ts(end)]);
                plot([0,0],[y_lim(1),y_lim(2)],'k-');%align to a behavior event
                if strcmp(behEventAlign,'stim onset')
                    plot(0.5*ones(1,2),[y_lim(1),y_lim(2)],'k-');%align to a behavior event
                elseif strcmp(behEventAlign,'delay onset')
                    plot(-0.5*ones(1,2),[y_lim(1),y_lim(2)],'k-');%align to a behavior event
                end
                %         xlabel(['time(s) from ',behEventAlign],'FontSize',12);
                if nResult==1
                    ylabel('licks/s');
                end
                box off;
                %set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
                set(gca,'Ylim',y_lim,'xtick',[-floor(frameNumTime(1)):1:frameNumTime(2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
                %         set(gca,'FontName','Arial','FontSize',14);
            end
            
            subplot(size(trialType,2)+1,2*size(trialType,1),2*size(trialType,1)*size(trialType,2)+1);
            if contains(trialTypeStr,'stimuli')
                %         h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','','contra hard','','contra easy'},'Location','best');
                h=legend(curve_meanTrace(:),num2str(Data_extract.Stimuli(:)),'Location','best');
            elseif contains(trialTypeStr,'difficulty')
                h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','contra hard','contra easy'},'Location','best');
            elseif contains(trialTypeStr,'choice')
                h=legend(curve_meanTrace(:),{'ipsi choice','contra choice'},'Location','best');
            elseif contains(trialTypeStr,'sensory')
                h=legend(curve_meanTrace(:),{'low click rate','high click rate'},'Location','best');
            end
            set(h,'box','off');
            clear curve_meanTrace curve_lick_meanTrace;
            saveas(A,[CurrFolder,'\',savefolder,'\',sessionName,'ROI-',num2str(roiNo),'-algin to ',behEventAlign,'-sort ',behEventSort,'colorplot.jpg'],'jpg');
            
            % set(A,'PaperPosition',[0,0,4,2]);
            % saveas(A,['F:\2P\example\',sessionName,'ROI-',num2str(roiNo),'-algin to ',behEventAlign,'-sort ',behEventSort,'colorplot.pdf'],'pdf');
            % set(figMeanTrace,'PaperPosition',[0,0,4,1.5]);
            % saveas(figMeanTrace,['H:\2P\example\',sessionName,'ROI-',num2str(roiNo),'-algin to ',behEventAlign,'-sort ',behEventSort,'mean trace.pdf'],'pdf');
            close all;
            if flagclearyrange==1
                clear y_lim;
            end
        end
        
    end
end
%%
function []=fLabelDffIncludedTrialRange(dff,filepath,ind_1stFrame,trial2include)
if exist(filepath,'file')
    return;
else
    figDff=figure;
    plot(dff,'k-');hold on;
    yrange=get(gca,'Ylim');
    for i=1:size(trial2include,1)
        tempx=[ind_1stFrame(trial2include(i,1)),ind_1stFrame(trial2include(i,2))];
        xpatch=[tempx,fliplr(tempx)];
        ypatch=[yrange(1),yrange(1),yrange(2),yrange(2)];
        patch(xpatch,ypatch,'b','FaceAlpha',0.5);
    end
    saveas(figDff,filepath,'jpg');
end
end
