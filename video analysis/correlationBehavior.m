dbstop if error;

close all;
clear;
filepath='F:\video tracking\M2 imaging video';
savepath=[filepath,filesep,'summary figures'];
summaryFile=[filepath,filesep,'imaging_data_summary.xlsx'];
[~,~,temp]=xlsread(summaryFile,1);
dataSummaryT=cell2table(temp(2:end,:));
dataSummaryT.Properties.VariableNames =temp(1,:);
DLCiteration=1;

for iSession=7%can be a loop
led_file=[filepath,filesep,dataSummaryT.OLEDFileName{iSession},'.csv'];
behdata=[filepath,filesep,dataSummaryT.session{iSession},'_Virables.mat'];
switch DLCiteration
    case 1
        file_trace=[filepath,filesep,'iteration1',filesep,dataSummaryT.DLCFileName1{iSession},'.csv'];
    case 2
        file_trace=[filepath,filesep,'iteration2',filesep,dataSummaryT.DLCFileName2{iSession},'.csv'];
end

%导入数据
name_file_trace=strsplit(file_trace,'.');
if ~exist([name_file_trace{1},'.mat'],'file')
    [dcnum,dctxt,~] = xlsread(file_trace,1);%this process is time consuming, so not do very often
    save([name_file_trace{1},'.mat'],'dcnum','dctxt');
else
    load([name_file_trace{1},'.mat'])
end

%% from .beh file calculated a vector of trial start, for trial start from video use
fr=24;
load(behdata);
timeTrialStartBeh=cellfun(@fTrialStartStr2double,Data_extract.Time_trialStart);% time(s) relative to start of session day
name_file_led=strsplit(led_file,filesep);
name_file_led2=strsplit(name_file_led{end},'.');
if ~exist([filepath,filesep,name_file_led2{1},'.mat'],'file')
    T=readtable(led_file);%different sessions may need mannually set parameters and save result separately
    frameTrialStartVideo=fTrialStartByOLEDrefArduino( T, fr ,timeTrialStartBeh);
    save([filepath,filesep,name_file_led2{1},'.mat'],'frameTrialStartVideo');
else
    load([filepath,filesep,name_file_led2{1},'.mat'])
end
frameTrialStartVideo=reshape(frameTrialStartVideo,1,[]);
% figure;
% plot(frameTrialStartVideo);

%
%% global settings
%coordinates of which body part
bodyparts={'Tongue','Tongue'};%{'Tongue','LeftHandFingerTip','RightHandFingerTip','LeftHandFingerRoot','RightHandFingerRoot','Nose','LeftWhiskerTip','RightWhiskerTip','LeftLickPort','RightLickPort'};%{'Tongue','Tongue','LeftHandFingerTip','LeftHandFingerTip','RightHandFingerTip','RightHandFingerTip','Nose','Nose','LeftWhiskerTip','LeftWhiskerTip','RightWhiskerTip','RightWhiskerTip'};%{Tongue,LeftHandFingerTip,LeftHandFingerRoot,LeftLickPort,Nose,LeftWhiskerTip,LeftWhiskerRoot, etc.}
coordinates={'x','y'};%{'likelihood','likelihood','likelihood','likelihood','likelihood','likelihood','likelihood','likelihood','likelihood','likelihood'};%{'x','y','x','y','x','y','x','y','x','y','x','y'};%{x,y,likelihood};
nbodyparts=2;%length(bodyparts);
%align to behavior event
behEventAlign='first lick';%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
masklick='no';
i_selectivity=4;%*********variable**************
selectivitystr={'stimuli','difficulty','sensory','choice'};%sensory means grouping difficulties;
trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
outcomeType='combineCorErr';%'combineCorErr','divideCorErr',
%plotting settings
if strcmp(behEventAlign,'stimOnset')
    frameNumTime=[0.5,2];%from 5s before align point to 5s after align point
else
    frameNumTime=[1,2];%from 5s before align point to 5s after align point
end

%% get behavior event from .beh
ind_tr_1=1;%using data from trial 1
frameNum=double(round(frameNumTime*fr));
[behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( Data_extract, frameTrialStartVideo, 1000/fr ,ind_tr_1);%get behavior event time
[trialType,rule,trialTypeStr] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix','left',outcomeType);%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
%% settings for plot
% yrange=[ 0 , 1 ];
%plot color plot,4 column(correct/error/miss/violation)
nfigcol=size(trialType,1)-2;%4 column(correct/error/miss/violation),companied with 4 lick raster
A=figure;
set(gcf, 'position', [0 0 250*nfigcol 300]);
figMeanTrace=figure;%plot mean trace
%set(gcf, 'position', [0 0 250*nfigcol 200*nbodyparts]);
set(gcf, 'PaperPosition', [1 1 2.5*nfigcol 2*nbodyparts]);
for iBodyPart=1:nbodyparts
    %% get coordinates of body parts, analogy to dff
    indcol1=cellfun(@(x) strcmp(bodyparts{iBodyPart},x),dctxt(2,:));
    indcol2=cellfun(@(x) strcmp(coordinates{iBodyPart},x),dctxt(3,:));
    indcol_likelihood=cellfun(@(x) strcmp('likelihood',x),dctxt(3,:));
    indcol_co=find(indcol1.*indcol2);
    indcol_like=find(indcol1.*indcol_likelihood);
    bodyco=dcnum(:,indcol_co);
    bodycoli=dcnum(:,indcol_like);
%     bodyco(bodycoli<0.95)=nan;%rule out those low likelihood data
    bodyco=reshape(bodyco,1,[]);
    bodyco=bodyco-nanmean(bodyco);%calculate pixel shift
    %% aligned the coordinates change of body parts to behavior events
    if strcmp(behEventAlign,'stimOnset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
        [ bodyco_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( bodyco, behEventFrameIndex,  frameNum );
    else
        [ bodyco_aligned, behEvent_aligned, licking_aligned ] = fAlignSigBehEvent( bodyco, behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
    end

    %% plot coordinates change, ref plotDffPSTH
    
    ColLimit = prctile(bodyco,90);
    if strcmp(outcomeType,'divideCorErr')
        titlestr={'Correct','Error','Miss','Violation'};
    else
        titlestr={'Do','Miss','Violation'};
    end
    titlestr=strcat(bodyparts{iBodyPart},'-',coordinates{iBodyPart},'-',titlestr);
    if size(trialType,2)==6
        color_mean_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};
    elseif size(trialType,2)==2
        color_mean_trace={[0 0 1],[1 0 0]};
    end
    if strcmp(rule,'low click rate-right')
        color_mean_trace=fliplr(color_mean_trace);
    end
    for nStim=1:size(trialType,2) %for each stimulus
        for  nResult=1:nfigcol
            selectedTrialInd=trialType(nResult,nStim,:);
            selectedTrialInd=logical(squeeze(selectedTrialInd))';
            sot=behEvent_aligned.stimOnset(selectedTrialInd);% stim onset, white
            flt=behEvent_aligned.lickFirst(selectedTrialInd);% first lick, black
            flt_l=behEvent_aligned.lickFirst_left(selectedTrialInd);%left lick, magenta dots
            flt_r=behEvent_aligned.lickFirst_right(selectedTrialInd);%right lick, red dots
            rwt=behEvent_aligned.rewTime(selectedTrialInd);%reward, cyan dots
            go=behEvent_aligned.go(selectedTrialInd);%go cue,white
            neuralActivity=bodyco_aligned(selectedTrialInd,:);
            temp=find(selectedTrialInd);
            leftLick=cell(1,length(temp));
            rightLick=cell(1,length(temp));
            for i=1:length(temp)
                leftLick{i}= licking_aligned.leftLick{temp(i)};
                rightLick{i}=licking_aligned.rightLick{temp(i)};
            end
            %decide which behavior event to be sort
            if strcmp(behEventSort,'first lick')
                [B,I]=sort(flt);
            elseif strcmp(behEventSort,'reward')
                [B,I]=sort(rwt);
            elseif strcmp(behEventSort,'go cue')
                [B,I]=sort(go);
            end
            if iBodyPart==1
                go=go(I);
                sot=sot(I);
                flt=flt(I);
                flt_l=flt_l(I);
                flt_r=flt_r(I);
                rwt=rwt(I);
                figure(A);
                subplot(size(trialType,2),2*nfigcol,2*nResult-1+2*nfigcol*(nStim-1));
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
                    set(gca,'xtick',[round(fr*(frameNumTime(1)-floor(frameNumTime(1))))+1:round(fr):size(neuralActivity,2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
                    xlabel(['time(s) from ',behEventAlign],'FontName','Arial','FontSize',14);
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
                    ylabel([num2str(Data_extract.Stimuli(nStim))]); %' clicks/s'
                end
                xlim(x_lim);
                
                
                
                subplot(size(trialType,2),2*nfigcol,2*nResult+2*nfigcol*(nStim-1));
                %             subplot(size(trialType,2),nfigcol,nResult+nfigcol*(nStim-1));%only plot licking raster
                for i=1:length(leftLick)%each trial as a row
                    for jl=1:length(leftLick{I(i)})
                        line([leftLick{I(i)}(jl),leftLick{I(i)}(jl)],[i-0.5,i+0.5],'color','b','linewidth',1);%left lick
                        hold on;
                    end
                    for jr=1:length(rightLick{I(i)})
                        line([rightLick{I(i)}(jr),rightLick{I(i)}(jr)],[i-0.5,i+0.5],'color','r','linewidth',1);%right lick
                        hold on;
                    end
                end
                set(gca,'ytick',sum(selectedTrialInd),'yticklabel',sum(selectedTrialInd),'ydir','normal');
                %         set(gca,'ytick',[]);%licking raster have corresponding color plot label, so not need
                if nStim==size(trialType,2)
                    set(gca,'xtick',[round(fr*(frameNumTime(1)-floor(frameNumTime(1))))+1:round(fr):size(neuralActivity,2)],'xticklabel',ceil(-frameNumTime(1)):1:frameNumTime(2));
                    xlabel(['time(s) from ',behEventAlign],'FontName','Arial','FontSize',12);
                else
                    set(gca,'xtick',[]);
                end
                xlim(x_lim);
                if sum(selectedTrialInd)>0
                    ylim([0.5 sum(selectedTrialInd)+0.5]);
                end
            end
            
            %plot mean trace
            figure(figMeanTrace);%save mean trace
            %subplot(size(trialType,2)+3,2*nfigcol,2*nResult-1+2*nfigcol*(size(trialType,2)+2));
            subplot(nbodyparts,nfigcol,nResult+(iBodyPart-1)*nfigcol);
            [neuralActivityMean, neuralActivityCI]=fMean_CI(neuralActivity,0.05);
            %plot CI
            %
            if size(trialType,2)==2%if only two stimuli, plot CI, otherwise, do not plot for the reason of simplicity
                xpatch=[1:size(neuralActivity,2), fliplr(1:size(neuralActivity,2))];
                ypatch=[neuralActivityCI(1,:),fliplr(neuralActivityCI(2,:))];
                p=patch(xpatch,ypatch,color_mean_trace{nStim});%plot confidence interval
                p.FaceAlpha=0.1;
                p.EdgeColor=color_mean_trace{nStim};%'none';
                hold on;
            end
            %}
            %plot mean
            curve_meanTrace(nStim)=plot(1:size(neuralActivity,2),neuralActivityMean,'Color',color_mean_trace{nStim},'linewidth',2);
            hold on;
        end
    end
    for nResult=1:nfigcol
        subplot(nbodyparts,nfigcol,nResult+(iBodyPart-1)*nfigcol);
        title(titlestr{nResult});
        if ~exist('yrange','var')
            y_lim=get(gca,'Ylim');
        else
            y_lim=yrange;
        end
        xlim(x_lim);
        plot([round(frameNumTime(1)*fr),round(frameNumTime(1)*fr)]+1,[y_lim(1),y_lim(2)],'k-');%align to a behavior event
        hold on;
        if strcmp(behEventAlign,'stimOnset')
            plot(((round(frameNumTime(1)+0.5)*fr)+1)*ones(1,2),[y_lim(1),y_lim(2)],'k-');%align to a behavior event
        end
        if iBodyPart==nbodyparts
            xlabel(['time(s) from ',behEventAlign]);
        end
        if nResult==1
            if strcmp(coordinates{iBodyPart},'likelihood')
                ylabel('likelihood');
            else
                ylabel('\it\DeltaC');
            end
        end
        set(gca,'xtick',[round(fr*(frameNumTime(1)-floor(frameNumTime(1))))+1:round(fr):size(neuralActivity,2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
%         set(gca,'xtick',[round(fr*(frameNumTime(1)-floor(frameNumTime(1))))+1:round(fr):size(neuralActivity,2)],'xticklabel',0:1:frameNumTime(2));
        set(gca,'FontName','Arial','FontSize',12);
        set(gca,'Ylim',y_lim);
    end

end
figure(figMeanTrace);
subplot(nbodyparts,nfigcol,1);
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
saveas(figMeanTrace,[savepath,filesep,'likelihood-algin to ',behEventAlign,'-sort ',behEventSort,'mean trace.pdf'],'pdf');
%}
end

%% auxiliary fuction
function [timeTrialStart]=fTrialStartStr2double(strTrialStart)
%transform the trial start time from string form to double form
t=strsplit(strTrialStart,'_');
t=str2double(t);
timeTrialStart=t(4)*3600+t(5)*60+t(6)+t(7)/1000; %time(s) relative to trial1
end