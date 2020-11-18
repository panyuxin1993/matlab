function [ session_name ] = fPSTHandEndingPointTrajectoryOneSession( rootpath,varargin )
%FPSTHANDENDINGPOINTTRAJECTORY plot PSTH and trajectories of ending point
%of each trials. The point is to see whether these two lines overlap
%   input: rootpath, the folder of session to plot
%   output: session name of plotted session

%%
%clear
% close all
% clc
if ~exist('rootpath','var')
    rootpath='F:\FP\0';
end
if ~isempty(varargin)
    position=varargin{1};
else
    position={[0,700,700,200],[700,700,700,200]};
end
temp=strsplit(rootpath,'\');
session_name=temp(end);
SavingFolder=rootpath;
DatePath=[];
DatePath{1}=rootpath;
files = dir(strcat(rootpath,'\*Virables.mat'));
if length(files)==1
    behaviorFile=files(1).name;
else
    warning('Multiple beh files');
end
%%
% the session to be plotted
for n_path=1:length(DatePath)
    clearvars -except n_path SavingFolder DatePath behaviorFile position
    cd(DatePath{n_path});
    %
    FrameInfo = dir('*.log');
%         [TrialCount,TrialStart_FrameCount]=textread(FrameInfo.name,'%d %d','headerlines',12);
    fileID = fopen(FrameInfo.name);
    C=textscan(fileID,'%d %d','HeaderLines',16);
    fclose(fileID);
    TrialCount=C{1,1};
    TrialStart_FrameCount=C{1,2};
    nTrial = length(TrialCount);
%     ImagingSetup
    ImagingSetup=FrameInfo.name(1:end-4);
    disp(ImagingSetup);
    ImagingSetup(ImagingSetup=='_')='-';
    ImagingSetup_L=sprintf('%s_Left', ImagingSetup);
    ImagingSetup_R=sprintf('%s_Right', ImagingSetup);
    ImagingSetup_Ctrl=sprintf('%s_410Ctrl', ImagingSetup);
    
    Session470LEDstartFrame=1;% usually 205
    Session410LEDstartFrame=2;

    fiberstr={'Left','Right'};
    FrameRate=40;
    FrameTime=1000/FrameRate;
    %% compute dff or just load
    if exist('dff_temp.mat','file')
        load('dff_temp');
    else
        dff=fGetFPdff(SavingFolder);
    end
    %% load Beh mat data and extract behavior events
    load(behaviorFile);%load behavior data
    [trialType,behrule] = fGetTrialType( Data_extract);%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
    frT=FrameTime*2;%2 channel, so framerate and frametime should be half
    frameNumTime=[1,2.5];%from 2s before align point to 5s after align point
    frameNum=double(round(frameNumTime*1000/frT));
    [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime(Data_extract,round(double(TrialStart_FrameCount')/2),frT);
    %% for each dff, plot mean trace, color plot, etc.
    dffstr={'dff baseline correction first','dff motion correction first','dff470','dff410'};
    behEventAlign='stimOnset';%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
    behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
    masklick='yes';%'yes' if mask lickings or 'no' when want to see all activity
    for nfiber=1:2
        for ndff=[1] %plot for each dff, see which is better and test whether 410 signals have similar trend
            if strcmp(behEventAlign,'stimOnset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
                [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff{1,ndff}(nfiber,:), behEventFrameIndex,  frameNum );          
            else
                [ dff_aligned, behEvent_aligned, licking_aligned ] = fAlignSigBehEvent( dff{1,ndff}(nfiber,:), behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
            end

            figMeanTrace=figure;%plot mean trace
            set(gcf, 'position', position{nfiber});
            titlestr={'Correct','Error','Miss','Violation'};
%             titlestr=strcat(ImagingSetup,'-',fiberstr{nfiber},'-',titlestr);
            color_mean_trace=fMeanTraceColor(behrule,size(trialType,2));
            for nStim=1:size(trialType,2) %for each stimulus%[1,6]%just 2 end trials 
                for  nResult=1:size(trialType,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
                    selectedTrialInd=trialType(nResult,nStim,:);
                    selectedTrialInd=logical(squeeze(selectedTrialInd))';
%                     sot=behEvent_aligned.stimOnset(selectedTrialInd);% stim onset, white
%                     flt=behEvent_aligned.lickFirst(selectedTrialInd);% first lick, black
%                     flt_l=behEvent_aligned.lickFirst_left(selectedTrialInd);%left lick, magenta dots
%                     flt_r=behEvent_aligned.lickFirst_right(selectedTrialInd);%right lick, red dots
%                     rwt=behEvent_aligned.rewTime(selectedTrialInd);%reward, cyan dots
%                     go=behEvent_aligned.go(selectedTrialInd);%go cue,white
                    neuralActivity=dff_aligned(selectedTrialInd,:);
                    temp=find(selectedTrialInd);
%                     leftLick=cell(1,length(temp));
%                     rightLick=cell(1,length(temp));
%                     for i=1:length(temp)
%                         leftLick{i}= licking_aligned.leftLick{temp(i)};
%                         rightLick{i}=licking_aligned.rightLick{temp(i)};
%                     end

                    %plot mean trace
                    figure(figMeanTrace);%save mean trace
                    subplot(1,size(trialType,1),nResult);
                    [neuralActivityMean, neuralActivityCI]=fMean_CI(neuralActivity,0.05);
                    %shadow as ci
                    %{
                    xpatch=[1:size(neuralActivity,2), fliplr(1:size(neuralActivity,2))];
                    ypatch=[neuralActivityCI(1,:),fliplr(neuralActivityCI(2,:))];
                    p=patch(xpatch,ypatch,color_mean_trace{nStim});%plot confidence interval
                    p.FaceAlpha=0.1;
                    p.EdgeColor=color_mean_trace{nStim};%'none';
                    hold on;
                    %}
                    curve_meanTrace(nStim)=plot(1:size(neuralActivity,2),neuralActivityMean,'Color',color_mean_trace{nStim},'linewidth',2);
                    hold on;
                    title(titlestr{nResult});
                    y_lim=get(gca,'Ylim');
                    plot([frameNumTime(1)*round(1000/frT),frameNumTime(1)*round(1000/frT)],[y_lim(1),y_lim(2)],'k-');%align to a behavior event
                    xlabel(['time(s) from ',behEventAlign]);
                    %plot endpoint dff
                    [ ~,~,neuralActivityGroupedtemp ] = fGetEndDff( neuralActivity ,'smooth','bin',0.3*1000/frT);
                    [ timepoint,endDff,neuralActivityGrouped ] = fGetEndDff( neuralActivityGroupedtemp ,'raw');
                    for i=1:size(neuralActivityGrouped,1)
                        plot(neuralActivityGrouped(i,:),'Color',color_mean_trace{nStim},'linewidth',1.5);
                    end
                    scatter(timepoint,endDff,50,color_mean_trace{nStim},'MarkerFaceColor','w');
                    if nResult==1
                        ylabel('\it\DeltaF/F');
                    end
                    set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
                    x_lim=get(gca,'Xlim');
                    xlim([x_lim(1),size(neuralActivity,2)]);
%                     ylim([nanmin(nanmin(neuralActivityGrouped)),nanmax(nanmax(neuralActivityGrouped))]);
%                     set(gca,'FontName','Arial','FontSize',14);    
                    box off;
                end
            end
            figure(figMeanTrace);
            subplot(1,size(trialType,1),1);
            h=legend(curve_meanTrace(:),num2str(Data_extract.Stimuli(:)),'Location','best');
            set(h,'box','off');
            if strcmp(masklick,'yes')
                suptitle(strcat(ImagingSetup,'-',fiberstr{nfiber},'-',dffstr{ndff},'-masklick'));
            else
                suptitle(strcat(ImagingSetup,'-',fiberstr{nfiber},'-',dffstr{ndff}));
            end
            saveas(figMeanTrace,[ImagingSetup,'-',fiberstr{nfiber},dffstr{ndff},'-algin to ',behEventAlign,'-sort ',behEventSort,'-masklick',masklick,'mean trace with end point dff.fig'],'fig');
%             close all;
        end
    end
    
end



end
