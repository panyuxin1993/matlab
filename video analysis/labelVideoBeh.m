%note: must be run in matlab2018b
filepath='F:\video tracking\M2 imaging video';
savepath=[filepath,filesep,'summary figures'];
summaryFile=[filepath,filesep,'imaging_video_data_summary.xlsx'];
[~,~,temp]=xlsread(summaryFile,1);
dataSummaryT=cell2table(temp(2:end,:));
dataSummaryT.Properties.VariableNames =temp(1,:);
nSession=size(dataSummaryT,1);
%choosing right data or exluded trash data
[~,datasource,~]=xlsread(summaryFile,2);%%%%%%%%%%%%%%%choosing which DLC model
%parameters for trial type
i_selectivity=4;%*********variable**************
selectivitystr={'stimuli', 'difficulty','sensory','choice'};%sensory means grouping difficulties;
trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
outcomeType='combineCorErr';%'combineCorErr','divideCorErr',
%parameters for time bin
nonOverlappingBin=1;
fr=24/nonOverlappingBin;
treshold4likelihood=0.1;
bodyparts={'Tongue','Nose'};
coordinates={'x','x'};
for iSession=7:nSession%can be a loop
    led_file=[filepath,filesep,dataSummaryT.OLEDFileName{iSession},'.csv'];
    behdata=[filepath,filesep,dataSummaryT.session{iSession},'_Virables.mat'];
    
    %% from .beh file calculated a vector of trial start, for trial start from video use
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
    frameTrialStartVideo=floor(frameTrialStartVideo/nonOverlappingBin);
    ind_tr_1=1;%using data from trial 1
    [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( Data_extract, frameTrialStartVideo, 1000/fr ,ind_tr_1);%get behavior event time
    [trialType,rule,trialTypeStr] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix','left',outcomeType);%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
    
    nbodyparts=length(bodyparts);
    TcoCell=cell(size(bodyparts));
    for iBodyPart=1:nbodyparts
        % import DLC result
        col_datasource=cellfun(@(x) strcmp(bodyparts{iBodyPart},x),datasource(1,:));
        row_datasource=cellfun(@(x) strcmp(dataSummaryT.session{iSession},x),datasource(:,1));
        DLCiteration=datasource{row_datasource,col_datasource};
        switch DLCiteration
            case 'iteration-1'
                file_trace=[filepath,filesep,'iteration1',filesep,dataSummaryT.DLCFileName1{iSession},'.csv'];
                flagSetNAN=0;
            case 'iteration-2'
                file_trace=[filepath,filesep,'iteration2',filesep,dataSummaryT.DLCFileName2{iSession},'.csv'];
                flagSetNAN=0;
            case 'iteration-3'
                file_trace=[filepath,filesep,'iteration3',filesep,dataSummaryT.DLCFileName3{iSession},'.csv'];
                flagSetNAN=0;
            case 'iteration-4'
                file_trace=[filepath,filesep,'iteration4',filesep,dataSummaryT.DLCFileName4{iSession},'.csv'];
                flagSetNAN=0;
            case 'iteration-5'
                file_trace=[filepath,filesep,'iteration5',filesep,dataSummaryT.DLCFileName5{iSession},'.csv'];
                flagSetNAN=0;
            otherwise %no data, skip this loop
                file_trace=[filepath,filesep,'iteration2',filesep,dataSummaryT.DLCFileName2{iSession},'.csv'];%for  iteration-1~3
                %                 file_trace=[filepath,filesep,'iteration4',filesep,dataSummaryT.DLCFileName4{iSession},'.csv'];%for  iteration-4+
                
                flagSetNAN=1;
        end
        name_file_trace=strsplit(file_trace,'.');
        if ~exist([name_file_trace{1},'.mat'],'file')
            [dcnum,dctxt,~] = xlsread(file_trace,1);%this process is time consuming, so not do very often
            save([name_file_trace{1},'.mat'],'dcnum','dctxt');
        else
            load([name_file_trace{1},'.mat'])
        end
        file_video=[name_file_trace{1},'_labeled.mp4'];
        %% get coordinates of body parts, analogy to dff
        indcol1=cellfun(@(x) strcmp(bodyparts{iBodyPart},x),dctxt(2,:));
        indcol2=cellfun(@(x) strcmp(coordinates{iBodyPart},x),dctxt(3,:));
        indcol_likelihood=cellfun(@(x) strcmp('likelihood',x),dctxt(3,:));
        indcol_co=find(indcol1.*indcol2);
        indcol_like=find(indcol1.*indcol_likelihood);
        bodyco=dcnum(:,indcol_co);
        bodycoli=dcnum(:,indcol_like);
        bodyco(bodycoli<treshold4likelihood)=nan;%rule out those low likelihood data
        if ~strcmp(coordinates{iBodyPart},'likelihood')
            if contains(bodyparts{iBodyPart},'LickPort')
                bodyco=fBaselineCorrection(bodyco,5*fr);%40s as span
            end
            bodyco=bodyco-nanmean(bodyco);%calculate pixel shift, if it is likelihood, no need for normalization
        end
        %bin tongue data with 2 frames, non-overlapping
        bodyco=bodyco(1:floor(length(bodyco)/nonOverlappingBin)*nonOverlappingBin,1);
        bodyco=reshape(bodyco,nonOverlappingBin,[]);
        bodyco=nanmean(bodyco,1);
        %         %replace nan as zeros
        %         bodyco(isnan(bodyco))=0;
        %remove whole session if low likelihood
        if flagSetNAN==1
            bodyco(:)=nan;
        end
    end
    nfigcol=size(trialType,1)-2;%4 column(correct/error/miss/violation),companied with 4 lick raster
    file_video=[name_file_trace{1},'_labeled.mp4'];
    ReaderObj = VideoReader(file_video);
    for nStim=1:size(trialType,2) %for each stimulus
        %open a new video to write
        file_videoOut=[name_file_trace{1},'_labeled_beh_longDelay-stim',num2str(nStim),'.avi'];
        WriterObj=VideoWriter(file_videoOut);%待合成的视频(不仅限于avi格式)的文件路径
        WriterObj.FrameRate=ReaderObj.FrameRate;
        open(WriterObj);
        for  nResult=1:nfigcol
            selectedTrialInd=trialType(nResult,nStim,:);
            selectedTrialInd=logical(squeeze(selectedTrialInd))';
            longDelayInd=Data_extract.Delay_duration>800;%select long delay trials (>800ms) and only label video made of these trials
            selectedTrialInd=logical(selectedTrialInd.*longDelayInd);
            listSelectedTrialInd=find(selectedTrialInd);
            for iTrial=1:length(listSelectedTrialInd)
                indFrame=frameTrialStartVideo(listSelectedTrialInd(iTrial));
                if isnan(indFrame)
                    continue;
                else
                    ReaderObj.CurrentTime=(indFrame-1)/ReaderObj.FrameRate;%The 1st frame in a video is 0
                    %                 i=ReaderObj.CurrentTime*ReaderObj.FrameRate+1;
                    while hasFrame(ReaderObj) && indFrame<frameTrialStartVideo(min(listSelectedTrialInd(iTrial)+1,length(frameTrialStartVideo))) %frome current trial to the next and not the end of the session
                        frame = readFrame(ReaderObj);  % 读取每一帧
                        imshow(frame);                 % 显示每一帧
                        temp=find(frameTrialStartVideo>=indFrame);
                        if ~isempty(temp)
                            indTrial=temp(1)-1;
                        else
                            indTrial=length(frameTrialStartVideo);
                        end
                        text(0.8,0.95,['Trial-',num2str(indTrial)],'FontName','Arial','FontSize',20,'Unit','normalized','Color','y');
                        if sum((behEventFrameIndex.stimOffset<=indFrame).*(behEventFrameIndex.go>=indFrame))>0 %in delay of one trial
                            text(0.8,0.9,['delay',num2str((behEventFrameIndex.go(indTrial)-behEventFrameIndex.stimOffset(indTrial))/fr)],'FontName','Arial','FontSize',20,'Unit','normalized','Color','y');
                        end
                        frameNew=getframe(gca);
                        writeVideo(WriterObj,frameNew);
                        indFrame=indFrame+1;
                    end
                end
            end
        end
        close(WriterObj);
    end
end

%% auxiliary fuction
function [timeTrialStart]=fTrialStartStr2double(strTrialStart)
%transform the trial start time from string form to double form
t=strsplit(strTrialStart,'_');
t=str2double(t);
timeTrialStart=t(4)*3600+t(5)*60+t(6)+t(7)/1000; %time(s) relative to trial1
end