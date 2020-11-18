%to examine the raw likelihood and coordinates shifts
dbstop if error;
close all;
clear;
filepath='F:\video tracking\M2 imaging video';
savepath=[filepath,filesep,'summary figures'];
summaryFile=[filepath,filesep,'imaging_video_data_summary.xlsx'];
[~,~,temp]=xlsread(summaryFile,1);
dataSummaryT=cell2table(temp(2:end,:));
dataSummaryT.Properties.VariableNames =temp(1,:);
DLCiteration=1:5;%DLC model from which iteration
fr=24;%frame rate of the video

%coordinates of which body part
bodyparts={'Tongue','LeftHandFingerTip','RightHandFingerTip','LeftHandFingerRoot','RightHandFingerRoot','Nose','LeftWhiskerTip','RightWhiskerTip','LeftLickPort','RightLickPort'};%{'Tongue','Tongue','LeftHandFingerTip','LeftHandFingerTip','RightHandFingerTip','RightHandFingerTip','Nose','Nose','LeftWhiskerTip','LeftWhiskerTip','RightWhiskerTip','RightWhiskerTip'};%{Tongue,LeftHandFingerTip,LeftHandFingerRoot,LeftLickPort,Nose,LeftWhiskerTip,LeftWhiskerRoot, etc.}

coordinates={'x','y'};%{'likelihood','likelihood','likelihood','likelihood','likelihood','likelihood','likelihood','likelihood','likelihood','likelihood'};%{'x','y','x','y','x','y','x','y','x','y','x','y'};%{x,y,likelihood};
colorCo={[1,0,0],[0,0,1],[0.3,0.3,0.3]};%colorCo{3} for low likelihood dots
nbodyparts=2;%length(bodyparts);
threshold4likelihood=0.1;
for i_iter=5
for iBody=6%1:length(bodyparts)
    for iSession=1:size(dataSummaryT,1)%can be a loop
        switch DLCiteration(i_iter)
            case 1
                file_trace=[filepath,filesep,'iteration1',filesep,dataSummaryT.DLCFileName1{iSession},'.csv'];
            case 2
                file_trace=[filepath,filesep,'iteration2',filesep,dataSummaryT.DLCFileName2{iSession},'.csv'];
            case 3
                file_trace=[filepath,filesep,'iteration3',filesep,dataSummaryT.DLCFileName3{iSession},'.csv'];
            case 4
                file_trace=[filepath,filesep,'iteration4',filesep,dataSummaryT.DLCFileName4{iSession},'.csv'];
            case 5
                file_trace=[filepath,filesep,'iteration5',filesep,dataSummaryT.DLCFileName5{iSession},'.csv'];
        end
        name_file_trace=strsplit(file_trace,'.');
        if ~exist([name_file_trace{1},'.mat'],'file')
            [dcnum,dctxt,~] = xlsread(file_trace,1);%this process is time consuming, so not do very often
            save([name_file_trace{1},'.mat'],'dcnum','dctxt');
        else
            load([name_file_trace{1},'.mat'])
        end

        indcolbody=cellfun(@(x) strcmp(bodyparts{iBody},x),dctxt(2,:));
        indcol_likelihood=cellfun(@(x) strcmp('likelihood',x),dctxt(3,:));
        indcol_like=find(indcolbody.*indcol_likelihood);
        bodycoli=dcnum(:,indcol_like);
          figCoLi=figure;
          figHisto=figure;
          coorCell=cell(2,2); %1d {x,y},2d high/low-likelihood
        % get coordinates of body parts, analogy to dff
        for iCoordinate=1:length(coordinates)
            figure(figCoLi);
            subh1=subplot(2,2,1);%plot coordinates shift
            indcolcoor=cellfun(@(x) strcmp(coordinates{iCoordinate},x),dctxt(3,:));
            indcol_co=find(indcolbody.*indcolcoor);
            bodyco=dcnum(:,indcol_co);
            bodyco_high=bodyco;
            bodyco_low=bodyco;
            bodyco_high(bodycoli<threshold4likelihood)=nan;%rule out those low likelihood data
            bodyco_low(bodycoli>=threshold4likelihood)=nan;
            coorCell{iCoordinate,1}=bodyco_high;
            coorCell{iCoordinate,2}=bodyco_low;
            span=10*fr;%20s as a span
%             bodyco_high=fBaselineCorrection(bodyco_high,span);
            bodyco_baseline=nanmean(bodyco_high);
            bodyco_low=bodyco_low-bodyco_baseline;%calculate pixel shift
            bodyco_high=bodyco_high-bodyco_baseline;%calculate pixel shift
            ts=(1:length(bodyco))/fr;
            ts_high=ts;
            ts_low=ts;
            ts_high(bodycoli<threshold4likelihood)=nan;
            ts_low(bodycoli>=threshold4likelihood)=nan;
            curve(iCoordinate)=plot(ts_high,bodyco_high,'color',colorCo{iCoordinate});hold on;
            plot(ts_low,bodyco_low,'color',colorCo{3});hold on;
            xlabel('time(s)');
            set(gca,'Box','off');
            subh2=subplot(2,2,2);%histogram
            histogram(bodyco_high,'FaceColor',colorCo{iCoordinate},'EdgeColor',colorCo{iCoordinate},'Orientation','horizontal');hold on;
            histogram(bodyco_low,'FaceColor',colorCo{3},'EdgeColor',colorCo{3},'Orientation','horizontal');
            set(gca,'Box','off');
        end
        figure(figHisto);
%         scatter(coorCell{1,1},coorCell{2,1},'r');hold on;
%         scatter(coorCell{1,2},coorCell{2,2},'k');
        
%         subplot(1,2,1);
        h_high=histogram2(coorCell{1,1},coorCell{2,1},'DisplayStyle','tile','ShowEmptyBins','on');
       
%         set(gca,'Xlim',[0,640],'Ylim',[0,480]);
%         subplot(1,2,1);
%         h_low=histogram2(coorCell{1,2},coorCell{2,2},'DisplayStyle','tile','ShowEmptyBins','on');
%         set(gca,'Xlim',[0,640],'Ylim',[0,480]);
        figure(figCoLi);
        subplot(2,2,1);%plot coordinates shift
        title(bodyparts{iBody});
        ylabel('\it\Delta C');
        legend(curve,'x','y');
        
        figure(figCoLi);
        subh3=subplot(2,2,3);%plot likelihood
        plot(ts,bodycoli);
        set(gca,'Box','off');
        subh4=subplot(2,2,4);%plot likelihood
        histogram(bodycoli,'BinWidth',0.1,'Orientation','horizontal');
        set(gca,'Ylim',[0,1]);
        ylabel('likelihood');
        set(gca,'Box','off');
        
        set(subh3,'Position',[0.1,0.1,0.75,0.2]);
        set(subh4,'Position',[0.87,0.1,0.05,0.2]);
        set(subh1,'Position',[0.1,0.4,0.75,0.5]);
        set(subh2,'Position',[0.87,0.4,0.05,0.5]);
        suptitle(strcat(dataSummaryT.videoName(iSession),'iteration-',num2str(i_iter)));
    end
end
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