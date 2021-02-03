function [ session_name ] = fSummaryFPOneSession( rootpath )
%FSUMMARYFPONESESSION Summary of this function goes here

%FPLOTFPEXAMPLE plot example site from FP recording data
% color plot and meantrace are in same fig

% position={[0.1,0.1,4,6],[4.1,0.1,4,6];[8.1,0.1,4,6],[12.1,0.1,4,6]};
% legposition={[0.91/4 0.1 0.02 0.08],[0.91/2 0.1 0.02 0.08];[0.91*3/4 0.1 0.02 0.08],[0.91 0.1 0.02 0.08]};
% close all;
position={[0,50,400,600],[400,50,400,600];[800,50,400,600],[1200,50,400,600]};
positionBeh=[1600,100,200,200];
temp=strsplit(rootpath,'\');
session_name=temp(end);
cd(rootpath);
SavingFolder=rootpath;
behaviorFile=dir([rootpath,'\*Virables.mat']);
FrameInfo = dir([rootpath,'\*.log']);
fileID = fopen([rootpath,'\',FrameInfo.name]);
C=textscan(fileID,'%d %d','HeaderLines',16);
%     ImagingSetup
ImagingSetup=FrameInfo.name(1:end-4);
disp(ImagingSetup);
ImagingSetup(ImagingSetup=='_')='-';
fclose(fileID);
TrialCount=C{1,1};
TrialStart_FrameCount=C{1,2};
fiberstr={'Left','Right'};

%% compute dff or just load
if exist('dff_temp.mat','file')
    load('dff_temp');
else
    dff=fGetFPdff(SavingFolder);
end

%% load Beh mat data and extract behavior events
load(behaviorFile.name);%load behavior data
[trialType,behrule] = fGetTrialType( Data_extract);%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
trialType(3:4,:,:)=[];%do not plot vio and miss
FrameRate=40;
FrameTime=1000/FrameRate;
frT=FrameTime*2;%2 channel, so framerate and frametime should be half
frameNumTime=[1,3];%from 2s before align point to 5s after align point
frameNum=double(round(frameNumTime*1000/frT));
[behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime(Data_extract,round(double(TrialStart_FrameCount')/2),frT);
%% compare different dff('dff baseline correction first','dff motion correction first','dff470','dff410')
%for each dff, plot mean trace, color plot, etc.
dffstr={'dff baseline correction first','dff motion correction first','dff470','dff410'};
behEventSort='go cue';% string can be in{'first lick','reward','go cue'};

for nevent=1:2 
    behEventAlign={'delay onset','first lick'};%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
    masklick={'yes','no'};%'yes' if mask lickings or 'no' when want to see all activity
    behEventAlign=behEventAlign{nevent};
    masklick=masklick{nevent};
    for nfiber=1:2
        %         A=figure;
        %         set(gcf, 'position', [0 0 1400 600]);
        %         figMeanTrace=figure;%plot mean trace
        %         set(gcf, 'position', [0 0 1400 300]);
        
        for ndff=1:1 %plot for each dff, see which is better and test whether 410 signals have similar trend
            if strcmp(behEventAlign,'delay onset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
                [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff{1,ndff}(nfiber,:), behEventFrameIndex,  frameNum );
            else
                [ dff_aligned, behEvent_aligned, licking_aligned ] = fAlignSigBehEvent( dff{1,ndff}(nfiber,:), behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
            end
            fig=figure;
            set(fig, 'Position', position{nevent,nfiber});
%             set(fig, 'PaperPosition', position{nevent,nfiber});
            ColLimit = prctile(dff{1,ndff}(nfiber,:)',98);
            titlestr={'Correct','Error','Miss','Violation'};
            %         titlestr=strcat(ImagingSetup,'-',fiberstr{nfiber},'-',titlestr);
            color_mean_trace=fMeanTraceColor(behrule,size(trialType,2));
            for nStim=1:size(trialType,2) %for each stimulus%[1,6]%just 2 end trials
                for  nResult=1:size(trialType,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
                    selectedTrialInd=trialType(nResult,nStim,:);
                    selectedTrialInd=logical(squeeze(selectedTrialInd))';
                    sot=behEvent_aligned.stimOnset(selectedTrialInd);% stim onset, white
                    flt=behEvent_aligned.lickFirst(selectedTrialInd);% first lick, black
                    flt_l=behEvent_aligned.lickFirst_left(selectedTrialInd);%left lick, magenta dots
                    flt_r=behEvent_aligned.lickFirst_right(selectedTrialInd);%right lick, red dots
                    rwt=behEvent_aligned.rewTime(selectedTrialInd);%reward, cyan dots
                    go=behEvent_aligned.go(selectedTrialInd);%go cue,white
                    neuralActivity=dff_aligned(selectedTrialInd,:);
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
                    go=go(I);
                    sot=sot(I);
                    flt=flt(I);
                    flt_l=flt_l(I);
                    flt_r=flt_r(I);
                    rwt=rwt(I);
                    figure(fig);
                    
                    subplot(size(trialType,2)+1,2*size(trialType,1),2*nResult-1+2*size(trialType,1)*(nStim-1));
                    hold on;
                    imagesc(neuralActivity(I,:));
                    %x_lim=get(gca,'Xlim');
                    x_lim=[0,sum(frameNum)];
                    %plot behavior event
                    plot(sot,[1:sum(selectedTrialInd)],'w.','markersize',8);
                    hold on;
                    plot(go,[1:sum(selectedTrialInd)],'w.','markersize',8);
                    plot(rwt,[1:sum(selectedTrialInd)],'c.','markersize',8);
                    hold on;
                    set(gca,'clim',[0 ColLimit]);
                    set(gca,'ytick',sum(selectedTrialInd),'yticklabel',sum(selectedTrialInd),'ydir','normal');
                    set(gca,'xtick',[]);
                    ColBar = colorbar;
                    set(ColBar,'position',[0.91 0.2400 0.02 0.1445]);
                    if sum(selectedTrialInd)>0
                        ylim([0.5 sum(selectedTrialInd)+0.5]);
                    end
                    if nStim == 1
                        title(titlestr{nResult});
                    end
                    if nResult==1
                        ylabel([num2str(Data_extract.Stimuli(nStim))]);% ' clicks/s'
                    end
                    xlim(x_lim);
                    subplot(size(trialType,2)+1,2*size(trialType,1),2*nResult+2*size(trialType,1)*(nStim-1));
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
                    set(gca,'ytick',[]);%licking raster have corresponding color plot label, so not need
                    set(gca,'xtick',[]);
                    xlim(x_lim);
                    if sum(selectedTrialInd)>0
                        ylim([0.5 sum(selectedTrialInd)+0.5]);
                    end
                    %plot mean trace
                    subplot(size(trialType,2)+1,2*size(trialType,1),2*nResult-1+2*size(trialType,1)*size(trialType,2));
                    %subplot(1,size(trialType,1),nResult);
                    [neuralActivityMean, neuralActivityCI]=fMean_CI(neuralActivity,0.05);
                    %                 xpatch=[1:size(neuralActivity,2), fliplr(1:size(neuralActivity,2))];
                    %                 ypatch=[neuralActivityCI(1,:),fliplr(neuralActivityCI(2,:))];
                    %                 p=patch(xpatch,ypatch,color_mean_trace{nStim});%plot confidence interval
                    %                 p.FaceAlpha=0.1;
                    %                 p.EdgeColor=color_mean_trace{nStim};%'none';
                    %                 hold on;
                    curve_meanTrace(nStim)=plot(1:size(neuralActivity,2),neuralActivityMean,'Color',color_mean_trace{nStim},'linewidth',2);
                    hold on;
%                     set(gca,'Ylim',[min(neuralActivityMean),max(neuralActivityMean)]);
                    y_lim=get(gca,'Ylim');
                    xlim(x_lim);
                    plot([frameNumTime(1)*round(1000/frT),frameNumTime(1)*round(1000/frT)],[y_lim(1),y_lim(end)],'k-');%align to a behavior event
                    xlabel(['time(s) from ',behEventAlign]);
                    if nResult==1
                        ylabel('\it\DeltaF/F');
                    end
                    set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
%                     set(gca,'ylim',[0,0.6]);
                    %                 set(gca,'FontName','Arial','FontSize',14);
                    box off;
                end
            end
            subplot(size(trialType,2)+1,2*size(trialType,1),2*size(trialType,1)*(size(trialType,2)+1)-1);
            leg=legend(curve_meanTrace(:),num2str(Data_extract.Stimuli(:)));
            set(leg,'position',[0.91 0.1 0.02 0.08]);
            legend boxoff;
            filestr=strcat('Summary of ',session_name,'-align-to-',behEventAlign,'-fiber-',fiberstr{nfiber});
            print(fig,filestr{1},'-dpdf');

            %         if strcmp(masklick,'yes')
            %             suptitle([dffstr{ndff},'-masklick']);
            %         else
            %             suptitle(dffstr{ndff});
            %         end
            %         suptitle(dffstr{ndff});
            
        end
    end
end
figBeh=figure;%plot curve
set(gcf, 'position',positionBeh);%控制fig尺寸
[animal_name,dataChoice]=fFindChoice(rootpath,1);
 p_ctrl=find(dataChoice(1,:,1)==0);%这里有bug，当坐满1000trial时是空值
 if isempty(p_ctrl) 
     p_ctrl=1000;
 end
 choice_ctrl=squeeze(dataChoice(1,1:p_ctrl-1,:));  
[correctProbe_ctrl,toneOct_ctrl,choicep_ctrl,f_ctrl,x_ctrl,yfit_ctrl]=fProbePerformance(choice_ctrl);
ctrl_curve=fPlotProbePerformance(figBeh,'k',correctProbe_ctrl,toneOct_ctrl,choicep_ctrl,x_ctrl,yfit_ctrl,animal_name);
% save([rootpath,'\SummaryOneSession.mat']);
end
function [curve] = fPlotProbePerformance(fig,color,correctProbe,toneOct,choicep,x,yfit,animal_name)
%x1 = linspace(min(toneOct)-0.1, max(toneOct)+0.1, 100);
    figure(fig);
    box off;
    points_legend=scatter(2.^choicep(:,1)*20,choicep(:,2)*100,correctProbe(1:size(correctProbe,1)-3,6),color,'filled');
    hold on;
    curve=plot(x ,yfit*100 ,color,'LineWidth',2);%标记该曲线，为legend提供名称
    hold on;
%     ylabel('% Choice of High Click Rate Side','FontName','Arial','FontSize',14);
    ylabel('% Choice High');
    xlabel('clicks/s');
%     title(['Psychometric curve'],'FontName','Arial','FontWeight','Bold','FontSize',16);
%     set(gca,'FontName','Arial','FontSize',14);
    set(gca, 'YLim',[0 100]);
    set(gca, 'XLim',[x(1) x(end)]); 
    set(gca,'XScale','log');
    set(gca,'XTick',[20,125],'LineWidth',1);
    set(gca,'XTickLabel',[20,125]);
    set(gca,'YTick',[0,50,100]);
    set(gca,'YTickLabel',[0,50,100]);
    plot([x(1),x(end)],[50,50],'--k','LineWidth',2);
    hold on;   
%     plot([50,50],[0,100],'--k','LineWidth',1);
    hold on;
end
function [correctProbe,toneOct,choicep,f,x,yfit] = fProbePerformance(choice)
    %click=[20 28 40 62 90 125];
    click=unique(choice(:,1));
    if length(click)==2
        warning('Not enough stimuli');
    end
    freq_n=length(click);
    correctProbe=zeros(freq_n+3,6);%每行表示一种频率，每列分别表示频率，正确率，miss率，正确率的std，violation率，correct+error的trial数量之和；最后另外加三行，分别表示end trial,和probe trial, all trials

    pclick=zeros(size(choice,1),freq_n);%用于记录每个click rate的trial位置
    for i=1:freq_n
        pclick(:,i)=double(choice(:,1)==click(i));
    end
    pclick(:,freq_n+1)=pclick(:,1)+pclick(:,freq_n);
    pclick(:,freq_n+2)=sum(pclick(:,2:freq_n-1),2);%对行求和
    pclick(:,freq_n+3)=pclick(:,freq_n+1)+pclick(:,freq_n+2);
    pperformance=zeros(size(choice,1),4);%用于记录各种类型的trial位置，包括correct,error,miss三列
    performance=[1 2 3 4];%分别表示CORRECT,ERROR,MISS,VIOLATION
    for i=1:4
        pperformance(:,i)=double(choice(:,2)==performance(i));
    end
    correctProbe(1:freq_n,1)=click;%最后两行由于是end 和probez，暂时不赋值了
    for i=1:freq_n+3
        do=(pperformance(:,1)+pperformance(:,2)).*pclick(:,i);
        notMiss=(do+pperformance(:,4)).*pclick(:,i);
        correctProbe(i,2)=sum(pperformance(:,1).*pclick(:,i))/sum(do);
        correctProbe(i,3)=sum(pperformance(:,3).*pclick(:,i))/sum(pclick(:,i));
        %计算各个频率的正确率的std
        correctNoMiss=pperformance(:,1).*pclick(:,i);
        correctNoMiss(~notMiss)=nan;
        correctProbe(i,4)=nanstd(correctNoMiss);
        correctProbe(i,5)=sum(pperformance(:,4).*pclick(:,i))/sum(notMiss);
        correctProbe(i,6)=sum(do);
    end
    
    %拟合
    toneOct  = log2(correctProbe(1:freq_n,1)/20);
    choicep =correctProbe(1:freq_n,1:2);
    choicep(:,1)=toneOct ;
    choicep(1:freq_n/2,2)=1-correctProbe(1:freq_n/2,2);
    ffun=fittype('a+b./(1+exp(-k.*(x-c)))','independent','x');
    f =fit(choicep(:,1),choicep(:,2),ffun,'Startpoint',[0,1,1,1]);
    slope=f.b*f.k/4;
    x =toneOct(1)-0.5:0.01:toneOct(end)+0.5;
    yfit =f(x);
    x=2.^x*20;
end
