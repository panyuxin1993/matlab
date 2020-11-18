%This code choose some sessions and Concatenate data together to a larger
%session and plot PSTH and Ending point trajectories, the purpose here is
%to increase trial number.
% the limitation here is only data from same animal can be grouped
dbstop if error
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
%set the criteria to choose sessions
animal_name='pyx';
experiment={'SC vglut2','SC vgat'};
neuralActivityCmp=cell(1,length(experiment));
for i_exp=1:length(experiment)
    ind_session=(T.delay_max==1.5).*(T.delay_min==0.3).*strcmp(T.experiment,experiment{i_exp}).*strcmp(T.used_as_data,'yes').*(T.retract==0);
    % ind_session=strcmp(T.animal,animal_name).*strcmp(T.used_as_data,'yes');.*(strcmp(T.brain_region,'bilateral
    % SC'))
    ind_session=find(ind_session);
    n_session=length(ind_session);
    animal_unique=unique(T.animal(ind_session));
    n_animal=length(animal_unique);
    n_site=zeros(2,1);%store n sites for soma and terminal
    dff_aligned_cat=cell(2,1);%two fibers so two cells to store,for multiple sessions, store 1-soma,2-terminal data
    neuralActivityCmp{1,i_exp}=cell(2,size(trialType_cat,1));%niber-by-size(trialType_cat,1)
    for nfiber=1:2
        trialType_cat=[];
        behEvent_aligned_cat=[];
        licking_aligned_cat=[];
        for i=1:n_animal
            ind_row=find(strcmp(T.animal,animal_unique{i}));
            if nfiber==1%soma
                if strcmp( T.brain_region(ind_row),'bilateral SC')
                    n_site(nfiber)=n_site(nfiber)+2;
                else
                    n_site(nfiber)=n_site(nfiber)+1;
                end
            elseif nfiber==2%terminal
                if strcmp( T.brain_region(ind_row),'bilateral SC')
                    n_site(nfiber)=n_site(nfiber);
                else
                    n_site(nfiber)=n_site(nfiber)+1;
                end
            end
        end
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
            rootpath=rootpath{1};
            SavingFolder=rootpath;
            DatePath=rootpath;
            files = dir(strcat(rootpath,'\*Virables.mat'));
            if length(files)==1
                behaviorFile=files(1).name;
            else
                warning('Multiple beh files');
            end
            cd(DatePath);
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
            
            Session470LEDstartFrame=1;% usually 205
            Session410LEDstartFrame=2;
            
            fiberstr={'Soma','Terminal'};
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
            for i_fiber=1:length(fiberSide)
                [trialType,behrule,trialTypeStr] = fGetTrialType( Data_extract,[],3,'matrix',fiberSide{i_fiber});%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
                frT=FrameTime*2;%2 channel, so framerate and frametime should be half
                frameNumTime=[1,2.5];%from 2s before align point to 5s after align point
                frameNum=double(round(frameNumTime*1000/frT));
                [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime(Data_extract,round(double(TrialStart_FrameCount')/2),frT);
                trialType_cat=cat(3,trialType_cat,trialType);
                %% concatenate dff, trial type maxtrix and behavior event, etc.
                dffstr={'dff baseline correction first','dff motion correction first','dff470','dff410'};
                behEventAlign='stimOnset';%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
                behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
                masklick='yes';%'yes' if mask lickings or 'no' when want to see all activity
                for ndff=[1] %plot for each dff, see which is better and test whether 410 signals have similar trend
                    if strcmp(behEventAlign,'stimOnset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
                        [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff{1,ndff}(ind_fiber(i_fiber),:), behEventFrameIndex,  frameNum );
                    else
                        [ dff_aligned, behEvent_aligned, licking_aligned ] = fAlignSigBehEvent( dff{1,ndff}(ind_fiber(i_fiber),:), behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
                    end
                    dff_aligned_cat{nfiber}=cat(1,dff_aligned_cat{nfiber},dff_aligned);
                    if strcmp(fiberSide{i_fiber},'right')
                        [ behEvent_aligned.lickFirst_left, behEvent_aligned.lickFirst_right]=deal(behEvent_aligned.lickFirst_right,behEvent_aligned.lickFirst_left);%swarp value
                        [licking_aligned.leftLick,licking_aligned.rightLick]=deal(licking_aligned.rightLick,licking_aligned.leftLick);
                    end
                    behEvent_aligned_cat=fMergeStruct(behEvent_aligned_cat,behEvent_aligned);
                    licking_aligned_cat=fMergeStruct(licking_aligned_cat,licking_aligned);
                end
            end
        end
        %plot raster and mean trace for large session
        for ndff=[1]
            if isempty(dff_aligned_cat{nfiber})%if no data, then end this loop
                break;
            end
            %         fig=figure;%plot raster
            figMeanTrace=figure;%plot mean trace
            set(gcf, 'position', [0 0 1400 300]);
            if size(trialType_cat,2)==2%miss && grouped with choice rather than stimuli
                figMeanTraceDiff=figure;%plot difference between mean trace
                set(gcf, 'position', [0 0 1400 300]);
            end
            titlestr={'Correct','Error','Miss','Violation'};
            ColLimit = prctile(dff_aligned_cat{nfiber}',98,'all');%here data is a matrix rather than a vector, so use 'all' to find the prctile of all data
            %         titlestr=strcat(animal_name,'-',fiberstr{nfiber},'-',titlestr);
            for nResult=1:size(trialType_cat,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
                if size(trialType_cat,2)==2%miss && grouped with choice rather than stimuli
                    neuralActivityMeanDiff=cell(1,2);
                end
                for nStim=1:size(trialType_cat,2) %for each stimulus%[1,6]%just 2 end trials
                    if nResult==3 && size(trialType_cat,2)==2%miss && grouped with choice rather than stimuli
                        color_mean_trace={[0,0,0],[0,0,0]};
                    else
                        color_mean_trace=fMeanTraceColor(behrule,size(trialType_cat,2));
                    end
                    selectedTrialInd=trialType_cat(nResult,nStim,:);
                    selectedTrialInd=logical(squeeze(selectedTrialInd))';
                    sot=behEvent_aligned_cat.stimOnset(selectedTrialInd);% stim onset, white
                    flt=behEvent_aligned_cat.lickFirst(selectedTrialInd);% first lick, black
                    flt_l=behEvent_aligned_cat.lickFirst_left(selectedTrialInd);%left lick, magenta dots(ipsi lick)
                    flt_r=behEvent_aligned_cat.lickFirst_right(selectedTrialInd);%right lick, red dots(contra lick)
                    rwt=behEvent_aligned_cat.rewTime(selectedTrialInd);%reward, cyan dots
                    go=behEvent_aligned_cat.go(selectedTrialInd);%go cue,white
                    neuralActivity=dff_aligned_cat{nfiber}(selectedTrialInd,:);
                    temp=find(selectedTrialInd);
                    leftLick=cell(1,length(temp));
                    rightLick=cell(1,length(temp));
                    for i=1:length(temp)
                        leftLick{i}= licking_aligned_cat.leftLick{temp(i)};
                        rightLick{i}=licking_aligned_cat.rightLick{temp(i)};
                    end
                    %                 %plot raster
                    %                 if strcmp(behEventSort,'first lick')
                    %                     [B,I]=sort(flt);
                    %                 elseif strcmp(behEventSort,'reward')
                    %                     [B,I]=sort(rwt);
                    %                 elseif strcmp(behEventSort,'go cue')
                    %                     [B,I]=sort(go);
                    %                 end
                    %                 go=go(I);
                    %                 sot=sot(I);
                    %                 flt=flt(I);
                    %                 flt_l=flt_l(I);
                    %                 flt_r=flt_r(I);
                    %                 rwt=rwt(I);
                    %                 figure(fig);
                    %                 subplot(size(trialType_cat,2)+1,2*size(trialType_cat,1),2*nResult-1+2*size(trialType_cat,1)*(nStim-1));
                    %                 imagesc(neuralActivity(I,:));
                    %                 hold on;
                    %                 x_lim=[0,sum(frameNum)];
                    %                 %plot behavior event
                    %                 plot(sot,1:sum(selectedTrialInd),'w.','markersize',8);
                    %                 hold on;
                    %                 plot(go,1:sum(selectedTrialInd),'w.','markersize',8);
                    %                 plot(rwt,1:sum(selectedTrialInd),'c.','markersize',8);
                    %                 hold on;
                    %                 set(gca,'clim',[0 ColLimit]);
                    %                 set(gca,'ytick',sum(selectedTrialInd),'yticklabel',sum(selectedTrialInd),'ydir','normal');
                    %                 set(gca,'xtick',[]);
                    %                 ColBar = colorbar;
                    %                 set(ColBar,'position',[0.91 0.2400 0.02 0.1445]);
                    %                 if sum(selectedTrialInd)>0
                    %                     ylim([0.5 sum(selectedTrialInd)+0.5]);
                    %                 end
                    %                 if nStim == 1
                    %                     title(titlestr{nResult});
                    %                 end
                    %                 if nResult==1&&strcmp(trialTypeStr,'stimuli')
                    %                     ylabel([num2str(Data_extract.Stimuli(nStim))]);% ' clicks/s'
                    %                 elseif nResult==1&&strcmp(trialTypeStr,'first lick')
                    %                     ylabelstr={'ipsi lick first','contra lick first'};
                    %                     ylabel(ylabelstr{nStim});
                    %                 end
                    %                 xlim(x_lim);
                    %                 subplot(size(trialType,2)+1,2*size(trialType,1),2*nResult+2*size(trialType,1)*(nStim-1));
                    %                 for i=1:length(leftLick)%each trial as a row
                    %                     for jl=1:length(leftLick{I(i)})
                    %                         line([leftLick{I(i)}(jl),leftLick{I(i)}(jl)],[i-0.5,i+0.5],'color','b','linewidth',1);%left lick
                    %                         hold on;
                    %                     end
                    %                     for jr=1:length(rightLick{I(i)})
                    %                         line([rightLick{I(i)}(jr),rightLick{I(i)}(jr)],[i-0.5,i+0.5],'color','r','linewidth',1);%right lick
                    %                         hold on;
                    %                     end
                    %                 end
                    %                 set(gca,'ytick',[]);%licking raster have corresponding color plot label, so not need
                    %                 set(gca,'xtick',[]);
                    %                 xlim(x_lim);
                    %                 if sum(selectedTrialInd)>0
                    %                     ylim([0.5 sum(selectedTrialInd)+0.5]);
                    %                 end
                    %plot mean trace
                    figure(figMeanTrace);%save mean trace
                    subplot(1,size(trialType_cat,1),nResult);
                    [neuralActivityMean, neuralActivityCI]=fMean_CI(neuralActivity,0.05);
                    if exist('neuralActivityMeanDiff','var')
                        neuralActivityMeanDiff{nStim}=neuralActivityMean;
                    end
                    %shadow as ci
                    
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
                    %
                    [ ~,~,neuralActivityGroupedtemp ] = fGetEndDff( neuralActivity ,'smooth','bin',0.3*1000/frT);
                    [ timepoint,endDff,neuralActivityGrouped ] = fGetEndDff( neuralActivityGroupedtemp ,'raw');
                    for i=1:size(neuralActivityGrouped,1)
                        plot(neuralActivityGrouped(i,:),'Color',color_mean_trace{nStim},'linewidth',1.5);
                    end
                    scatter(timepoint,endDff,50,color_mean_trace{nStim},'MarkerFaceColor','w');
                    %}
                    if nResult==1
                        ylabel('\it\DeltaF/F');
                    end
                    set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
                    x_lim=get(gca,'Xlim');
                    xlim([x_lim(1),size(neuralActivity,2)]);
                    %                     ylim([nanmin(nanmin(neuralActivityGrouped)),nanmax(nanmax(neuralActivityGrouped))]);
                    set(gca,'FontName','Arial','FontSize',14);
                    box off;
                end
                if exist('neuralActivityMeanDiff','var')
                    neuralActivityCmp{1,i_exp}{nfiber,nResult}=neuralActivityMeanDiff{2}-neuralActivityMeanDiff{1};
                end
                if size(trialType_cat,2)==2%miss && grouped with choice rather than stimuli
%                     figure(figMeanTraceDiff);
%                     subplot(1,size(trialType_cat,1),nResult);
%                     curve_meanTraceDiff=plot(1:size(neuralActivityMeanDiff{1},2),neuralActivityMeanDiff{2}-neuralActivityMeanDiff{1},'Color',[0,0,0],'linewidth',2);
%                     hold on;
%                     y_lim=get(gca,'Ylim');
%                     plot([frameNumTime(1)*round(1000/frT),frameNumTime(1)*round(1000/frT)],[y_lim(1),y_lim(2)],'k-');%align to a behavior event
%                     xlabel(['time(s) from ',behEventAlign]);
%                     title(titlestr{nResult});
                end
            end
            figure(figMeanTrace);
            subplot(1,size(trialType_cat,1),3);%miss plot has most space
            if contains(trialTypeStr,'stimuli')
                h=legend(curve_meanTrace(:),num2str(Data_extract.Stimuli(:)),'Location','best');
            elseif contains(trialTypeStr,'first lick')
                h=legend(curve_meanTrace(:),{'ipsi lick first','contra lick first'},'Location','best');
            end
            set(h,'box','off');
            x_lim=get(gca,'Xlim');
            y_lim=get(gca,'Ylim');
            text(x_lim(1),y_lim(end),strcat('n=',num2str(n_session),' sessions,',num2str(n_site(nfiber)),' site,',num2str(n_animal),' animal'),'FontName','Arial','FontSize',14);
            if n_animal==1%only one animal
                if strcmp(masklick,'yes')
                    suptitle(strcat(animal_name,'-',fiberstr{nfiber},'-',dffstr{ndff},'-masklick'));
                else
                    suptitle(strcat(animal_name,'-',fiberstr{nfiber},'-',dffstr{ndff}));
                end
            else%use common feature, e.g. experiment/group name
                if strcmp(masklick,'yes')
                    suptitle(strcat(experiment,'-',fiberstr{nfiber},'-',dffstr{ndff},'-masklick'));
                else
                    suptitle(strcat(experiment,'-',fiberstr{nfiber},'-',dffstr{ndff}));
                end
            end
            saveas(figMeanTrace,[ImagingSetup,'-',fiberstr{nfiber},dffstr{ndff},'-algin to ',behEventAlign,'-sort ',behEventSort,'-masklick',masklick,'mean trace with end point dff.fig'],'fig');
            %             close all;
        end
    end
end
figDffCmp=figure;
color={[1,0,0],[0,0,1]};
for nfiber=1:2
    for i_exp=1:length(experiment)
        for nResult=1:4
            figure(figDffCmp);
            subplot(2,4,nResult+nfiber*4-4);
            curve_cmpMeanTraceDiff(i_exp)=plot(1:size(neuralActivityCmp{1,i_exp}{nfiber,nResult},2),neuralActivityCmp{1,i_exp}{nfiber,nResult},'Color',color{i_exp},'linewidth',2);
            hold on;
            y_lim=get(gca,'Ylim');
            plot([frameNumTime(1)*round(1000/frT),frameNumTime(1)*round(1000/frT)],[y_lim(1),y_lim(2)],'k-');%align to a behavior event
            xlabel(['time(s) from ',behEventAlign]);
            title(titlestr{nResult});
            set(gca,'FontName','Arial','FontSize',14);
            set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
            x_lim=get(gca,'Xlim');
            xlim([x_lim(1),size(neuralActivity,2)]);
        end
    end
end
subplot(2,4,4);
 h=legend(curve_cmpMeanTraceDiff(:),{'vglut2','vgat'},'Location','best');
set(h,'box','off');
subplot(2,4,1);
ylabel('{\it\delta} activities(contra-ipsi,{\it\DeltaF/F}) of soma');
subplot(2,4,5);
ylabel('{\it\delta} activities(contra-ipsi,{\it\DeltaF/F}) of terminal');

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