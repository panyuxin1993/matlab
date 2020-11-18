
fn = 'AlineCaSig_DfModef0SubPreTone_h37_20161113';

saveFig = 0;  % 0: don't save fig            1: save fig
savePPT = 0;  % 0: don't save picture to ppt            1: save picture to ppt

load(fn);
%%  this block use CaSig data
temp = dff_mode; %dff_subtr_preSoundmean; % choose calcium signal or spikes data ##### dff_subtr_preSoundmean   Spikes
                                                                     %##### nSpikes_Smooth
                                                                     
length_frames_T = length( temp{1}(1,[F_num_T(1,1):F_num_T(1,2)]));
length_frames_A = length( temp{1}(1,[F_num_A(1,1):F_num_A(1,2)]));
for nROI = 1:nROIs
    max_temp = 0;
    min_temp = 0;
    for triN = 1:TrialN
        temp_dff_1{nROI}(triN,:) = temp{nROI}(triN,[F_num_T(triN,1):F_num_T(triN,2)]);  % get the same frames before and after tone onset  
        if ~Miss_Ind(triN)
            temp_dff_2{nROI}(triN,:) = temp{nROI}(triN,[F_num_A(triN,1):F_num_A(triN,2)]);  % get the same frames before and after Action  
        else
            temp_dff_2{nROI}(triN,:) = NaN;
        end
    end    
    ii = 0;
    tt = 1;
    while ii<=0
        col_lim(nROI) = prctile(reshape(temp_dff_1{nROI},1,[]),90+(tt-1));  
        ii = col_lim(nROI);
        tt = tt+1;
    end
    max_temp = max([max(temp_dff_1{nROI}) max_temp]);
    min_temp = min([min(temp_dff_1{nROI}) min_temp]);    
    maxCaS(nROI) = max_temp;
    minCaS(nROI) = min_temp; 
end

%%  this block use Spike data

for N = 1:size(Spikes,2)
    temp{N}(:,:) = Spikes(:,N,:);
end
                                                                     
length_frames_T = length( temp{1}(1,[F_num_T(1,1):F_num_T(1,2)]));
length_frames_A = length( temp{1}(1,[F_num_A(1,1):F_num_A(1,2)]));
for nROI = 1:size(Spikes,2)
    max_temp = 0;
    min_temp = 0;
    for triN = 1:size(Spikes,1)
        temp_dff_1{nROI}(triN,:) = temp{nROI}(triN,[F_num_T(triN,1):F_num_T(triN,2)]);  % get the same frames before and after tone onset  
        if ~Miss_Ind(triN)
            temp_dff_2{nROI}(triN,:) = temp{nROI}(triN,[F_num_A(triN,1):F_num_A(triN,2)]);  % get the same frames before and after Action  
        else
            temp_dff_2{nROI}(triN,:) = NaN;
        end
    end    
    ii = 0;
    tt = 1;
    while ii<=0
        col_lim(nROI) = prctile(reshape(temp_dff_1{nROI},1,[]),90+(tt-1));  
        ii = col_lim(nROI);
        tt = tt+1;
    end
    max_temp = max([max(temp_dff_1{nROI}) max_temp]);
    min_temp = min([min(temp_dff_1{nROI}) min_temp]);    
    maxCaS(nROI) = max_temp;
    minCaS(nROI) = min_temp; 
end
%%
Trials_sum(1,:) = sum(TrialType(1:8,:)');
Trials_sum(2,:) = sum(TrialType(9:16,:)');
Trials_sum(3,:) = sum(TrialType(17:24,:)');
blank_val = round(sum(Trials_sum)/15);

Y_step(1,:) =  Trials_sum(1,:) + blank_val;
Y_step(2,:) =  Y_step(1,:) + Trials_sum(2,:);
Y_step(3,:) =  Y_step(2,:) + blank_val;
Y_step(4,:) =  Y_step(3,:) + Trials_sum(3,:);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting color plot sorted by answer
savepath = 'F:\Process DATA\2P\h37\20161113\ColorPlots\SortByAciton';
saveFig = 0;  % 0: don't save fig            1: save fig
savePPT = 0;  % 0: don't save picture to ppt            1: save picture to ppt

% isOpen = exportToPPTX();
% if ~isempty(isOpen)
%    exportToPPTX('close'); 
% end

% exportToPPTX('open',[fn '.pptx']);

for nROI = 60%:nROIs    
    fig1 = figure;
    set(gcf,'color','w');
    YlabelName = {'Corr','Wro','Mis'}; 
    set(gcf,'position',[10 100 1600 400],'color','w');
    for p = 1:length(Fre) 
        pl(p) = subplot(1,length(Fre),p);
        title([num2str(Fre(p)) ' Hz']);
        hold on;
        set(gca,'fontsize',10,'fontweight','bold');
        if p ==1
            text(-120,1,(YlabelName{1}),'fontsize',15,'fontweight','bold','rotation',90);
            text(-120,Y_step(1,p),(YlabelName{2}),'fontsize',15,'fontweight','bold','rotation',90);
            text(-120,Y_step(3,p),(YlabelName{3}),'fontsize',15,'fontweight','bold','rotation',90);
        end
        if p == 1
              xlabel('(T/s)','fontsize',12,'fontweight','bold');
        end
        patch([0 length_frames_T length_frames_T 0],[Trials_sum(1,p)+0.5 Trials_sum(1,p)+0.5 Y_step(1,p)+0.5 Y_step(1,p)+0.5],'w','edgecolor','k');
        patch([0 length_frames_T length_frames_T 0],[Y_step(2,p)+0.5 Y_step(2,p)+0.5 Y_step(3,p)+0.5 Y_step(3,p)+0.5],'w','edgecolor','k');     
        
        
        for nn = 1:3
            CaS_temp_1 = temp_dff_1{nROI}(TrialType(p+(nn-1)*length(Fre),:),:);
            [SortT_temp{nn} OnsetInd] = sort(AnswFramAlineByTone(TrialType(p+(nn-1)*length(Fre),:)));
            CaS_temp_2{nn} = CaS_temp_1(OnsetInd,:); % Find CaSig sort by answer time
            if nn == 1
                imagesc(CaS_temp_2{1}); 
                h2 = plot(SortT_temp{nn},[1:1:Trials_sum(nn,p)],'w.','linestyle','none','markersize',6);
            elseif nn == 2
                imagesc([1:length_frames_T],[Y_step(1,p)+1:Y_step(2,p)],CaS_temp_2{2});
                h2 = plot(SortT_temp{nn},[Y_step(1,p)+1:1:Y_step(2,p)],'w.','linestyle','none','markersize',8);
            elseif nn == 3
                imagesc([1:length_frames_T],[Y_step(3,p)+1:Y_step(4,p)],CaS_temp_2{3}); 
            end
        end
        h1 = line([MinOnsFrame MinOnsFrame],[0 Y_step(4,p)+1],'color','w','linewidth',1);
    
        ylim([0.5 Y_step(4,p)+0.5]);
        xlim([0 length_frames_T])       

        colormap jet;
        set(gca,'xtick',[0:1000/FrameTime:7000/FrameTime],'XTickLabel',[0:1:7],'ytick',Y_step(4,p),'YTickLabel',sum(Trials_sum(:,p)),'clim',[0 col_lim(nROI)]);
    end    
    
    if isempty(Data_extract.Animal_Name_Setting)
        suptitle([Data_extract.Animal_Name_File(4:6)  '-' Data_extract.Experiment_date '-Sorted by Action- ROI-' num2str(nROI)]);
    else
        suptitle([Data_extract.Animal_Name_Setting  '-' Data_extract.Experiment_date '-Sorted by Action- ROI-' num2str(nROI)]);
    end
    if p == length(Fre)
        PO = get(pl(p),'pos');
        set(pl(p),'pos',PO);
        hcolbar = colorbar;
        get(hcolbar,'position');
        set(hcolbar,'position',[0.94 0.12 0.015 0.12]);
    end
%     if saveFig == 1
%         saveas(gcf,[savepath '\' fn([8:end]) '-ColorPlot_SorSub-' num2str(nROI,'%03d') '.fig']);
%     end
%     if savePPT == 1
%         exportToPPTX('switchslide',nROI);
%         exportToPPTX('addpicture',fig1,'position',[2 0 8 2]);       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     end
%     close;

end

% exportToPPTX('saveandclose');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting trce plot alined by tone

fn = 'dff_ZL_h37_20161113_field1_d167_3x_Same_ROI_for_2nd_of_3_consecutive_days';
saveFig = 0;  % 0: don't save fig            1: save fig
savePPT = 0;  % 0: don't save picture to ppt            1: save picture to ppt
SmoBand = 21; % set a smooth bin size   1: means don't smooth
plotSEM = 0;  % whether plot sem
LinWidth = 3; % trace linewidth

% if exist('MeanTracePlots','dir') == 0
%     mkdir('MeanTracePlots');   
% end
% cd ('MeanTracePlots');
% 
% if exist('MeanTraceAlinedByTone') == 0
%     mkdir('MeanTraceAlinedByTone');
% end
% cd ('MeanTraceAlinedByTone');
% PathName = pwd;
% cd ..
% cd ..
% isOpen = exportToPPTX();
% if ~isempty(isOpen)
%    exportToPPTX('close'); 
% end
close;
% exportToPPTX('open',[fn '.pptx']);
TracelabelName = {'Right','Left','Mis'}; 
for nROI = 60%:nROIs    
    fig1 = figure;
    set(gcf,'color','w');    
    set(gcf,'position',[10 100 1600 200],'color','w');
    for n = 1:length(Fre) 
        for m = 1:3
            CaS_temp{n}{m} = temp_dff_1{nROI}(TrialType(n+(m-1)*length(Fre),:),:);
            if Trials_sum(m,n)>3
                max_1(n,m) = max(smooth(mean(CaS_temp{n}{m}),SmoBand));
                min_1(n,m) = min(smooth(mean(CaS_temp{n}{m}),SmoBand));
            else
                max_1(n,m) = NaN;
                min_1(n,m) = NaN;  
            end
        end  
    end
    for p = 1:length(Fre) 
        if p <= length(Fre)/2 
            ColorName = {'k','r','b'};
        else
            ColorName = {'r','k','b'};
        end
        pl(p) = subplot(1,length(Fre),p);
        title([num2str(Fre(p)) ' Hz']);
        hold on;
        set(gca,'fontsize',10,'fontweight','bold');
        patch([MinOnsFrame MinOnsFrame+round(300/FrameTime) MinOnsFrame+round(300/FrameTime) MinOnsFrame],...
            [min(min(min_1)) min(min(min_1)) max(max(max_1)) max(max(max_1))],[.8 .8 .8],'edgecolor',[.8 .8 .8]);
        
        if p == 1
              xlabel('(T/s)','fontsize',12,'fontweight','bold');
        end
        for nn = 1:3
           Sem_temp_T(nn,:) = std(CaS_temp{p}{nn})/(Trials_sum(nn,p))^.5;
           Mean_T(nn,:) = smooth(mean(CaS_temp{p}{nn}),SmoBand);
           Sem_T(nn,:) = [Mean_T(nn,:)-Sem_temp_T(nn,:) fliplr(Mean_T(nn,:)+Sem_temp_T(nn,:))];
        end
        for nn = 1:3
            if Trials_sum(nn,p) > 3 % only plot while trials is more than 3
                if plotSEM
                    PatchPlot(nn) = patch([1:length(Mean_T(nn,:)) fliplr(1:length(Mean_T(nn,:)))],Sem_T(nn,:),ColorName{nn},'edgecolor','none');
                end
                PL_T(nn) = plot(Mean_T(nn,:),'color',ColorName{nn},'linewidth',LinWidth);
            end            
        end  
        alpha(.3);
        set(gca,'xtick',[0:1000/FrameTime:7000/FrameTime],'XTickLabel',[0:1:7]);
        ylim([min(min(min_1)) max(max(max_1))]);
        xlim([0 length_frames_T]) 
    end    
    if p == length(Fre)
        legend(PL_T,TracelabelName,'Position',[0.89 0.8 0.05 0.05])
        legend('boxoff')
    end
       
    if isempty(Data_extract.Animal_Name_Setting)
        suptitle([Data_extract.Animal_Name_File(4:6)  '-' Data_extract.Experiment_date '-Alined by Tone- ROI-' num2str(nROI)]);
    else
        suptitle([Data_extract.Animal_Name_Setting  '-' Data_extract.Experiment_date '-Alined by Tone- ROI-' num2str(nROI)]);
    end

%     if saveFig == 1
%         saveas(gcf,[PathName '\' fn([8:end]) '-TraceByT' num2str(nROI,'%03d') '.fig']);
%     end
%     if savePPT == 1
%         exportToPPTX('switchslide',nROI);
%         exportToPPTX('addpicture',gcf,'position',[2 2 8 1]);       
%     end
%     close;

end

exportToPPTX('saveandclose');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting trce plot alined by action
fn = 'dff_ZL_h37_20161113_field1_d167_3x_Same_ROI_for_2nd_of_3_consecutive_days';
saveFig = 1;  % 0: don't save fig            1: save fig
savePPT = 1;  % 0: don't save picture to ppt            1: save picture to ppt
SmoBand = 11; % set a smooth bin size   1: means don't smooth
plotSEM = 0;  % whether plot sem
LinWidth = 4; % trace linewidth

if exist('MeanTracePlots','dir') == 0
    mkdir('MeanTracePlots');   
end
cd ('MeanTracePlots');

if exist('MeanTraceAlinedByAction') == 0
    mkdir('MeanTraceAlinedByAction');
end
cd ('MeanTraceAlinedByAction');
PathName = pwd;
cd ..
cd ..

isOpen = exportToPPTX();
if ~isempty(isOpen)
   exportToPPTX('close'); 
end
 
exportToPPTX('open',[fn '.pptx']);
TracelabelName = {'Right','Left'}; 
for nROI = 1:nROIs    
    fig1 = figure;
    set(gcf,'color','w');

    set(gcf,'position',[10 100 1600 200],'color','w');
    for n = 1:length(Fre) 
        for m = 1:2
            CaS_temp{n}{m} = temp_dff_2{nROI}(TrialType(n+(m-1)*length(Fre),:),:);
            if Trials_sum(m,n)>3
                max_1(n,m) = max(mean(CaS_temp{n}{m}));
                min_1(n,m) = min(mean(CaS_temp{n}{m}));
            else
                max_1(n,m) = NaN;
                min_1(n,m) = NaN;  
            end            
        end  
    end
    for p = 1:length(Fre) 
        if p <= length(Fre)/2 
            ColorName = {'k','r','b'};
        else
            ColorName = {'r','k','b'};
        end
        pl(p) = subplot(1,length(Fre),p);
        title([num2str(Fre(p)) ' Hz']);
        hold on;
        set(gca,'fontsize',10,'fontweight','bold');
        line([MinActFrame MinActFrame],...
            [min(min(min_1)) max(max(max_1))],'color',[.8 .8 .8],'linewidth',2);
        alpha(.5);
        if p == 1
              xlabel('(T/s)','fontsize',12,'fontweight','bold');
        end
        for nn = 1:2
           Sem_temp_A(nn,:) = std(CaS_temp{p}{nn})/(Trials_sum(nn,p))^.5;
           Mean_A(nn,:) = smooth(mean(CaS_temp{p}{nn}),SmoBand);
           Sem_A(nn,:) = [Mean_A(nn,:)-Sem_temp_A(nn,:) fliplr(Mean_A(nn,:)+Sem_temp_A(nn,:))];
        end
        for nn = 1:2
            if Trials_sum(nn,p) > 3 % only plot while trials is more than 3
                if plotSEM
                    PatchPlot(nn) = patch([1:length(Mean_A(nn,:)) fliplr(1:length(Mean_A(nn,:)))],Sem_A(nn,:),ColorName{nn},'edgecolor','none');
                end
                PL_A(nn) = plot(Mean_A(nn,:),'color',ColorName{nn},,'linewidth',LinWidth);
            end            
        end
        alpha 0.3

        set(gca,'xtick',[0:1000/FrameTime:7000/FrameTime],'XTickLabel',[0:1:7]);
        ylim([min(min(min_1)) max(max(max_1))]);
        xlim([0 length_frames_A]) 
    end    
    if p == length(Fre)
        legend(PL_A,TracelabelName,'Position',[0.89 0.8 0.05 0.05])
        legend('boxoff')
    end
       
    if isempty(Data_extract.Animal_Name_Setting)
        suptitle([Data_extract.Animal_Name_File(4:6)  '-' Data_extract.Experiment_date '-Alined by Action- ROI-' num2str(nROI)]);
    else
        suptitle([Data_extract.Animal_Name_Setting  '-' Data_extract.Experiment_date '-Alined by Action- ROI-' num2str(nROI)]);
    end

    if saveFig == 1
        saveas(fig1,[PathName '\' fn([8:end]) '-TraceByA' num2str(nROI,'%03d') '.fig']);
    end
    if savePPT == 1
        exportToPPTX('switchslide',nROI);
        exportToPPTX('addpicture',fig1,'position',[2 3 8 1]);       
    end
    close;

end
exportToPPTX('saveandclose');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting color plot unsorted 
fn = 'dff_ZL_h37_20161113_field1_d167_3x_Same_ROI_for_2nd_of_3_consecutive_days';
saveFig = 1;  % 0: don't save fig            1: save fig
savePPT = 1;  % 0: don't save picture to ppt            1: save picture to ppt

if exist('ColorPlots','dir') == 0
    mkdir('ColorPlots');   
end
cd ('ColorPlots');

if exist('UnsortedColPlot') == 0
    mkdir('UnsortedColPlot');
end
cd ('UnsortedColPlot');
PathName = pwd;
cd ..
cd ..


isOpen = exportToPPTX();
if ~isempty(isOpen)
   exportToPPTX('close'); 
end
exportToPPTX('open',[fn '.pptx']);

for nROI = 1:nROIs    
    fig1 = figure;
    set(gcf,'color','w');
    YlabelName = {'Corr','Wro','Mis'}; 
    set(gcf,'position',[10 100 1600 400],'color','w');
    for p = 1:length(Fre) 
        pl(p) = subplot(1,length(Fre),p);
        title([num2str(Fre(p)) ' Hz']);
        hold on;
        set(gca,'fontsize',10,'fontweight','bold');
        if p ==1
            text(-120,1,(YlabelName{1}),'fontsize',15,'fontweight','bold','rotation',90);
            text(-120,Y_step(1,p),(YlabelName{2}),'fontsize',15,'fontweight','bold','rotation',90);
            text(-120,Y_step(3,p),(YlabelName{3}),'fontsize',15,'fontweight','bold','rotation',90);
        end
        if p == 1
              xlabel('(T/s)','fontsize',12,'fontweight','bold');
        end
        patch([0 length_frames_T length_frames_T 0],[Trials_sum(1,p)+0.5 Trials_sum(1,p)+0.5 Y_step(1,p)+0.5 Y_step(1,p)+0.5],'w','edgecolor','k');
        patch([0 length_frames_T length_frames_T 0],[Y_step(2,p)+0.5 Y_step(2,p)+0.5 Y_step(3,p)+0.5 Y_step(3,p)+0.5],'w','edgecolor','k');     
        
        for nn = 1:3
            CaS_temp_3 = temp_dff_1{nROI}(TrialType(p+(nn-1)*length(Fre),:),:);
            AnswF = AnswFramAlineByTone(TrialType(p+(nn-1)*length(Fre),:));
            if nn == 1
                imagesc(CaS_temp_3);      
                plot(AnswF,[1:Trials_sum(nn,p)],'w.','linestyle','none','markersize',6);
            elseif nn == 2
                imagesc([1:length_frames_T],[Y_step(1,p)+1:Y_step(2,p)],CaS_temp_3);
                plot(AnswF,[Y_step(1,p)+1:1:Y_step(2,p)],'w.','linestyle','none','markersize',6);
            elseif nn == 3
                imagesc([1:length_frames_T],[Y_step(3,p)+1:Y_step(4,p)],CaS_temp_3); 
            end
            
        end        

        h1 = line([MinOnsFrame MinOnsFrame],[0 Y_step(4,p)+1],'color','w','linewidth',1);
        
        ylim([0.5 Y_step(4,p)+0.5]);
        xlim([0 length_frames_T])       

        colormap jet;
        set(gca,'xtick',[0:1000/FrameTime:7000/FrameTime],'XTickLabel',[0:1:7],'ytick',Y_step(4,p),'YTickLabel',sum(Trials_sum(:,p)),'clim',[0 col_lim(nROI)]);
    end    
    
    if isempty(Data_extract.Animal_Name_Setting)
        suptitle([Data_extract.Animal_Name_File(4:6)  '-' Data_extract.Experiment_date '-UnSorted-ROI-' num2str(nROI)]);
    else
        suptitle([Data_extract.Animal_Name_Setting  '-' Data_extract.Experiment_date '-UnSorted-ROI-' num2str(nROI)]);
    end
    if p == length(Fre)
        PO = get(pl(p),'pos');
        set(pl(p),'pos',PO);
        hcolbar = colorbar;
        get(hcolbar,'position');
        set(hcolbar,'position',[0.94 0.12 0.015 0.12]);
    end
    if saveFig == 1
        saveas(fig1,[PathName '\' fn([8:end]) '-ColorPlot-UnSorSub' num2str(nROI,'%03d') '.fig']);
    end
    if savePPT == 1
        exportToPPTX('switchslide',nROI);
        exportToPPTX('addpicture',fig1,'position',[2 4 8 2]);       
    end
    close;
end

exportToPPTX('saveandclose');



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting trial sequency with color plot 
saveFig = 1;  % 0: don't save fig            1: save fig
savePPT = 1;  % 0: don't save picture to ppt            1: save picture to ppt

if exist('ColorPlots','dir') == 0
    mkdir('ColorPlots');   
end
cd ('ColorPlots');

PathName = pwd;
cd ..


isOpen = exportToPPTX();
if ~isempty(isOpen)
   exportToPPTX('close'); 
end
exportToPPTX('open',[fn '.pptx']);  

fig1 = figure;       
set(gcf,'color','w','position',[10 100 400 400],'color','w');
YlabelName = {'Corr','Wro','Mis'}; 
for p = 1:length(Fre) 
    pl(p) = subplot(1,length(Fre),p);
    hold on;
    set(gca,'fontsize',10,'fontweight','bold');
    if p ==1
        text(-120,1,(YlabelName{1}),'fontsize',15,'fontweight','bold','rotation',90);
        text(-120,Y_step(1,p),(YlabelName{2}),'fontsize',15,'fontweight','bold','rotation',90);
        text(-120,Y_step(3,p),(YlabelName{3}),'fontsize',15,'fontweight','bold','rotation',90);
    end
    patch([0.5 1.5 1.5 0.5],[Trials_sum(1,p)+0.5 Trials_sum(1,p)+0.5 Y_step(1,p)+0.5 Y_step(1,p)+0.5],'w','edgecolor','k');
    patch([0.5 1.5 1.5 0.5],[Y_step(2,p)+0.5 Y_step(2,p)+0.5 Y_step(3,p)+0.5 Y_step(3,p)+0.5],'w','edgecolor','k');            
    for nn = 1:3
        CaS_temp_4 = find(TrialType(p+(nn-1)*length(Fre),:)==1);
        if nn == 1
            imagesc(1,[1 Trials_sum(1,p),],CaS_temp_4');
        elseif nn == 2
            imagesc(1,[Y_step(1,p)+1:Y_step(2,p)],CaS_temp_4');
        elseif nn == 3
            imagesc(1,[Y_step(3,p)+1:Y_step(4,p)],CaS_temp_4'); 
        end
    end           
    ylim([0.5 Y_step(4,p)+0.5]);
    xlim([0.5 1.5]);       
    set(gca,'xtick',[0:1000/FrameTime:7000/FrameTime],'XTickLabel',[0:1:7],'ytick',Y_step(4,p),'YTickLabel',sum(Trials_sum(:,p)),'clim',[1 TrialN]);
    xlabel([num2str(round(Fre(p)/1000))]);
end    
colormap hsv;
suptitle('Trial Sequency');

for nROI = 1:nROIs 
    if savePPT == 1
        exportToPPTX('switchslide',nROI);
        exportToPPTX('addpicture',fig1,'position',[0 5 2 2]);       
    end    
end

if saveFig == 1
    saveas(fig1,[PathName '\' fn([8:end]) '-TrialSeq-' num2str(nROI,'%03d') '.fig']);
end
close;

exportToPPTX('saveandclose');

