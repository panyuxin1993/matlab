% files = dir('*.mat');

% for ff = 1:length(files)
Alig2Tone = 1; % 1: aline to tone    0: aline to first lick time
LineW = 1;
%=================================================
% fn = files(ff).name;
% saveName = fn(12:end);
% load(fn);
load('D:\xulab\imaging_data\pyx083-20180522\im_data_reg_cpu\result_save\AlineCaSig_DfmodeF0_dff_pyx083_20180522_920nm_power50_zoom4x_dft.mat');
if Alig2Tone
    SeleData = AlineCaSigData.CaSigAlineTone;
    SeleOnsFrame = MinOnsFrame;
else
    SeleData = AlineCaSigData.CaSigAlineAct;
    SeleOnsFrame = MinActFrame;
end
[ CorrTrial,WroTrial,MisTrial] = TrialType_fun_pyx(Data_extract,[],1);
Col = {'k',[0.4 .4 .4],[0.6 0.6 0.6],[0.8 0.8 0.8],[1 .9 0],[1 .6 0],[1 .3 0],'r'};
%%
% if ~isempty(exportToPPTX('query'))
%     exportToPPTX('close');
% end
% if exist([saveName '.pptx'])
% %     exportToPPTX('open',saveName);
% continue;
% else
%     exportToPPTX('new','Dimensions',[16 9]);
% end

Fre = unique(Data_extract.Stimulus);
Lick_out = FirstLickAndRewTime_fun_pyx(Data_extract,1);
% [ FirstLickTime_out,RewardTime_out ] = FirstLickAndRewTime_fun(Data_extract,TestTrialNum);
FirstLickTime = round(Lick_out.FirstLickTime(3,:)/FrameTime) + double(MinOnsFrame);
RewardTime = round(Lick_out.RewardT(2,:)/FrameTime) + double(MinOnsFrame);
for nROI = 1:size(SeleData,2)
    ROI_data = squeeze(SeleData(:,nROI,:))*100; 
    temp=mean(ROI_data(:,1:MinOnsFrame),2);
    ROI_data= ROI_data-repmat(temp,1,size(ROI_data,2));
    fig1 = figure;set(gcf,'position',[1950 0 450 1100],'color','w'); %colormap jet; 
    fig2 = figure;set(gcf,'position',[2400 0 450 1100],'color','w'); %colormap jet; 
    fig3 = figure;set(gcf,'position',[2850 0 450 1100],'color','w'); %colormap jet; 
%     fig4 = figure;set(gcf,'position',[3300 0 1000 1100],'color','w'); 
%     fig5 = figure;set(gcf,'position',[3000 700 550 400],'color','w'); line([0 1000],[0 0],'linewidth',1,'color','k');
%     set(gca,'xtick',[SeleOnsFrame:round(1000/FrameTime):size(ROI_data,2)],'xticklabel',[0:10],'fontweight','bold'); title(['ROI-' num2str(nROI) ' Correct']);
%     fig6 = figure;set(gcf,'position',[3300 700 550 400],'color','w'); line([0 1000],[0 0],'linewidth',1,'color','k');
%     set(gca,'xtick',[SeleOnsFrame:round(1000/FrameTime):size(ROI_data,2)],'xticklabel',[0:10],'fontweight','bold'); title(['ROI-' num2str(nROI) ' Wrong']);
%     fig7 = figure;set(gcf,'position',[3300 700 550 400],'color','w'); line([0 1000],[0 0],'linewidth',1,'color','k');
%     set(gca,'xtick',[SeleOnsFrame:round(1000/FrameTime):size(ROI_data,2)],'xticklabel',[0:10],'fontweight','bold'); title(['ROI-' num2str(nROI) ' Miss']);    
%     
    

   ColLimit = prctile(ROI_data(:),90);   
   ColLimit_1 = prctile(ROI_data(:),20);
   ColLimit_2 = prctile(ROI_data(:),95);  
   for i = 1:size(CorrTrial,1)
       CorrTrialInd = CorrTrial(i,:);
       WroTrialInd = WroTrial(i,:);
       MisTrialInd = MisTrial(i,:);
       figure(fig1);   % color plot for correct trials       
       if i == 1
           subplot(size(CorrTrial,1)+2,1,[1:2]); 
       elseif i == size(CorrTrial,1)
           subplot(size(CorrTrial,1)+2,1,[size(CorrTrial,1)+1:size(CorrTrial,1)+2]);
       else
           subplot(size(CorrTrial,1)+2,1,i+1);
       end
       hold on;
       FT = FirstLickTime(CorrTrialInd);
       RT = RewardTime(CorrTrialInd);
       [B,I] = sort(FT);
       tempData = ROI_data(CorrTrialInd,:);
       imagesc(tempData(I,:));
       plot(B,[1:sum(CorrTrialInd)],'m.','markersize',8);
       plot(RT(I),[1:sum(CorrTrialInd)],'c.','markersize',8);
       line([SeleOnsFrame SeleOnsFrame],[0.5 sum(CorrTrialInd)+0.5],'linestyle','--','linewidth',LineW*2,'color','w');
%        line([SeleOnsFrame SeleOnsFrame]+round(500/FrameTime),[0.5 sum(CorrTrialInd)+0.5],'linestyle','--','linewidth',LineW*2,'color','w');
%        line([SeleOnsFrame SeleOnsFrame]+round(1400/FrameTime),[0.5 sum(CorrTrialInd)+0.5],'linestyle','-','linewidth',LineW*3,'color',[.7 .7 .7]);       
       set(gca,'clim',[0 ColLimit],'ytick',sum(CorrTrialInd),'yticklabel',sum(CorrTrialInd),'ydir','normal','xtick',[SeleOnsFrame:round(1000/FrameTime):size(ROI_data,2)],'xticklabel',[0:10]);
%        xlim([0 size(ROI_data,2)]);
       if sum(CorrTrialInd)>0
        ylim([0.5 sum(CorrTrialInd)+0.5]);
       end
       if i<size(CorrTrial,1)
           set(gca,'xcolor','none');
       end
       ylabel([num2str(Fre(i)) ' clicks/s']);
       if i == size(CorrTrial,1)
           xlabel('Time(s)');
           ColBar = colorbar;
           set(ColBar,'position',[0.91 0.1100 0.04 0.1445]);
       end
       if i == 1
          title(['ROI-' num2str(nROI) ' Correct']); 
       end
       figure(fig2);  % color plot for wrong trials
       subplot(size(CorrTrial,1),1,i);
      
       hold on;
       FT = FirstLickTime(WroTrialInd); 
       [B,I] = sort(FT);
       tempData = ROI_data(WroTrialInd,:);
       imagesc(tempData(I,:));
       plot(B,[1:sum(WroTrialInd)],'m.','markersize',8);
       line([SeleOnsFrame SeleOnsFrame],[0.5 sum(WroTrialInd)+0.5],'linestyle','--','linewidth',LineW*2,'color','w');
%        line([SeleOnsFrame SeleOnsFrame]+round(500/FrameTime),[0.5 sum(WroTrialInd)+0.5],'linestyle','--','linewidth',LineW*2,'color','w');
%        line([SeleOnsFrame SeleOnsFrame]+round(1400/FrameTime),[0.5 sum(WroTrialInd)+0.5],'linestyle','-','linewidth',LineW*3,'color',[.7 .7 .7]);           
       set(gca,'clim',[0 ColLimit],'ytick',sum(WroTrialInd),'yticklabel',sum(WroTrialInd),'ydir','normal','xtick',[SeleOnsFrame:round(1000/FrameTime):size(ROI_data,2)],'xticklabel',[0:10]);
       if sum(WroTrialInd)>0
       ylim([0.5 sum(WroTrialInd)+0.5]);
       end
       if i<size(CorrTrial,1)
           set(gca,'xcolor','none');
       end   
       ylabel([num2str(Fre(i)) ' clisks/s']);
       if i == size(CorrTrial,1)
           xlabel('Time(s)');
       end     
       if i == 1
          title(['ROI-' num2str(nROI) ' Wrong']); 
       end
       figure(fig3); hold on;  % color plot for miss trials
       subplot(size(CorrTrial,1),1,i);
       imagesc(ROI_data(MisTrialInd,:));     
       line([SeleOnsFrame SeleOnsFrame],[0.5 sum(MisTrialInd)+0.5],'linestyle','--','linewidth',LineW*2,'color','w');
%        line([SeleOnsFrame SeleOnsFrame]+round(300/FrameTime),[0.5 sum(MisTrialInd)+0.5],'linestyle','--','linewidth',LineW*2,'color','w');
%        line([SeleOnsFrame SeleOnsFrame]+round(800/FrameTime),[0.5 sum(MisTrialInd)+0.5],'linestyle','-','linewidth',LineW*3,'color',[.7 .7 .7]);       
       set(gca,'clim',[0 ColLimit],'ytick',sum(MisTrialInd),'yticklabel',sum(MisTrialInd),'ydir','normal','xtick',[SeleOnsFrame:round(1000/FrameTime):size(ROI_data,2)],'xticklabel',[0:10]);
       if sum(MisTrialInd)>0
       ylim([0.5 sum(MisTrialInd)+0.5]);
       end       
       if i<size(CorrTrial,1)
           set(gca,'xcolor','none');
       end   
       ylabel([num2str(Fre(i)) ' clisks/s']);
       if i == size(CorrTrial,1)
           xlabel('Time(s)');
       end
       if i == 1
          title(['ROI-' num2str(nROI) ' Miss']); 
       end       
       figure(fig4);   % Trace plot 
       MeanTrace_Corr = mean(ROI_data(CorrTrialInd,:),1);
       MeanTrace_Wro = mean(ROI_data(WroTrialInd,:),1);
       MeanTrace_Mis = mean(ROI_data(MisTrialInd,:),1);

       MeanTrace_Corr_sem = std(ROI_data(CorrTrialInd,:))/(sum(CorrTrialInd))^0.5;
       if sum(WroTrialInd)>0
           MeanTrace_Wro_sem = std(ROI_data(WroTrialInd,:))/(sum(WroTrialInd))^0.5;
       else
           MeanTrace_Wro_sem = [];
       end
       if sum(MisTrialInd)>0          
           MeanTrace_Mis_sem = std(ROI_data(MisTrialInd,:))/(sum(MisTrialInd))^0.5;
       else          
           MeanTrace_Mis_sem = [];
       end
       subplot(size(CorrTrial,1)/2,2,i);  hold on;
       line([0 size(ROI_data,2)],[0 0],'linewidth',0.5,'color','k');
       patch([1:length(MeanTrace_Corr) fliplr(1:length(MeanTrace_Corr))],[MeanTrace_Corr+MeanTrace_Corr_sem fliplr(MeanTrace_Corr-MeanTrace_Corr_sem)],'k','edgecolor','none');
       if sum(WroTrialInd)>0
            patch([1:length(MeanTrace_Corr) fliplr(1:length(MeanTrace_Corr))],[MeanTrace_Wro+MeanTrace_Wro_sem fliplr(MeanTrace_Wro-MeanTrace_Wro_sem)],'r','edgecolor','none');
       end
       if sum(MisTrialInd)>0
            patch([1:length(MeanTrace_Corr) fliplr(1:length(MeanTrace_Corr))],[MeanTrace_Mis+MeanTrace_Mis_sem fliplr(MeanTrace_Mis-MeanTrace_Mis_sem)],'c','edgecolor','none');
       end
       alpha 0.3;       
       plot(MeanTrace_Mis,'color','c','linewidth',1);
       plot(MeanTrace_Wro,'color','r','linewidth',1);
       plot(MeanTrace_Corr,'color','k','linewidth',1);

       line([SeleOnsFrame SeleOnsFrame],[ColLimit_1 ColLimit_2],'linestyle','--','linewidth',LineW,'color','k');
       line([SeleOnsFrame SeleOnsFrame]+round(300/FrameTime),[ColLimit_1 ColLimit_2],'linestyle','--','linewidth',LineW,'color','k');
       line([SeleOnsFrame SeleOnsFrame]+round(800/FrameTime),[ColLimit_1 ColLimit_2],'linestyle','-','linewidth',LineW*2,'color',[.7 .7 .7]);
       ylim([ColLimit_1 ColLimit_2]);
       xlim([1 size(ROI_data,2)]);
       if i<size(CorrTrial,1)-1 
           set(gca,'xcolor','none');
       end   
       set(gca,'xtick',[SeleOnsFrame:round(1000/FrameTime):size(ROI_data,2)],'xticklabel',[0:10]);
       ylabel([num2str(Fre(i)) ' kHz']);
       if i == size(CorrTrial,1)
           xlabel('Time(s)');
       end
       if i == 1
           text(5,ColLimit_2*0.9,'Correct','color','k');
           text(5,ColLimit_2*0.8,'Wrong','color','r');
           text(5,ColLimit_2*0.7,'Miss','color','c');
           title(['ROI-' num2str(nROI)]); 
       end   
%        figure(fig5);  %  correct trace
%        hold on;
% %        patch([1:length(MeanTrace_Corr) fliplr(1:length(MeanTrace_Corr))],[MeanTrace_Corr+MeanTrace_Corr_sem fliplr(MeanTrace_Corr-MeanTrace_Corr_sem)],Col{i},'edgecolor','none');
%        plot(MeanTrace_Corr,'color',Col{i},'linewidth',2.5);       
%        alpha 0.3;  
%        xlim([1 size(ROI_data,2)]);
%        ylim([ColLimit_1 ColLimit_2]);
%        if i == 8
%           xlabel('Time(s)');
%           ylabel('\DeltaF/F(%)');   
%        line([SeleOnsFrame SeleOnsFrame],[ColLimit_1 ColLimit_2],'linestyle','--','linewidth',LineW,'color','k');
%        line([SeleOnsFrame SeleOnsFrame]+round(300/FrameTime),[ColLimit_1 ColLimit_2],'linestyle','--','linewidth',LineW,'color','k');
%        line([SeleOnsFrame SeleOnsFrame]+round(800/FrameTime),[ColLimit_1 ColLimit_2],'linestyle','-','linewidth',LineW*2,'color',[.7 .7 .7]);          
%        end
%        text(5,ColLimit_2-ColLimit_2*0.05*i,[num2str(Fre(i)) ' kHz'],'color',Col{i});
%        figure(fig6);  %  wrong trace
%        hold on;
% %        patch([1:length(MeanTrace_Corr) fliplr(1:length(MeanTrace_Corr))],[MeanTrace_Corr+MeanTrace_Corr_sem fliplr(MeanTrace_Corr-MeanTrace_Corr_sem)],Col{i},'edgecolor','none');
%        plot(MeanTrace_Wro,'color',Col{i},'linewidth',2.5);       
%        alpha 0.3;  
%        xlim([1 size(ROI_data,2)]);
%        ylim([ColLimit_1 ColLimit_2]);
%        if i == 8
%           xlabel('Time(s)');
%           ylabel('\DeltaF/F(%)');   
%        line([SeleOnsFrame SeleOnsFrame],[ColLimit_1 ColLimit_2],'linestyle','--','linewidth',LineW,'color','k');
%        line([SeleOnsFrame SeleOnsFrame]+round(300/FrameTime),[ColLimit_1 ColLimit_2],'linestyle','--','linewidth',LineW,'color','k');
%        line([SeleOnsFrame SeleOnsFrame]+round(800/FrameTime),[ColLimit_1 ColLimit_2],'linestyle','-','linewidth',LineW*2,'color',[.7 .7 .7]);          
%        end
%        text(5,ColLimit_2-ColLimit_2*0.05*i,[num2str(Fre(i)) ' kHz'],'color',Col{i});      
%        figure(fig7);  %  Miss trace
%        hold on;
% %        patch([1:length(MeanTrace_Corr) fliplr(1:length(MeanTrace_Corr))],[MeanTrace_Corr+MeanTrace_Corr_sem fliplr(MeanTrace_Corr-MeanTrace_Corr_sem)],Col{i},'edgecolor','none');
%        plot(MeanTrace_Mis,'color',Col{i},'linewidth',2.5);       
%        alpha 0.3;  
%        xlim([1 size(ROI_data,2)]);
%        ylim([ColLimit_1 ColLimit_2]);
%        if i == 8
%           xlabel('Time(s)');
%           ylabel('\DeltaF/F(%)');   
%        line([SeleOnsFrame SeleOnsFrame],[ColLimit_1 ColLimit_2],'linestyle','--','linewidth',LineW,'color','k');
%        line([SeleOnsFrame SeleOnsFrame]+round(300/FrameTime),[ColLimit_1 ColLimit_2],'linestyle','--','linewidth',LineW,'color','k');      
%        end
%        text(5,ColLimit_2-ColLimit_2*0.05*i,[num2str(Fre(i)) ' kHz'],'color',Col{i});           
   end
%    exportToPPTX('addslide');
%    exportToPPTX('addpicture',fig1,'position',[0 0 3 7]);
%    exportToPPTX('addpicture',fig2,'position',[3 0 3 7]);
%    exportToPPTX('addpicture',fig3,'position',[6 0 3 7]);
%    exportToPPTX('addpicture',fig5,'position',[0 6.5 3 2.5]);
%    exportToPPTX('addpicture',fig6,'position',[3 6.5 3 2.5]);
%    exportToPPTX('addpicture',fig7,'position',[6 6.5 3 2.5]);
%    exportToPPTX('addpicture',fig4,'position',[8.5 0 7 7]);
%    exportToPPTX('addtext',['ROI-' num2str(nROI) '(' saveName ')'],[0 0.1]);


saveas(fig1,['D:\xulab\imaging_data\fig\ROI-',num2str(nROI),'-Cor.jpg'],'jpg');
saveas(fig2,['D:\xulab\imaging_data\fig\ROI-',num2str(nROI),'-Err.jpg'],'jpg');
saveas(fig3,['D:\xulab\imaging_data\fig\ROI-',num2str(nROI),'-Mis.jpg'],'jpg');
   close all;
end
%    exportToPPTX('saveandclose',saveName);
% end















