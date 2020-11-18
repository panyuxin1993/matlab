%%
clear
close all
clc
SavingFolder='F:\FP\pyx297_20200708\pyx297_20200708_mat';
DatePath=[];
DatePath{1}='F:\FP\pyx297_20200708\pyx297_20200708_mat';

%%
% the session to be plotted
for n_path=1:length(DatePath)
    clearvars -except n_path SavingFolder DatePath behaviorFile
    %
    cd(DatePath{n_path});
    MatFiles = dir('*MMStack*All.mat');

    for n_mat=1:length(MatFiles)
        load(MatFiles(n_mat).name);
   
        Session470LEDstartFrame=1;% usually 205
        Session410LEDstartFrame=2;
        Session560LEDstartFrame=3;
        cd(SavingFolder);
        
        
        if exist('ROIdata_all')==0 && size(ROIdata_all,2)==2
            Data_Left=(ROIdata(:,1));
            Data_Right=(ROIdata(:,2));
        elseif size(ROIdata_all,2)==2
            Data_Left=(ROIdata_all(:,1));
            Data_Right=(ROIdata_all(:,2));
            fiberstr={'Left','Right'};
        elseif size(ROIdata_all,2)==3
            Data_Left=(ROIdata_all(:,1));
            Data_Right=(ROIdata_all(:,2));
            Data_peri=(ROIdata_all(:,3));
            fiberstr={'Left','Right','Peri'};
        elseif size(ROIdata_all,2)==4
            Data_Left=(ROIdata_all(:,1));
            Data_Right=(ROIdata_all(:,2));
            Data_peri=(ROIdata_all(:,3));
            Data_stim=(ROIdata_all(:,4));
            fiberstr={'Left','Right','Peri','Stimuli'};
        end
      
        DarkFrame_L=Data_Left(1:200);
        DarkFrame_R=Data_Right(1:200);

        wavelengthstr={'470','410'};
        a=figure;   %plot fluorescence distribution
        set(gcf,'Position',[0,0,600,600]);
        for nfiber=1:size(ROIdata_all,2)
            Data470_temp=ROIdata_all(1:2:end,nfiber);
            Data410_ctrl=ROIdata_all(2:2:end,nfiber);
            for i=1:2
                subplot(2,size(ROIdata_all,2),2*i+nfiber-2);
                histogram(ROIdata_all(i:2:end,nfiber));
                title([fiberstr{nfiber},wavelengthstr{i}],'FontSize',18);
                if i==2
                    xl=xlim;
                elseif i==4
                    xlim(xl);
                end           
                set(gca,'FontName','Arial','FontSize',14);
            end
%             subplot(5,2,8+nfiber);
%             boxplot([ROIdata(2:4:end,nfiber),ROIdata(4:4:end,nfiber)],'Notch','on','Labels',{'410 control','410 test'});
            
        end
        saveas(gcf,['histogram of ',MatFiles(n_mat).name,'.fig'],'fig');
        b=figure;%plot raw trace
        set(gcf,'Position',[0,0,600,400]);
        for nfiber=1:size(ROIdata_all,2)
            subplot(1,size(ROIdata_all,2),nfiber);
            for i=1:size(ROIdata_all,2)
                plot(ROIdata_all(i:2:end,nfiber));
                hold on;
            end
        end
        saveas(gcf,['raw trace of ',MatFiles(n_mat).name,'.fig'],'fig');
    end
end