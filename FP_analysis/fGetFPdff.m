function [ dff ] = fGetFPdff( rootpath )
%FGETFPDFF using '*All.mat' file to get dff for following analyses
%   Input is rootpath(also can define in function to use like script
%   output is dff, also will save as mat file.

close all
clc
if ~exist('rootpath','var')
    rootpath='F:\FP\pyx237_20191214';
end
DatePath=[];
DatePath{1}=rootpath;

% the session to be plotted
for n_path=1:length(DatePath)
    clearvars -except n_path SavingFolder DatePath behaviorFile
    %
    cd(DatePath{n_path});
    MatFiles = dir('*All.mat');
    for n_mat=1:length(MatFiles)
        load(MatFiles(n_mat).name);
    end
    
    
    Session470LEDstartFrame=1;% usually 205
    Session410LEDstartFrame=2;

    if exist('ROIdata_all','var')==0
        Data_Left=(ROIdata(:,1));
        Data_Right=(ROIdata(:,2));
    else
        Data_Left=(ROIdata_all(:,1));
        Data_Right=(ROIdata_all(:,2));
    end
    
    DarkFrame_L=Data_Left(1:200);
    DarkFrame_R=Data_Right(1:200);
    fiberstr={'Left','Right'};
    figure;   %plot fluorescence distribution
    for nfiber=1:2
        Data470_temp=ROIdata_all(1:2:end,nfiber);
        Data410_temp=ROIdata_all(2:2:end,nfiber);
        subplot(2,2,nfiber)
        histogram(Data470_temp);
        title([fiberstr{nfiber},'470'],'FontSize',18);
        xl=xlim;
        set(gca,'FontName','Arial','FontSize',14);
        subplot(2,2,nfiber+2)
        histogram(Data410_temp);
        title([fiberstr{nfiber},'410'],'FontSize',18);
        %xlim(xl);
        set(gca,'FontName','Arial','FontSize',14);
    end
    saveas(gcf,'Fluorescence_distribution.fig','fig');
    
    %movement correction and baseline correction
    SessionLEDstartFrame=max(Session470LEDstartFrame,Session410LEDstartFrame);
    FP_equalEnd=(floor(length(Data_Left(SessionLEDstartFrame:end))/2))*2;
    Data470_L=Data_Left(Session470LEDstartFrame:2:FP_equalEnd);
    Data410_L=Data_Left(Session410LEDstartFrame:2:FP_equalEnd);
    Data470_R=Data_Right(Session470LEDstartFrame:2:FP_equalEnd);
    Data410_R=Data_Right(Session410LEDstartFrame:2:FP_equalEnd);
    Data2fiber2wavelength=cell(2,2);
    Data2fiber2wavelength{1,1}=Data470_L';%transform to 1*n vector
    Data2fiber2wavelength{2,1}=Data410_L';
    Data2fiber2wavelength{1,2}=Data470_R';
    Data2fiber2wavelength{2,2}=Data410_R';
    frameNum=cellfun(@(x) length(x),Data2fiber2wavelength);
    FrameMin=min(min(frameNum));%determine the min frame number
    % after polyfit, then baseline correction, offset then dff OR baseline correction then polyfit for motion correction
    Data_offsetted=cell(2,2);%2fiber2wavelength
    Data_motionCor=zeros(2,FrameMin);%1d-two fiber, 2d-470signals that are motion corrected with 410fitted signal, then baseline correted and offset
    dff_1m2b=zeros(2,FrameMin);%1d-two fiber, 2d-dff; first motion corrected then baseline corrected
    dff_1b2m=zeros(2,FrameMin);%1d-two fiber, 2d-dff; first baseline corrected then motion corrected
    dff470=zeros(2,FrameMin);%1d-two fiber, 2d-dff; no motion correction
    dff410=zeros(2,FrameMin);
    fighist=figure;%plot hist
    for nfiber=1:2
        %method one for motion correction- fitting 410 with 470
        fit_index_410=polyfit(Data2fiber2wavelength{2,nfiber},Data2fiber2wavelength{1,nfiber},1);
        fited_410=polyval(fit_index_410,Data2fiber2wavelength{2,nfiber});
        motionCor=Data2fiber2wavelength{1,nfiber}-fited_410+mean(fited_410);%adding mean back to prevent data being negative
        %method two for motion correction- 470-410
        %         D=mode_distance(Data2fiber2wavelength{1,nfiber},Data2fiber2wavelength{2,nfiber});
        %         motionCor=Data2fiber2wavelength{1,nfiber}-Data2fiber2wavelength{2,nfiber}-D;
        Data470_fitleveled=baseline_correction(Data2fiber2wavelength{1,nfiber},500);
        Data410_fitleveled=baseline_correction(Data2fiber2wavelength{2,nfiber},500);
        Data_motionCor(nfiber,:)=baseline_correction(motionCor,500);
        %         Data470_offseted=Data470_fitleveled-min(min(Data470_fitleveled),min(Data410_fitleveled));%so all data are positive
        %         Data410_offseted=Data410_fitleveled-min(min(Data470_fitleveled),min(Data410_fitleveled));
        %         Data_motionCor(nfiber,:)=Data_motionCor(nfiber,:)-min(Data_motionCor(nfiber,:));
        Data_offsetted{1,nfiber}=Data470_fitleveled;% Data470_offseted;
        Data_offsetted{2,nfiber}=Data410_fitleveled;% Data410_offseted;
        % 470 mode as F baseline and calculate dff
        fit_index=polyfit(Data_offsetted{2,nfiber},Data_offsetted{1,nfiber},1);
        fited_410_afterCor=polyval(fit_index,Data_offsetted{2,nfiber});
        F_MotionCorr=Data_offsetted{1,nfiber}(1:FrameMin)-fited_410_afterCor(1:FrameMin)+mean(fited_410_afterCor(1:FrameMin));%add mean back to prevent data being negative
        %         D=mode_distance(Data_offsetted{1,nfiber},Data_offsetted{2,nfiber});
        %         F_MotionCorr=Data_offsetted{1,nfiber}-Data_offsetted{2,nfiber}-D;
        dff_1b2m(nfiber,:) = dff_mode(F_MotionCorr);        %get dff
        dff_1m2b(nfiber,:) = dff_mode(Data_motionCor(nfiber,:));%get dff
        dff470(nfiber,:)=dff_mode(Data470_fitleveled);
        dff410(nfiber,:)=dff_mode(Data410_fitleveled);
        figdff2=figure;
        plot(dff_1b2m(nfiber,:));
        hold on;
        plot(dff_1m2b(nfiber,:));
        plot(dff470(nfiber,:));
        plot(dff410(nfiber,:));
        title('dff-mode');
        legend('dff baseline correction first','dff motion correction first','dff470','dff410');
        saveas(gcf,['Dff comparison_of_',fiberstr{nfiber},'_fiber.fig'],'fig');
        figure(fighist);
        subplot(2,2,nfiber)
        histogram(Data470_fitleveled);
        title([fiberstr{nfiber},'470'],'FontSize',18);
        xl=xlim;
        set(gca,'FontName','Arial','FontSize',14);
        subplot(2,2,nfiber+2)
        histogram(Data410_fitleveled);
        title([fiberstr{nfiber},'410'],'FontSize',18);
        %xlim(xl);
        set(gca,'FontName','Arial','FontSize',14);
        if nfiber==2
            saveas(gcf,'Fluorescence_distribution_after_baseline_correction.fig','fig');
        end
        
        %plot 1m2b process
        fig1m2d=figure;
        subplot(2,3,1);
        plot(Data2fiber2wavelength{1,nfiber});
        xlabel('time(frames)');
        ylabel('fluorescence');
        title('470&410');
        hold on;
        plot(Data2fiber2wavelength{2,nfiber});
        legend('470','410');
        subplot(2,3,2);
        plot(Data2fiber2wavelength{1,nfiber});
        hold on;
        plot(fited_410);
        title('470&fitted410');
        xlabel('time(frames)');
        ylabel('fluorescence');
        legend('470','fitted410');
        subplot(2,3,3);
        plot(motionCor);
        title('Motion corrected 470');
        xlabel('time(frames)');
        ylabel('fluorescence');
        subplot(2,3,4);
        scatter(Data2fiber2wavelength{2,nfiber},Data2fiber2wavelength{1,nfiber});
        ylabel('470');
        xlabel('410');
        title('relationship of 410 and 470 signals');
        subplot(2,3,5);
        plot(Data_motionCor(nfiber,:));
        title('baseline correted');
        xlabel('time(frames)');
        ylabel('fluorescence');
        subplot(2,3,6);
        plot(dff_1m2b(nfiber,:));
        title('dff');
        xlabel('time(frames)');
        ylabel('fluorescence');
        saveas(gcf,['Fluorescence_of_',fiberstr{nfiber},'_fiber_1m2b.fig'],'fig');
        %plot 1b2m process
        fig1b2m=figure;
        subplot(2,3,1);
        plot(Data2fiber2wavelength{1,nfiber});
        xlabel('time(frames)');
        ylabel('fluorescence');
        title('470&410');
        hold on;
        plot(Data2fiber2wavelength{2,nfiber});
        legend('470','410');
        subplot(2,3,2);
        plot(Data470_fitleveled);
        hold on;
        plot(Data410_fitleveled);
        title('baseline corrected');
        xlabel('time(frames)');
        ylabel('fluorescence');
        legend('470','410');
        subplot(2,3,3);
        plot(Data470_fitleveled);
        hold on;
        plot(fited_410_afterCor);
        title('470&fitted410');
        xlabel('time(frames)');
        ylabel('fluorescence');
        legend('470','fitted 410');
        subplot(2,3,4);
        scatter(Data470_fitleveled,Data410_fitleveled);
        ylabel('470');
        xlabel('410');
        title('relationship of 410 and 470 signals after baseline correction');
        subplot(2,3,5);
        plot(F_MotionCorr);
        title('motion correted');
        xlabel('time(frames)');
        ylabel('fluorescence');
        subplot(2,3,6);
        plot(dff_1b2m(nfiber,:));
        title('dff');
        xlabel('time(frames)');
        ylabel('fluorescence');
        saveas(gcf,['Fluorescence_of_',fiberstr{nfiber},'_fiber_1b2m.fig'],'fig');
    end
    dff=cell(1,4);
    dff{1}=dff_1b2m;
    dff{2}=dff_1m2b;
    dff{3}=dff470;
    dff{4}=dff410;
    save('dff_temp','dff');
    save('raw data.mat','Data2fiber2wavelength');
    close all;
end
end
%% assistant function part
function f=baseline_correction(data, span)
totalFr=size(data,2);%how many frames
ind_x = 0;
x = [];
for i = 1:totalFr
    %     ind_x = ind_x + 1;
    ind1 = i- span;
    ind2 = i+span;
    if ind1 <= 0
        ind1 = 1;
    end
    if ind2 > totalFr
        ind2 = totalFr;
    end
    x(i) = prctile(data(1,ind1:ind2),5);
end
f = data - x + mean(x);
% B=figure;
% set(gcf, 'position', [0 0 1500 500]);
% plot(data,'k-');
% hold on;
% plot(x,'c-');
% plot(f,'g-');
% legend('f','moving f mean','f baseline correction');
% set(gca,'FontName','Arial','FontSize',14);

end
function dff=dff_mode(f)
[N,edges,bin] = histcounts(f,100);
f0 = edges(N == max(N));
if length(f0)==2
    f0=mean(f0);
end
dff = (f - f0)/f0;
end
function D=mode_distance(data1,data2)
[N,edges,bin] = histcounts(data1,100);
mode1 = edges(N == max(N));
if length(mode1)==2
    mode1=mean(mode1);
end
[N,edges,bin] = histcounts(data2,100);
mode2 = edges(N == max(N));
if length(mode2)==2
    mode2=mean(mode2);
end
D=mode1-mode2;
end
