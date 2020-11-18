%%
clear
close all
clc
% SavingFolder='G:\Data_FiberPhotometry\#FiberPhotometry_DataPlot\plot_2AFCdayfresh';
SavingFolder='G:\Data_FiberPhotometry\#FiberPhotometry_DataPlot\plot_MeanTrace_StimOct';

DatePath=[];

% DatePath{1}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\FI07\tsh_FI07_FP_20180605_anterior_beh_1';
% DatePath{2}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\FI07\tsh_FI07_FP_20180524_posterior_beh_1';
% DatePath{3}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\FI08\tsh_FI08_FP_20180602_posterior_beh_1';
% DatePath{4}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\FI08\tsh_FI08_FP_20180605_anteriorL_beh_1';
% DatePath{5}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\FI09\tsh_FI09_FP_20180524_posterior_beh_1';
% DatePath{6}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\FI09\tsh_FI09_FP_20180601_anterior_beh_1';
% DatePath{7}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\FI06\tsh_FI06_20180208_10uW_10ms_1';
% DatePath{8}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\FI02\tsh_FI02_20180208_10uW_10ms_1';
% DatePath{9}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\FI02\tsh_FI02_20180206_10uW_10ms_1';
% DatePath{10}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\FI03\tsh_FI03_20180209_10uW_10ms_1';
% DatePath{11}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\FI01\tsh_FI01_20180208_10uW_10ms_1';
% DatePath{12}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\FI01\tsh_FI01_20180206_10uW_10ms_2';
% DatePath{13}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\FI06\tsh_FI06_20180206_10uW_10ms_1';
% DatePath{14}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\FI04\tsh_FI04_20180207_10uW_10ms_1';

% DatePath{1}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\F10\tsh_F10_FP_20180601_posterior_beh_1';
% DatePath{2}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\F10\tsh_F10_FP_20180604_anterior_beh_1';
% DatePath{3}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\F07\tsh_F07_FP_20180601_posterior_beh_1';
% DatePath{4}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\F07\tsh_F07_FP_20180604_anterior_beh_1';
% DatePath{5}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\F06\tsh_F06_FP_20180523_posterior_beh_1';
% DatePath{6}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\F06\tsh_F06_FP_20180530_anterior_beh_1';
% DatePath{7}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\F09\tsh_F09_FP_20180523_posterior_beh_1';
% DatePath{8}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\F09\tsh_F09_FP_20180530_anterior_beh_1';

DatePath{1}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\F10\tsh_F10_FP_20180601_posterior_beh_1';
DatePath{2}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\F10\tsh_F10_FP_20180604_anterior_beh_1';
DatePath{3}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\F07\tsh_F07_FP_20180601_posterior_beh_1';
DatePath{4}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\F07\tsh_F07_FP_20180604_anterior_beh_1';
DatePath{5}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\F06\tsh_F06_FP_20180523_posterior_beh_1';
DatePath{6}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\F06\tsh_F06_FP_20180530_anterior_beh_1';
DatePath{7}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\F09\tsh_F09_FP_20180523_posterior_beh_1';
DatePath{8}='G:\Data_FiberPhotometry\#FiberPhotometry_DataCooked\F09\tsh_F09_FP_20180530_anterior_beh_1';


%%
% the session to be plotted
for n_path=1:length(DatePath)
    clearvars -except n_path SavingFolder DatePath


%
cd(DatePath{n_path});

MatFiles = dir('*.mat');
FrameInfo = dir('*.log');
for n_mat=1:length(MatFiles)
load(MatFiles(n_mat).name);
end

% load BehFrame txt/mat (lastest data format)
[TrialCount,TrialStart_FrameCount]=textread(FrameInfo.name,'%d %d','headerlines',12);

% ImagingSetup
ImagingSetup=FrameInfo.name(1:end-4);
disp(ImagingSetup);
ImagingSetup_L=sprintf('%s_Left', ImagingSetup);
ImagingSetup_R=sprintf('%s_Right', ImagingSetup);
ImagingSetup_Ctrl=sprintf('%s_410Ctrl', ImagingSetup);

Session470LEDstartFrame=1;% usually 205
Session410LEDstartFrame=2;
nTrial=length(TrialCount);
cd(SavingFolder);

%% ImagingArduino BehFrame Data modification
% for previous data format
%{
TrialStart_Count=size(BehFrame,1);
TrialStart_FrameCount=zeros(1,TrialStart_Count);
TrialStart_FrameCount(1:end)=BehFrame(1:end,2);
% TrialStart_FrameCount(432:end)=BehFrame(432:end,2);

% TrialStart_FrameCount=BehFrame(:,2);
TrialStart_FrameTime=zeros(TrialStart_Count,1);
TrialStart_FrameTime(1:end)=BehFrame(1:end,3);
% TrialStart_FrameTime(432:end)=BehFrame(432:end,3);
% TrialStart_FrameTime=BehFrame(:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nTrial=TrialStart_Count;

%%
TrialStart_TimeInterval=(TrialStart_FrameTime(2:end)-TrialStart_FrameTime(1:end-1))';
TrialFrameNum=double(TrialStart_FrameCount(2:end)-TrialStart_FrameCount(1:end-1));
FrameTimeAll=TrialStart_TimeInterval(1:end)./TrialFrameNum;
figure
histogram(FrameTimeAll);

temp=polyfit(1:length(FrameTimeAll),FrameTimeAll,1);
FrameTime=temp(2);
FrameRate=1000/FrameTime;
%}

% for lastest data format
nTrial = length(TrialCount);
FrameRate=40;
FrameTime=1000/FrameRate;
TrialFrameNum=double(TrialStart_FrameCount(2:end)-TrialStart_FrameCount(1:end-1));

%%


if exist('ROIdata_all')==0
Data_Left=(ROIdata(:,1));
Data_Right=(ROIdata(:,2));
else
    Data_Left=(ROIdata_all(:,1));
Data_Right=(ROIdata_all(:,2));
end




% DarkFrame_L=Data_Left(1:200);
% DarkFrame_R=Data_Right(1:200);

Data470_L_temp=Data_Left(1:2:end);
Data410_L_temp=Data_Left(2:2:end);

figure
subplot(2,1,1)
histogram(Data470_L_temp); 
title('470','FontSize',18);
subplot(2,1,2)
histogram(Data410_L_temp);
title('410','FontSize',18);

SessionLEDstartFrame=max(Session470LEDstartFrame,Session410LEDstartFrame);

    FP_equalEnd=(floor(length(Data_Left(SessionLEDstartFrame:end))/2))*2;
    Data470_L=Data_Left(Session470LEDstartFrame:2:FP_equalEnd);
    Data410_L=Data_Left(Session410LEDstartFrame:2:FP_equalEnd);
    Data470_R=Data_Right(Session470LEDstartFrame:2:FP_equalEnd);
    Data410_R=Data_Right(Session410LEDstartFrame:2:FP_equalEnd);


    
%% after polyfit, then baseline correction, offset
%{
fit_index_410_L=polyfit(Data410_L,Data470_L,1);
fited_410_L=polyval(fit_index_410_L,Data410_L);
figure
plot(Data470_L);
hold on;
plot(fited_410_L);

Data470_fitleveled_L=baseline_correction(Data470_L,500);
Data410_fitleveled_L=baseline_correction(fited_410_L,500);
figure
plot(Data470_fitleveled_L);
hold on;
plot(Data410_fitleveled_L);
%%%%%%%%%%%%%%%%%%%%%%%%

fit_index_410_R=polyfit(Data410_R,Data470_R,1);
fited_410_R=polyval(fit_index_410_R,Data410_R);
figure
plot(Data470_R);
hold on;
plot(fited_410_R); 

Data470_fitleveled_R=baseline_correction(Data470_R,500);
Data410_fitleveled_R=baseline_correction(fited_410_R,500);
figure
plot(Data470_fitleveled_R);
hold on;
plot(Data410_fitleveled_R);
    
    
Data470_offseted_L=Data470_fitleveled_L-min(min(Data470_fitleveled_L),min(Data410_fitleveled_L));
Data410_offseted_L=Data410_fitleveled_L-min(min(Data470_fitleveled_L),min(Data410_fitleveled_L));

Data470_offseted_R=Data470_fitleveled_R-min(min(Data470_fitleveled_R),min(Data410_fitleveled_R));
Data410_offseted_R=Data410_fitleveled_R-min(min(Data470_fitleveled_R),min(Data410_fitleveled_R));
%}
%% no polyfit, directly baseline correction, offset  
% %{
Data470_leveled_L=baseline_correction_leveled(Data470_L,500);
Data410_leveled_L=baseline_correction_leveled(Data410_L,500);
Data470_leveled_R=baseline_correction_leveled(Data470_R,500);
Data410_leveled_R=baseline_correction_leveled(Data410_R,500);

figure
plot(Data470_leveled_L);
hold on;
plot(Data410_leveled_L);
figure
histogram(Data470_leveled_L);
hold on;
histogram(Data410_leveled_L);

figure
plot(Data470_leveled_R);
hold on;
plot(Data410_leveled_R);
figure
histogram(Data470_leveled_R);
hold on;
histogram(Data410_leveled_R);

% Data470_offseted_L=Data470_leveled_L-min(min(Data470_leveled_L),min(Data410_leveled_L));
% Data410_offseted_L=Data410_leveled_L-min(min(Data470_leveled_L),min(Data410_leveled_L));
% 
% Data470_offseted_R=Data470_leveled_R-min(min(Data470_leveled_R),min(Data410_leveled_R));
% Data410_offseted_R=Data410_leveled_R-min(min(Data470_leveled_R),min(Data410_leveled_R));
%}
%%
%{
FrameMin=min(length(Data470_offseted_L),length(Data410_offseted_L));
df_L=Data470_offseted_L(1:FrameMin)-Data410_offseted_L(1:FrameMin);
dff_L=(df_L./Data410_offseted_L(1:FrameMin))*100;
figure
plot(dff_L);
title('dff-L');
Data_dffLeft=zeros(size(Data_Left));
Data_dffLeft(Session470LEDstartFrame:2:FP_equalEnd)=dff_L;
Data_dffLeft(Session410LEDstartFrame:2:FP_equalEnd)=Data410_offseted_L;


FrameMin=min(length(Data470_offseted_R),length(Data410_offseted_R));
df_R=Data470_offseted_R(1:FrameMin)-Data410_offseted_R(1:FrameMin);
dff_R=(df_R./Data410_offseted_R(1:FrameMin))*100;
figure
plot(dff_R);
title('dff-R');
Data_dffRight=zeros(size(Data_Right));

Data_dffRight(Session470LEDstartFrame:2:FP_equalEnd)=dff_R;
Data_dffRight(Session410LEDstartFrame:2:FP_equalEnd)=Data410_offseted_R;
%}
%% 470 mode as F baseline
FrameMin=min(length(Data470_leveled_L),length(Data410_leveled_L));
F_MotionCorr_L=Data470_leveled_L(1:FrameMin)-Data410_leveled_L(1:FrameMin);
F_MotionCorr_L_offseted=F_MotionCorr_L-min(F_MotionCorr_L);
% dff_L=(F_MotionCorr_L./Data410_leveled_L(1:FrameMin))*100;

dff_L_mode=((F_MotionCorr_L_offseted-Mode_tsh(F_MotionCorr_L_offseted,500))/Mode_tsh(F_MotionCorr_L_offseted,500))*100;
figure
% plot(dff_L);
Data_dffLeft=zeros(size(Data_Left));
Data_dffLeft(Session470LEDstartFrame:2:FP_equalEnd)=dff_L_mode;
Data_dffLeft(Session410LEDstartFrame:2:FP_equalEnd)=Data410_leveled_L;
hold on
plot(dff_L_mode);
title('dff-L');

% legend('mode');

FrameMin=min(length(Data470_leveled_R),length(Data410_leveled_R));
F_MotionCorr_R=Data470_leveled_R(1:FrameMin)-Data410_leveled_R(1:FrameMin);
F_MotionCorr_R_offseted=F_MotionCorr_R-min(F_MotionCorr_R);
% dff_R=(F_MotionCorr_R./Data410_leveled_R(1:FrameMin))*100;
dff_R_mode=((F_MotionCorr_R_offseted-Mode_tsh(F_MotionCorr_R_offseted,500))/Mode_tsh(F_MotionCorr_R_offseted,500))*100;

figure
% plot(dff_R);
Data_dffRight=zeros(size(Data_Right));

Data_dffRight(Session470LEDstartFrame:2:FP_equalEnd)=dff_R_mode;
Data_dffRight(Session410LEDstartFrame:2:FP_equalEnd)=Data410_leveled_R;
hold on
plot(dff_R_mode);
title('dff-R');

% legend('410','mode');
%% load Beh mat data
Trial_FP_startTime=cellfun(@(x) x.Time_trialLED_On,SessionResults);
StimOnset_Time= cellfun(@(x) x.Time_stimOnset, SessionResults)-Trial_FP_startTime;
Action_Choice = double(cellfun(@(x) x.Action_choice,SessionResults));
Trial_StimType = cellfun(@(x) x.Trial_Type, SessionResults);
AnswerTime=cellfun(@(x) x.Time_answer,SessionResults)-Trial_FP_startTime;
% BehTrialFrameCount= cellfun(@(x) x.Trial_FrameCount, SessionResults);
Trial_StimFreq = cellfun(@(x) x.Stim_toneFreq, SessionResults);
Trial_Inds= cellfun(@(x) x.Trial_inds, SessionResults);
TrialNum = length(Trial_Inds);%Trial num from behavior file
if TrialNum ~= nTrial
    warning('behavior trial num not equal to imaging trial num');
end

%%
Time_trialStart_ab=cellfun(@(x) x.Time_trialStart , SessionResults,'UniformOutput',false);

for j = 1: TrialNum  % abstract lick time
    temp_m = regexp(cell2mat(Time_trialStart_ab(j)), '\_', 'split');
    temp_n = [];
    for k = 1:length(temp_m)
        if isempty(temp_m(k))
            %             temp3 = temp3;
        else
            temp_n = [temp_n str2num(temp_m{k})];
        end
    end

    trial_start_ab_Time{j} = temp_n;
end

trial_start_ab_Time_millis=cellfun(@(x) x(4)*3600000+x(5)*60000+x(6)*1000+x(7),trial_start_ab_Time);
TrialLength=trial_start_ab_Time_millis(2:end)-trial_start_ab_Time_millis(1:end-1);
trialFrameLength=(TrialLength*FrameRate)/1000;
% FrameResidue=TrialFrameNum-trialFrameLength(1:end-1)';
%%
if max(StimOnset_Time)>100
    StimOnset_Frame = round((double(StimOnset_Time)/1000)*FrameRate);
elseif max(StimOnset_Time)<2.5
    StimOnset_Frame = round(StimOnset_Time*FrameRate);
end

if max(AnswerTime) > 100
    AnswerFrame = round((double(AnswerTime)/1000)*FrameRate);
else
    AnswerFrame = round(AnswerTime*FrameRate);
end
%%
FramePlotBeforeStim=min(StimOnset_Frame);
FramePlotAfterStim=min(TrialFrameNum'-StimOnset_Frame(1:end-1));
nFrame=floor((FramePlotBeforeStim+FramePlotAfterStim)/2)*2;
StimOnset_FrameCount=(TrialStart_FrameCount'+StimOnset_Frame)';
TrialFrameCount_AlignToSound=[];
TrialFrameCount_AlignToSound(:,1)=StimOnset_FrameCount-FramePlotBeforeStim;
TrialFrameCount_AlignToSound(:,2)=StimOnset_FrameCount+nFrame-FramePlotBeforeStim;
TrialData_L=zeros(TrialNum,nFrame);
TrialData_R=zeros(TrialNum,nFrame); 

Trial_dff_L=zeros(TrialNum,nFrame/2);
TrialData410_L=zeros(TrialNum,nFrame/2);
Trial_dff_R=zeros(TrialNum,nFrame/2);
TrialData410_R=zeros(TrialNum,nFrame/2);
%%
for n_trial=1:TrialNum
    TrialData_L(n_trial,:)=Data_dffLeft(TrialFrameCount_AlignToSound(n_trial,1)+1:TrialFrameCount_AlignToSound(n_trial,2));
    TrialData_R(n_trial,:)=Data_dffRight(TrialFrameCount_AlignToSound(n_trial,1)+1:TrialFrameCount_AlignToSound(n_trial,2));
    
    if mod(StimOnset_FrameCount(n_trial),2)==1
        Trial_dff_L(n_trial,:)=TrialData_L(n_trial,1:2:end);
        TrialData410_L(n_trial,:)=TrialData_L(n_trial,2:2:end);
        
        Trial_dff_R(n_trial,:)=TrialData_R(n_trial,1:2:end);
        TrialData410_R(n_trial,:)=TrialData_R(n_trial,2:2:end);
    end
    if mod(StimOnset_FrameCount(n_trial),2)==0
        Trial_dff_L(n_trial,:)=TrialData_L(n_trial,2:2:end);
        TrialData410_L(n_trial,:)=TrialData_L(n_trial,1:2:end);
        
        Trial_dff_R(n_trial,:)=TrialData_R(n_trial,2:2:end);
        TrialData410_R(n_trial,:)=TrialData_R(n_trial,1:2:end);
        
    end
end

Im_L=Trial_dff_L(:);
clim_L=[];
clim_L(1)=prctile(Im_L,10);
clim_L(2)=prctile(Im_L,90);
Im_R=Trial_dff_R(:);
clim_R=[];
clim_R(1)=prctile(Im_R,10);
clim_R(2)=prctile(Im_R,90);

figure
subplot(2,1,1);
imagesc(Trial_dff_L,clim_L);
colorbar;
hold on;
subplot(2,1,2);
imagesc(Trial_dff_R,clim_R);
colorbar;
%%
close
aa=11;
figure
subplot(2,1,1);
plot(Trial_dff_L(aa:aa,:)');
hold on;
subplot(2,1,2);
plot(Trial_dff_R(aa:aa,:)');
title(['Trial-' num2str(aa)]);
%%
close
figure
subplot(2,1,1);
plot(Trial_dff_L');
hold on;
subplot(2,1,2);
plot(Trial_dff_R');
%%

figure

%%
Im_L=TrialData410_L(:);
clim_L=[];
clim_L(1)=prctile(Im_L,10);
clim_L(2)=prctile(Im_L,90);
Im_R=TrialData410_R(:);
clim_R=[];
clim_R(1)=prctile(Im_R,10);
clim_R(2)=prctile(Im_R,90);

figure
subplot(2,1,1);
imagesc(TrialData410_L,clim_L);
colorbar;
hold on;
subplot(2,1,2);
imagesc(TrialData410_R,clim_R);
colorbar;
%%
close
figure
subplot(2,1,1);
plot(TrialData410_L');
hold on;
subplot(2,1,2);
plot(TrialData410_R');
% %%
%     Trial_ImagingTime = nFrame*FrameTime;
%     PSNL=0;
%     LickingRasterPlot_SortedByAnswerTime(SessionResults, SavingFolder,Trial_ImagingTime, PSNL,ImagingSetup)

% %%
% StimOnset_TimeAlighed=ones(TrialNum,1)*(double(min(StimOnset_Time)));
% AnswerTimeAfterSoundOnset=AnswerTime-StimOnset_Time+min(StimOnset_Time);
% FrameRate_UniChannel=FrameRate/2;
% 
% ROIdata470=zeros(size(Trial_dff_L,1),2,size(Trial_dff_L,2));
ROIdata470=[];
ROIdata470(:,1,:)=Trial_dff_L;%Corresponded to the ROI colorplot annotation 
ROIdata470(:,2,:)=Trial_dff_R;
% colorplot_FP(ROIdata470,Trial_StimFreq,StimOnset_TimeAlighed,Action_Choice,FrameRate_UniChannel,AnswerTimeAfterSoundOnset,SavingFolder,ImagingSetup);% align to sound &sort by answer time
% % 
% %%
% ROIdata410=zeros(size(TrialData410_L,1),2,size(TrialData410_L,2));
% ROIdata410(:,1,:)=TrialData410_L;%Corresponded to the ROI colorplot annotation
% ROIdata410(:,2,:)=TrialData410_R;
% colorplot_FP(ROIdata410,Trial_StimFreq,StimOnset_TimeAlighed,Action_Choice,FrameRate_UniChannel,AnswerTimeAfterSoundOnset,SavingFolder,ImagingSetup_Ctrl);% align to sound &sort by answer time
% 
% % %%
% % ROIdata470(:,1,:)=Trial_dff_L;%Corresponded to the ROI colorplot annotation
% % ROIdata470(:,2,:)=Trial_dff_R;
% colorplot_MeanTrace_FP(ROIdata470,Trial_StimFreq,StimOnset_TimeAlighed,Action_Choice,FrameRate_UniChannel,AnswerTimeAfterSoundOnset,SavingFolder,ImagingSetup);% align to sound &sort by answer time
%%
TrialFrameNum=[TrialFrameNum;Mode_tsh(TrialFrameNum,50)]; 

close all

% Colorplot_AlighToSound_MeanTrace_Anova_LickingRaster_FP(ROIdata470, SessionResults, FrameRate/2, SavingFolder,ImagingSetup);
% ImagingHemi=0;
% h_all=plot_FP_MeanTrace_1stLickAfterAns_ANOVAN_PSNL(Data_dffLeft,TrialFrameNum,TrialStart_FrameCount,SessionResults, FrameRate/2, SavingFolder, ImagingSetup,DatePath{n_path},ImagingHemi);
ImagingHemi=0;
ROIdata470=[];
ROIdata470(:,1,:)=Trial_dff_L;
saveMeanAroundAns_FP(ROIdata470,TrialFrameNum,TrialStart_FrameCount,SessionResults, FrameRate/2, SavingFolder, ImagingSetup,DatePath{n_path},ImagingHemi);
ImagingHemi=1;
ROIdata470=[];
ROIdata470(:,1,:)=Trial_dff_R;
saveMeanAroundAns_FP(ROIdata470,TrialFrameNum,TrialStart_FrameCount,SessionResults, FrameRate/2, SavingFolder, ImagingSetup,DatePath{n_path},ImagingHemi);

% %%
%  
% aa=ones(1,nTrial)*50;
% figure
% plot(Data_dffLeft)
% hold on
% sz=ones(1,nTrial)*20;
% scatter((StimOnset_FrameCount-206),aa,sz);
% % scatter(BehTrialFrameCount,aa,sz);
% %%
% clf
% aa=10;
% % figure
% subplot(4,1,1);
% plot(Trial_dff_R(aa,:)');
% line([46 46],[-50 100],'color','r')
% 
% hold on;
% subplot(4,1,2);
% plot(Trial_dff_R(aa+1,:)');
% line([46 46],[-50 100],'color','r')
% 
% subplot(4,1,3);
% plot(Trial_dff_R(aa+2,:)');
% line([46 46],[-50 100],'color','r')
% 
% subplot(4,1,4);
% plot(Trial_dff_R(aa+3,:)');
% 

% line([46 46],[-50 100],'color','r')


close all
% disp('Session Plot Done!')
disp(ImagingSetup);
end
