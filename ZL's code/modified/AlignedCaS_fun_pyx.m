function [  ] = AlignedCaS_fun_pyx( dff_filename,DataType )
%UNTITLED Summary of this function goes here
%   dff_filename: use dff data
%   DataType:  % 0:dff_modeF0  1:dff_Modef0SubPreTone  2:dff_meanF0
    if ~exist('New_aligned_data')
        mkdir('New_aligned_data');
    end
    cd 'New_aligned_data';
    savepath = pwd;
    cd ..
    fn = 'D:\xulab\imaging_data\pyx083-20180522\im_data_reg_cpu\result_save\dff_pyx083_20180522_920nm_power50_zoom4x_dft.mat';%dff_filename;       %#####################################################
    SaveName = fn(1:end-4);
    load(fn);
    align_time = FirstLickAndRewTime_fun_pyx(Data_extract,1);

    OnsetFrame = round(align_time.onsetTime/FrameTime);
    FirstLickFrame = round(align_time.FirstLickTime(2,:)/FrameTime); % original frame
    RewardFrame = round(align_time.RewardT(1,:)/FrameTime); % original frame

    MinOnsFrame = min(OnsetFrame);
    MaxOnsFrame = max(OnsetFrame);
    MinActFrame = min(FirstLickFrame);
    MaxActFrame = max(FirstLickFrame);
    MinRewFrame = min(RewardFrame(align_time.OutcomeType==1));
    MaxRewFrame = max(RewardFrame);
    Frames = SavedCaTrials.nFrames; %size(dff_modeF0,3);
    stimTime=500;
    delay=900;
    drinkTime=2000;
    ITI=4000;
    maxFrameFromOnset=(stimTime+delay+drinkTime+ITI)/FrameTime;
    for triN = 1:length(OnsetFrame) %for each trials
       F_num_T(triN,:) = OnsetFrame(triN) + [-MinOnsFrame+1 maxFrameFromOnset];     
%        if ~isnan(FirstLickFrame(triN))
%            F_num_A(triN,:) = FirstLickFrame(triN) + [-MinActFrame+1 Frames-MaxActFrame];           
%        else
%            F_num_A(triN,:) = [1 MinActFrame+Frames-MaxActFrame];    
%        end
%        if ~isnan(RewardFrame(triN)) & RewardFrame(triN)>0
%            F_num_R(triN,:) = RewardFrame(triN) + [-MinRewFrame+1 Frames-MaxRewFrame];
%        else
%            F_num_R(triN,:) = [1 MinRewFrame+Frames-MaxRewFrame]; 
%        end
    end
%     if DataType == 0
        temp = dff_mode; %dff_modeF0;   % 0:dff_modeF0  1:dff_Modef0SubPreTone  2:dff_meanF0   #####################################################
        SavingName = ['AlineCaSig_DfmodeF0_' SaveName];
%     elseif DataType == 1
%         temp = dff_Modef0SubPreTone;
%         SavingName = ['AlineCaSig_DfModef0SubPreTone_' SaveName];
%     elseif DataType == 2
%         temp = dff_mean; %dff_meanF0;
%         SavingName = ['AlineCaSig_DfMeanF0_' SaveName];
%     end

    length_frames_A = MinActFrame+Frames-MaxActFrame;
    length_frames_R = MinRewFrame+Frames-MaxRewFrame;

    for triN = 1:length(OnsetFrame) 
        for nROIs=1:SavedCaTrials.nROIs
            CaSigAlineTone(triN,nROIs,:) = temp{nROIs}(triN,[F_num_T(triN,1):F_num_T(triN,2)]);  % get the same frames before and after tone onset
%         CaSigAlineAct(triN,:,:) = temp(triN,:,[F_num_A(triN,1):F_num_A(triN,2)]);  % get the same frames before and after choice 
%         CaSigAlineRew(triN,:,:) = temp(triN,:,[F_num_R(triN,1):F_num_R(triN,2)]);
        end
    end    

    AlineCaSigData.CaSigAlineTone = CaSigAlineTone;
%     AlineCaSigData.CaSigAlineAct = CaSigAlineAct;
%     AlineCaSigData.CaSigAlineRew = CaSigAlineRew;

%     save([savepath '\' SavingName],'AlineCaSigData','MinOnsFrame','MinActFrame','MinRewFrame','FrameTime','Trials_sum','Data_extract','TestTrialNum','align_time');
saveFileName=['D:\xulab\imaging_data\pyx083-20180522\im_data_reg_cpu\result_save\AlineCaSig_DfmodeF0_dff_pyx083_20180522_920nm_power50_zoom4x_dft.mat'];
save(saveFileName);
end
% function [align_time ]=FirstLickAndRewTime_fun(Data_extract)
% align_time=struct('onsetTime',[],'FirstLickTime',[],'RewardT',[],'OutcomeType',[]);
% align_time.onsetTime= Data_extract.Tone_onset_time;
% align_time.RewardT=Data_extract.Reward_time;
% align_time.OutcomeType= Data_extract.Reward_time>0;%zero when error/miss
% for i=1:length(Data_extract.Trial_number)
%     tempL=Data_extract.Left_lick_time{i}((Data_extract.Left_lick_time{i}-Data_extract.Tone_onset_time(i))>0);
%     align_time.LeftFistLickTime(i)=min(tempL);
%     tempR=Data_extract.Right_lick_time{i}((Data_extract.Right_lick_time{i}-Data_extract.Tone_onset_time(i))>0);
%     align_time.RightFistLickTime=min(tempR);
%     align_time.FirstLickTime=min(tempL,tempR);
%     if isempty(align_time.FirstLickTime)
%         align_time.FirstLickTime=nan;
%     end
% end
% end
