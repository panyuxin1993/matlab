function [  ] = AlignedCaS_fun( dff_filename,DataType )
%UNTITLED Summary of this function goes here
%   dff_filename: use dff data
%   DataType:  % 0:dff_modeF0  1:dff_Modef0SubPreTone  2:dff_meanF0
    if ~exist('New_aligned_data')
        mkdir('New_aligned_data');
    end
    cd 'New_aligned_data';
    savepath = pwd;
    cd ..
    fn = dff_filename;       %#####################################################
    SaveName = fn(1:end-4);
    load(fn);
    align_time = FirstLickAndRewTime_fun(Data_extract,1);

    OnsetFrame = round(align_time.onsetTime/FrameTime);
    FirstLickFrame = round(align_time.FirstLickTime(2,:)/FrameTime); % original frame
    RewardFrame = round(align_time.RewardT(1,:)/FrameTime); % original frame

    MinOnsFrame = min(OnsetFrame);
    MaxOnsFrame = max(OnsetFrame);
    MinActFrame = min(FirstLickFrame);
    MaxActFrame = max(FirstLickFrame);
    MinRewFrame = min(RewardFrame(align_time.OutcomeType==1));
    MaxRewFrame = max(RewardFrame);
    Frames = size(dff_modeF0,3);
    for triN = 1:length(OnsetFrame) 
       F_num_T(triN,:) = OnsetFrame(triN) + [-MinOnsFrame+1 Frames-MaxOnsFrame];     
       if ~isnan(FirstLickFrame(triN))
           F_num_A(triN,:) = FirstLickFrame(triN) + [-MinActFrame+1 Frames-MaxActFrame];           
       else
           F_num_A(triN,:) = [1 MinActFrame+Frames-MaxActFrame];    
       end
       if ~isnan(RewardFrame(triN)) & RewardFrame(triN)>0
           F_num_R(triN,:) = RewardFrame(triN) + [-MinRewFrame+1 Frames-MaxRewFrame];
       else
           F_num_R(triN,:) = [1 MinRewFrame+Frames-MaxRewFrame]; 
       end
    end
    if DataType == 0
        temp = dff_modeF0;   % 0:dff_modeF0  1:dff_Modef0SubPreTone  2:dff_meanF0   #####################################################
        SavingName = ['AlineCaSig_DfmodeF0_' SaveName];
    elseif DataType == 1
        temp = dff_Modef0SubPreTone;
        SavingName = ['AlineCaSig_DfModef0SubPreTone_' SaveName];
    elseif DataType == 2
        temp = dff_meanF0;
        SavingName = ['AlineCaSig_DfMeanF0_' SaveName];
    end

    length_frames_A = MinActFrame+Frames-MaxActFrame;
    length_frames_R = MinRewFrame+Frames-MaxRewFrame;

    for triN = 1:size(temp,1)
        CaSigAlineTone(triN,:,:) = temp(triN,:,[F_num_T(triN,1):F_num_T(triN,2)]);  % get the same frames before and after tone onset 
        CaSigAlineAct(triN,:,:) = temp(triN,:,[F_num_A(triN,1):F_num_A(triN,2)]);  % get the same frames before and after choice 
        CaSigAlineRew(triN,:,:) = temp(triN,:,[F_num_R(triN,1):F_num_R(triN,2)]);
    end    

    AlineCaSigData.CaSigAlineTone = CaSigAlineTone;
    AlineCaSigData.CaSigAlineAct = CaSigAlineAct;
    AlineCaSigData.CaSigAlineRew = CaSigAlineRew;

    save([savepath '\' SavingName],'AlineCaSigData','MinOnsFrame','MinActFrame','MinRewFrame','FrameTime','Trials_sum','Data_extract','TestTrialNum','align_time');

end

