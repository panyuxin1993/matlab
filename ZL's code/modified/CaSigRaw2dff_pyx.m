function [ output_args ] = CaSigRaw2dff_pyx( CaSigFileName,BehFileName )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Calcium signal data
%     load(CaSigFileName);
%     load(BehFileName);
load('D:\xulab\imaging_data\pyx083-20180522\im_data_reg_cpu\result_save\CaTrialsSIM_pyx083_20180522_920nm_power50_zoom4x_dftReg_.mat');
load('D:\xulab\behavior\pyx083\Data_Virables\2018_05_22_pyx083-imaging_Virables.mat');
    Ca_raw = SavedCaTrials.f_raw;
    Ca_raw_ring = SavedCaTrials.RingF;
    FrameTime = SavedCaTrials.FrameTime;
    TrialN = SavedCaTrials.TrialNum;
    Frames = SavedCaTrials.nFrames; 
    nROIs = SavedCaTrials.nROIs;
    % behavior data
    Miss_Ind = Data_extract.Miss_Ind;
    Action_choice = Data_extract.Action_choice;
%     Trial_type = Data_extract.Trial_type;
%     Fre = Data_extract.Frequencies;
    Fre =Data_extract.Stimuli;
    Probe_Ind = Data_extract.Probe_index;
%     Tone_frequency = Data_extract.Tone_frequency;
    Tone_frequency =Data_extract.Stimulus;
    OnsetT= Data_extract.Tone_onset_time;
    MinOnsFrame = min(OnsetT)/FrameTime;
    AnswerT = Data_extract.Answer_time;
    RewardT = Data_extract.Reward_time;
    if isfield(Data_extract,'IsRewInd')
        IsRewInd = Data_extract.IsRewInd;
    end
    
    MinPreTime = min(OnsetT);
    MinPreFrames = round(MinPreTime/FrameTime);
    ProFrames = Frames - max(OnsetT)/FrameTime;
    TimeScale = [1-MinPreFrames:ProFrames]*FrameTime;
    
    for triN = 1:length(Tone_frequency)
        L_lickT{triN} = Data_extract.Left_lick_time{triN}/FrameTime;
        R_lickT{triN} = Data_extract.Right_lick_time{triN}/FrameTime;
    end
    bound = length(Fre)/2;
    for n = 1:length(Fre)      
       if n <= bound  %low fre-->left, high fre --> right
           TrialType(n,:) = Tone_frequency == Fre(n) & Action_choice == 0;
           TrialType(n+length(Fre),:) = Tone_frequency == Fre(n) & Action_choice == 1;
           TrialType(n+2*length(Fre),:) = Tone_frequency == Fre(n) & Miss_Ind == 1;
       else
           TrialType(n,:) = Tone_frequency == Fre(n) & Action_choice == 1;
           TrialType(n+length(Fre),:) = Tone_frequency == Fre(n) & Action_choice == 0;
           TrialType(n+2*length(Fre),:) = Tone_frequency == Fre(n) & Miss_Ind == 1;
       end
    end
    % get dff  using min Signal of each trials  as f0
    temp_CaS = cell(1,SavedCaTrials.nROIs);
    dff_single = cell(1,SavedCaTrials.nROIs);
    dff_mean = cell(1,SavedCaTrials.nROIs);
    dff_single_ring = cell(1,SavedCaTrials.nROIs);
    dff_mode = cell(1,SavedCaTrials.nROIs);
    [ModeOfROI,AccumF,PreToneCaS,MeanF0s] = Find_mode_fun_pyx(Data_extract, SavedCaTrials);
    FrameN=cellfun(@(x) size(x,2),Ca_raw);
    maxFramN=max(FrameN);
    for nROIs = 1:SavedCaTrials.nROIs
%         temp_CaS{nROIs} = squeeze(Ca_raw(:,nROIs,:));
%         temp_CaS_ring{nROIs} = squeeze(Ca_raw_ring(:,nROIs,:));
        temp_CaS{nROIs} = zeros(TrialN,Frames);
        temp_CaS_ring{nROIs} = zeros(TrialN,Frames);
        for i=1:length(Ca_raw) %for each trials            
            temp_CaS{nROIs} (i,1:FrameN(i))=Ca_raw{i}(nROIs,:);
            temp_CaS{nROIs} (i,FrameN(i)+1:end)=nan;
            temp_CaS_ring{nROIs}(i,1:FrameN(i))=Ca_raw_ring{i}(nROIs,:);
            temp_CaS_ring{nROIs} (i,FrameN(i)+1:end)=nan;
        end
        dff_single{nROIs} = zeros(TrialN,Frames);
        dff_single_ring{nROIs} = zeros(TrialN,Frames);
        dff_mode{nROIs} = zeros(TrialN,Frames);
        dff_mode{nROIs} = (temp_CaS{nROIs} - ModeOfROI(nROIs))/ModeOfROI(nROIs);
        dff_mean{nROIs} = zeros(TrialN,Frames);
        
        max_temp = 0;
        min_temp = 0;
        for T = 1:TrialN
            dff_single{nROIs}(T,:) = (temp_CaS{nROIs}(T,:) - nanmin(temp_CaS{nROIs}(T,:)))/nanmin(temp_CaS{nROIs}(T,:));
            dff_single_ring{nROIs}(T,:) = (temp_CaS_ring{nROIs}(T,:) - nanmin(temp_CaS_ring{nROIs}(T,:)))/nanmin(temp_CaS_ring{nROIs}(T,:));
            if sum(nanmax(nanmax(dff_single{nROIs}(T,:)))>60)
                dff_single{nROIs}(T,:) = (temp_CaS{nROIs}(T,:) - ModeOfROI(nROIs))/ModeOfROI(nROIs);
            end
            dff_mean{nROIs}(T,:) = (temp_CaS{nROIs}(T,:) - MeanF0s{nROIs}(T))/MeanF0s{nROIs}(T);
        end
        
        max_temp = nanmax([nanmax(dff_single{nROIs}) max_temp]);
        min_temp = nanmin([nanmin(dff_single{nROIs}) min_temp]);    
        maxCaS(nROIs) = max_temp;
        minCaS(nROIs) = min_temp;    
    end
%     if isfield(Data_extract,'IsRewInd')
%         save(['dff_' CaSigFileName(13:end-8) '.mat'],'dff_single','dff_single_ring','dff_mode','dff_mean','Miss_Ind','Action_choice',...
%             'Fre','Probe_Ind','Tone_frequency','OnsetT','MinOnsFrame','AnswerT',...
%             'RewardT','L_lickT','R_lickT','TrialType','FrameTime','Frames','Trial_type','nROIs','MinPreFrames',...
%             'ProFrames','TimeScale','IsRewInd','a1','b1','c1','Data_extract','Note','Scores','toneOct');
%     else
%         save(['dff_' CaSigFileName(13:end-8) '.mat'],'dff_single','dff_single_ring','dff_mode','dff_mean','Miss_Ind','Action_choice',...
%             'Fre','Probe_Ind','Tone_frequency','OnsetT','MinOnsFrame','AnswerT',...
%             'RewardT','L_lickT','R_lickT','TrialType','FrameTime','Frames','Trial_type','nROIs','MinPreFrames',...
%             'ProFrames','TimeScale','a1','b1','c1','Data_extract','Note','Scores','toneOct');
%     end
saveFileName=['D:\xulab\imaging_data\pyx083-20180522\im_data_reg_cpu\result_save\dff_pyx083_20180522_920nm_power50_zoom4x_dft.mat'];
save(saveFileName);
end





