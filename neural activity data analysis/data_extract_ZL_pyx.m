%%//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
%%%/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
%%%     This function is available for          
%%%       //ZL master4 opto// to extract data

%  Cover the original files?  Yes: 1      No: 0

Cover_option = 1;
% CurrFolder = pwd;
% CurrFolder = 'F:\2P\pyx252_20200111\im_data_reg\result_save';
CurrFolder = 'D:\pyx354_20210423\im_data_reg\result_save';% rootpath='F:\FP\pyx241_20191130';
% savePath = [CurrFolder '\Data_Virables\'];
% if ~exist('Data_Virables')
%     mkdir('Data_Virables');
% end
savePath =[CurrFolder '\'];
files = dir(strcat(CurrFolder,'\*imaging.mat'));
for i = 1:length(files)
    clearvars -except 'files' 'i' 'savePath' 'Cover_option' 'CurrFolder'
    fn = files(i).name;
    load(strcat(CurrFolder,'\',fn));
    if isfield('animalName',SessionSettings{1})
        Animal_Name_Setting = SessionSettings{1}.animalName;
    else
        Animal_Name_Setting = [];
    end
    Animal_Name_File = fn; 
    Rig_name = SessionSettings{2}.Rig_name;
    L_R_WaterValveDuration = [num2str(SessionSettings{1}.leftWaterValveDuration) '/' num2str(SessionSettings{1}.rightWaterValveDuration)];
    TriType = cellfun(@(x) x.Trial_Type, SessionResults);
    TriChoice = cellfun(@(x) x.Action_choice, SessionResults);
    Trial_index = cellfun(@(x) x.Trial_inds, SessionResults);
    Trial_number = cellfun(@(x) x.Trial_Num, SessionResults);
    Answer_time = cellfun(@(x) x.Time_answer, SessionResults);
    Reward_time = cellfun(@(x) x.Time_reward, SessionResults);
    OptoStimOnsetTime = cellfun(@(x) x.Time_optoStimOnset, SessionResults);
    OptoStimOffTime = cellfun(@(x) x.Time_optoStimOffTime, SessionResults);

    
    Answer_Period = SessionSettings{1}.answerPeriod;
    Experiment_date = fn(8:15);                 %==========
    Miss_Ind = cellfun(@(x) x.Action_choice == 2, SessionResults);
    Vio_Ind=cellfun(@(x)  x.Action_choice == 3, SessionResults);
    Probe_index = cellfun(@(x) x.Trial_isProbeTrial == 1, SessionResults);
    training_opto = cellfun(@(x) x.Trial_isOptoTraingTrial == 1, SessionResults);
    probe_opto = cellfun(@(x) x.Trial_isOptoProbeTrial == 1, SessionResults);
    Opto_trial_index = zeros(1,length(training_opto));
    Opto_trial_index(training_opto | probe_opto) = 1;
    Tone_onset_time = cellfun(@(x) x.Time_stimOnset, SessionResults);
    L_lick_time = cellfun(@(x) x.Action_lickTimeLeft, SessionResults,'UniformOutput',false);
    R_lick_time = cellfun(@(x) x.Action_lickTimeRight, SessionResults,'UniformOutput',false);
    Time_trialStart = cellfun(@(x) x.Time_trialStart, SessionResults,'UniformOutput',false);
    
    MinOnsetTime = min(cellfun(@(x) x.Time_stimOnset, SessionResults));
    
    for j = 1:length(SessionResults)
        if isfield(SessionResults{j},'stimDuration') == 1
            if ~isempty(SessionResults{j}.stimDuration)
                ToneDurTime(j) = SessionResults{j}.stimDuration;
            else
                ToneDurTime(j)=NaN;
            end   
        end
        % Lick Time
        temp_1 = regexp(cell2mat(L_lick_time(j)), '\|', 'split');
        temp_2 = regexp(cell2mat(R_lick_time(j)), '\|', 'split');
        temp_3 = [];
        temp_4 = [];
        for k = 1:length(temp_1)
            if isempty(temp_1(k))
                    temp_3 = temp_3;
            else
                temp_3 = [temp_3 str2num(temp_1{k})];
            end
        end
        for k = 1:length(temp_2)
            if isempty(temp_2(k))
                    temp_4 = temp_4;
            else
                temp_4 = [temp_4 str2num(temp_2{k})];
            end
        end
        Left_lick_time{j} = temp_3;
        Right_lick_time{j} = temp_4;
%         % Tone Frequences
%         if isfield(SessionResults{j},'Stim_Probe_pureTone_freq') == 1
%            Tone_frequency(j) = SessionResults{j}.Stim_Probe_pureTone_freq;
%         else
%            Tone_frequency(j) = SessionResults{j}.Stim_toneFreq;
%         end
        %stimulus
        Stimulus(j)=SessionResults{j}.Stim_clickRate;
    end
%     Max_frequency = double(max(Tone_frequency));
%     Min_frequency = double(min(Tone_frequency));
%     DetBundary_frequency = 2^(log2(Min_frequency) + (log2(Max_frequency)-log2(Min_frequency))/2);
    
    Data_extract.Animal_Name_Setting = Animal_Name_Setting;
    Data_extract.Animal_Name_File = Animal_Name_File;
    
    Data_extract.Experiment_date = Experiment_date;
    Data_extract.Rig_name = Rig_name;
    Data_extract.L_R_WaterValveDuration = L_R_WaterValveDuration;
    Data_extract.Trial_index = Trial_index;
    Data_extract.Trial_number = Trial_number;
%     Data_extract.Trial_type = TriType;
    Data_extract.Answer_Period = Answer_Period;
    Data_extract.Miss_Ind = Miss_Ind;
    Data_extract.Vio_Ind = Vio_Ind;
    Data_extract.Action_choice = TriChoice;
    Data_extract.Opto_trial_index = Opto_trial_index;
    Data_extract.Probe_index = Probe_index;
%     Data_extract.Tone_frequency = Tone_frequency;
% if Stimulus>1000
%     Data_extract.Stimulus =floor(Stimulus/100);
% else
    Data_extract.Stimulus = Stimulus;%note that stim may be wrong(conbined with volume)
% end
%     Data_extract.Frequencies = unique(Tone_frequency);
    Data_extract.Stimuli = unique(Stimulus);
    Data_extract.Stim_onset_time = Tone_onset_time;
    Data_extract.Stim_offset_time = cellfun(@(x) x.Time_stimOffset, SessionResults);
    Data_extract.Trial_Type = cellfun(@(x) x.Trial_Type, SessionResults);
    Data_extract.Delay_duration = cellfun(@(x) x.Delay_duration, SessionResults);
    Data_extract.Time_delayOffset = cellfun(@(x) x.Time_delayOffset, SessionResults);
    Data_extract.Min_Onset_Time = MinOnsetTime;
    Data_extract.Tone_Duration_Time = ToneDurTime;
    Data_extract.Time_trialStart =Time_trialStart;
    Data_extract.Answer_time = Answer_time;
    Data_extract.Reward_time = Reward_time;
    Data_extract.Opto_Onset_Time = OptoStimOnsetTime;
    Data_extract.Opto_Off_Time = OptoStimOffTime;
    Data_extract.Left_lick_time = Left_lick_time;
    Data_extract.Right_lick_time = Right_lick_time;
    Data_extract.Go_time=cellfun(@(x) x.Time_delayOnset+x.Delay_duration, SessionResults);
    if isfield(SessionSettings{1, 1},'clickRate_L')%mat file that translated recently
        if SessionSettings{1, 1}.clickRate_L>SessionSettings{1, 1}.clickRate_R
            Data_extract.rule='low click rate-right';
        else
            Data_extract.rule='low click rate-left';
        end
    else
        if SessionSettings{1, 1}.clickRate_<50
            Data_extract.rule='low click rate-right';
        else
            Data_extract.rule='low click rate-left';
        end
    end
%     Data_extract.Determine_boundary = DetBundary_frequency;
    
    if (exist([savePath fn(1:end-4) '_Virables'],'file') == 0)
        save([savePath fn(1:end-4) '_Virables'],'Data_extract');
    elseif (exist([savePath fn(1:end-4) '_Virables'],'file') == 1)
        if Cover_option == 0
            continue;
        elseif Cover_option == 1
            save([savePath fn(1:end-4) '_Virables'],'Data_extract');
        end
    end    
end
