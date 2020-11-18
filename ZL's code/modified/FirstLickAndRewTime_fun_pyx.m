function [ Output ] = FirstLickAndRewTime_fun_pyx( behData ,useAllTrial)
%   =========== finish revision at 20171208 ==============
%   Detailed explanation goes here
%   behData: Data_extract
%     load(fn);
% useAllTrial: 1  use all trials
% useAllTrial: TestTrialNum
    if  length(useAllTrial)==1
        selected_trials = [1:length(behData.Tone_onset_time)];
    else
        selected_trials = useAllTrial;
    end

    LLickTime = behData.Left_lick_time(selected_trials);
    RLickTime = behData.Right_lick_time(selected_trials);
    rewardTime = behData.Reward_time(selected_trials);
    OnsetT = behData.Tone_onset_time(selected_trials);
    Act = behData.Action_choice(selected_trials);
    frequency = behData.Stimulus(selected_trials);
    Mis = behData.Miss_Ind(selected_trials);
    Freq = unique(frequency);
    if isfield(behData,'Opto_trial_index')
        Opto_Trainingtrial_index = behData.Opto_trial_index(selected_trials);
    end
    if isfield(behData,'Opto_Trainingtrial_index')
        Opto_Trainingtrial_index = behData.Opto_Trainingtrial_index(selected_trials);
    end
    
    tryType = frequency>Freq(length(Freq)/2);
    if exist('isReverse','var')
        if isReverse
            tryType = ~tryType;
        end
    end
    OutcomeType = zeros(1,length(tryType));
    OutcomeType(tryType==Act) = 1;
    OutcomeType(Mis) = 2;
        
    FirstLickTime = zeros(3,length(Act));  % 1st row is original choice time(include pre-stim),  2nd is is original choice time(after stim),  3rd is time related to stim
    RewardT = zeros(3,length(Act)); % 1st row is original reward time,  2nd is time related to stim,  3rd is time related to first lick(after stim)
    Note = '1st row is original choice time(include pre-stim),  2nd is is original choice time(after stim),  3rd is time related to stim';
    for i = 1:length(Act)
       if  Act(i) < 2
           if Act(i)==1
              temp = RLickTime{i};
           elseif Act(i)==0
              temp = LLickTime{i}; 
           end
           FirstLickTime(1,i) = temp(1);
           temp(temp<OnsetT(i)) = [];
           FirstLickTime(2,i) = temp(1);
           FirstLickTime(3,i) = temp(1)-OnsetT(i);
           RewardT(1,i) = rewardTime(i);
           RewardT(2,i) = rewardTime(i) - OnsetT(i);
           RewardT(3,i) = RewardT(1,i) - FirstLickTime(2,i);  
       elseif Act(i) == 2
           FirstLickTime(:,i) = NaN;
           RewardT(:,i) = NaN;      
       end
    end
    Output.onsetTime=OnsetT;
    Output.FirstLickTime=FirstLickTime;
    Output.RewardT=RewardT;
    Output.OutcomeType=OutcomeType;
    Output.tryType=tryType;
    Output.frequency=frequency;
    Output.ActChoice=Act;
    Output.Miss=Mis;
    Output.Freq=Freq;
    Output.Opto_Trainingtrial_index=Opto_Trainingtrial_index;
    Output.Note=Output;
    if exist('isReverse','var')
        if isReverse
           Output.isReverse=isReverse;
        end
    end
end
