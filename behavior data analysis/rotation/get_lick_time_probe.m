function [first_AnswerLick,first_lick,freq_type,licktime_c_Opto_delta,licktime_w_Opto_delta,licktime_m_Opto_delta,licktime_c_nonOpto_delta,licktime_w_nonOpto_delta,licktime_m_nonOpto_delta] = get_lick_time_probe(session_data)
% [optoLickRate, nonOptoLickRate] = get_lick_rate(session_data, side)
% Compute lickrate of trials with or without light stimulation.
% session_data is the loaded session data
% side, 'left' or 'right'
% Freq, 'High' or 'Low'
% Ac, 'c','w'or 'm'


% inds_freqHigh = cellfun(@(x) x.Trial_Type == 1, session_data);
% inds_freqLow = cellfun(@(x) x.Trial_Type == 0, session_data);
inds_correct = cellfun(@(x) x.Time_reward ~= 0, session_data);
inds_wrong = cellfun(@(x) x.Time_reward == 0, session_data);
inds_miss = cellfun(@(x) x.Action_choice == 2, session_data);

inds_OptoP = cellfun(@(x) x.Trial_isOptoTraingTrial == 1 ,session_data);
inds_OptoT = cellfun(@(x) x.Trial_isOptoProbeTrial == 1 ,session_data);
inds_nonOptoP = cellfun(@(x) x.Trial_isOptoTraingTrial == 0 ,session_data);
inds_nonOptoT = cellfun(@(x) x.Trial_isOptoProbeTrial == 0 ,session_data);

inds_Opto = inds_OptoP | inds_OptoT;
inds_nonOpto = inds_nonOptoP & inds_nonOptoT;

inds_c_o =  inds_correct & inds_Opto;
inds_w_o =  inds_wrong & inds_Opto;
inds_m_o =  inds_miss & inds_Opto;
inds_c_no = inds_correct & inds_nonOpto;
inds_w_no =  inds_wrong & inds_nonOpto;
inds_m_no = inds_miss & inds_nonOpto;

trial_isProbe_flag = cellfun(@(x) x.Trial_isProbeTrial == 1, session_data);
soundDur = session_data{1}.stimDuration;

for n=1:length(session_data);
%     if trial_isProbe_flag(n) == 1 ;
%             freq(n) = session_data{n}.Stim_Probe_pureTone_freq;
%     elseif trial_isProbe_flag(n) == 0 ;
%             freq(n) = session_data{n}.Stim_toneFreq;
%     end
    freq(n) = session_data{n}.Stim_toneFreq;
end

freq_type = unique(freq);
freq_type = sort(freq_type);
num_freq = length(freq_type);

for fq = 1:num_freq
for t = 1:length(session_data)
if freq(t) == freq_type(fq);
    inds_freq(fq,t) = true;
else
    inds_freq(fq,t) = false;
end
end
end
    
 

time_stimOnset = cellfun(@(x) x.Time_stimOnset, session_data);
%min_time_stimOnset = min(time_stimOnset);


 licktime=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
 licktime_inds = [];
 licktime_delta=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
 
 for f = 1:num_freq
     licktime_c_Opto_delta{f}=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
     licktime_c_nonOpto_delta{f}=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
     licktime_w_Opto_delta{f}=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
     licktime_w_nonOpto_delta{f}=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
     licktime_m_Opto_delta{f}=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
     licktime_m_nonOpto_delta{f}=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
 end
 
%  licktime_hc_Opto_delta=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
%  licktime_hw_Opto_delta=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
%  licktime_hm_Opto_delta=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
%  licktime_hc_nonOpto_delta=struct('Action  lickTimeRight',[],'ActionlickTimeLeft',[]);
%  licktime_hw_nonOpto_delta=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
%  licktime_hm_nonOpto_delta=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
%  licktime_lc_Opto_delta=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
%  licktime_lw_Opto_delta=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
%  licktime_lm_Opto_delta=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
%  licktime_lc_nonOpto_delta=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
%  licktime_lw_nonOpto_delta=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
%  licktime_lm_nonOpto_delta=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
 
 for i=1:length(session_data)
    delta = time_stimOnset(i); 
    LickRight_time_delta= [];
    LickLeft_time_delta= [];
    LickRight_time_answer = [];
    LickLeft_time_answer = [];
    
    LickRight_time_str= strsplit(session_data{i}.Action_lickTimeRight,'|');
    LickRight_time_str(:,1)=[];
    LickRight_time=str2double( LickRight_time_str );
    LickLeft_time_str = strsplit(session_data{i}.Action_lickTimeLeft,'|');
    LickLeft_time_str(:,1)=[];
    LickLeft_time=str2double( LickLeft_time_str );
    if isnan(LickLeft_time)
        LickLeft_time=[];
    end
    if isnan(LickRight_time)
        LickRight_time=[];
    end
    for j = 1:length(LickRight_time)
        LickRight_time_delta(j) = LickRight_time(j) - delta;
        LickRight_time_answer(j) = LickRight_time(j) - delta - soundDur;
    end
    for j = 1:length(LickLeft_time)
        LickLeft_time_delta(j) = LickLeft_time(j) - delta;
        LickLeft_time_answer(j) = LickLeft_time(j) - delta - soundDur;
    end
    licktime(i).ActionlickTimeRight=LickRight_time;
    licktime(i).ActionlickTimeLeft=LickLeft_time;
    
    licktime_delta(i).ActionlickTimeRight=LickRight_time_delta;
    licktime_delta(i).ActionlickTimeLeft=LickLeft_time_delta;  
    if length([LickRight_time_delta,LickLeft_time_delta]) ~= 0
        first_lick(i) = nanmin([LickRight_time_delta,LickLeft_time_delta]);
        first_AnswerLick(i) = nanmin([LickRight_time_answer,LickLeft_time_answer]);
    elseif length([LickRight_time_delta,LickLeft_time_delta]) == 0
        first_lick(i) = NaN;
        first_AnswerLick(i) = NaN;   
    end
 end
%   for f2 = 1:num_freq
%      licktime_c_Opto_delta{f2}=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
%      licktime_c_nonOpto_delta{f2}=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
%      licktime_w_Opto_delta{f2}=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
%      licktime_w_nonOpto_delta{f2}=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
%      licktime_m_Opto_delta{f2}=struct('ActionlickTimeRight',[],'ActionlickTimeLeft',[]);
%      licktime_m_nonOpto_delta{f2}=struct('ActionlickTimeRig  ht',[],'ActionlickTimeLeft',[]);
%   end
 for f3 = 1:num_freq
     licktime_c_Opto_delta{f3} = licktime_delta(inds_c_o & inds_freq(f3,:));
     licktime_w_Opto_delta{f3} = licktime_delta(inds_w_o & inds_freq(f3,:));
     licktime_m_Opto_delta{f3} = licktime_delta(inds_m_o & inds_freq(f3,:));
     licktime_c_nonOpto_delta{f3} = licktime_delta(inds_c_no & inds_freq(f3,:));
     licktime_w_nonOpto_delta{f3} = licktime_delta(inds_w_no & inds_freq(f3,:));
     licktime_m_nonOpto_delta{f3} = licktime_delta(inds_m_no & inds_freq(f3,:));
 end
