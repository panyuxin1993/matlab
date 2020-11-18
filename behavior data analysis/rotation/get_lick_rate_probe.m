function [freq_type,nonOptoLickRate_c_r,nonOptoLickRate_c_l,nonOptoLickRate_w_r,nonOptoLickRate_w_l,OptoLickRate_c_r,OptoLickRate_c_l,OptoLickRate_w_r,OptoLickRate_w_l] = get_lick_rate_probe(session_data,bin,Tlength,sample)

% bin = 100;
% length = 6000;
% sample = Tlength/bin +1;

inds_freqHigh = cellfun(@(x) x.Trial_Type == 1, session_data);
inds_freqLow = cellfun(@(x) x.Trial_Type == 0, session_data);
inds_correct = cellfun(@(x) x.Time_reward ~= 0, session_data);
inds_miss = cellfun(@(x) x.Action_choice == 2, session_data);
inds_wrong  = ~inds_correct & ~inds_miss;

inds_OptoP = cellfun(@(x) x.Trial_isOptoTraingTrial == 1 ,session_data);
inds_OptoT = cellfun(@(x) x.Trial_isOptoProbeTrial == 1 ,session_data);
inds_nonOptoP = cellfun(@(x) x.Trial_isOptoTraingTrial == 0 ,session_data);
inds_nonOptoT = cellfun(@(x) x.Trial_isOptoProbeTrial == 0 ,session_data);

inds_Opto = inds_OptoP | inds_OptoT;
inds_nonOpto = inds_nonOptoP & inds_nonOptoT;

inds_c_o = inds_correct & inds_Opto;
inds_w_o =  inds_wrong & inds_Opto;
inds_m_o =  inds_miss & inds_Opto;
inds_c_no =  inds_correct & inds_nonOpto;
inds_w_no =  inds_wrong & inds_nonOpto;
inds_m_no =  inds_miss & inds_nonOpto;


trial_isProbe_flag = cellfun(@(x) x.Trial_isProbeTrial == 1, session_data);
% inds_Opto = cellfun(@(x) x.Trial_isProbeTrial == 1, session_data);
% inds_nonOpto = cellfun(@(x) x.Trial_isProbeTrial == 0, session_data);

for n=1:length(session_data);
    if trial_isProbe_flag(n) == 1 ;
            freq(n) = session_data{n}.Stim_Probe_pureTone_freq;
    elseif trial_isProbe_flag(n) == 0 ;
            freq(n) = session_data{n}.Stim_toneFreq;
    end
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
 
 
 for i=1:length(session_data)
    delta = time_stimOnset(i);
    LickRight_time_delta= [];
    LickLeft_time_delta= [];
    
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
    end
    for j = 1:length(LickLeft_time)
        LickLeft_time_delta(j) = LickLeft_time(j) - delta;
    end
    licktime(i).ActionlickTimeRight=LickRight_time;
    licktime(i).ActionlickTimeLeft=LickLeft_time;
    
    licktime_delta(i).ActionlickTimeRight=LickRight_time_delta;
    licktime_delta(i).ActionlickTimeLeft=LickLeft_time_delta;   
 end
 
 for f3 = 1:num_freq
     licktime_c_Opto_delta{f3} = licktime_delta(inds_c_o & inds_freq(f3,:));
     licktime_w_Opto_delta{f3} = licktime_delta(inds_w_o & inds_freq(f3,:));
     licktime_m_Opto_delta{f3} = licktime_delta(inds_m_o & inds_freq(f3,:));
     licktime_c_nonOpto_delta{f3} = licktime_delta(inds_c_no & inds_freq(f3,:));
     licktime_w_nonOpto_delta{f3} = licktime_delta(inds_w_no & inds_freq(f3,:));
     licktime_m_nonOpto_delta{f3} = licktime_delta(inds_m_no & inds_freq(f3,:));
     lickNum_c_nonOpto_r{f3} = zeros(length(licktime_c_nonOpto_delta{f3}),sample);
     lickNum_w_nonOpto_r{f3} = zeros(length(licktime_w_nonOpto_delta{f3}),sample);
     lickNum_m_nonOpto_r{f3} = zeros(length(licktime_m_nonOpto_delta{f3}),sample);
     lickNum_c_Opto_r{f3} = zeros(length(licktime_c_Opto_delta{f3}),sample);
     lickNum_w_Opto_r{f3} = zeros(length(licktime_w_Opto_delta{f3}),sample);
     lickNum_m_Opto_r{f3} = zeros(length(licktime_m_Opto_delta{f3}),sample);
     lickNum_c_nonOpto_l{f3} = zeros(length(licktime_c_nonOpto_delta{f3}),sample);
     lickNum_w_nonOpto_l{f3} = zeros(length(licktime_w_nonOpto_delta{f3}),sample);
     lickNum_m_nonOpto_l{f3} = zeros(length(licktime_m_nonOpto_delta{f3}),sample);
     lickNum_c_Opto_l{f3} = zeros(length(licktime_c_Opto_delta{f3}),sample);
     lickNum_w_Opto_l{f3} = zeros(length(licktime_w_Opto_delta{f3}),sample);
     lickNum_m_Opto_l{f3} = zeros(length(licktime_m_Opto_delta{f3}),sample);
 end
 
 for f4 = 1:num_freq
 for a = 0:bin:Tlength
     
     for tn = 1:length(licktime_c_nonOpto_delta{f4})
     lickNum_c_nonOpto_r{f4}(tn,(a/bin+1)) = length(find(licktime_c_nonOpto_delta{f4}(tn).ActionlickTimeRight>(a-1000) & licktime_c_nonOpto_delta{f4}(tn).ActionlickTimeRight<(a-1000+bin)));
     lickNum_c_nonOpto_l{f4}(tn,(a/bin+1)) = length(find(licktime_c_nonOpto_delta{f4}(tn).ActionlickTimeLeft>(a-1000) & licktime_c_nonOpto_delta{f4}(tn).ActionlickTimeLeft<(a-1000+bin)));
     end
     for t1 = 1:length(licktime_w_nonOpto_delta{f4})
     lickNum_w_nonOpto_r{f4}(t1,(a/bin+1)) = length(find(licktime_w_nonOpto_delta{f4}(t1).ActionlickTimeRight>(a-1000) & licktime_w_nonOpto_delta{f4}(t1).ActionlickTimeRight<(a-1000+bin)));
     lickNum_w_nonOpto_l{f4}(t1,(a/bin+1)) = length(find(licktime_w_nonOpto_delta{f4}(t1).ActionlickTimeLeft>(a-1000) & licktime_w_nonOpto_delta{f4}(t1).ActionlickTimeLeft<(a-1000+bin)));
     end
     for t2 = 1:length(licktime_c_Opto_delta{f4})
     lickNum_c_Opto_r{f4}(t2,(a/bin+1)) = length(find(licktime_c_Opto_delta{f4}(t2).ActionlickTimeRight>(a-1000) & licktime_c_Opto_delta{f4}(t2).ActionlickTimeRight<(a-1000+bin)));
     lickNum_c_Opto_l{f4}(t2,(a/bin+1)) = length(find(licktime_c_Opto_delta{f4}(t2).ActionlickTimeLeft>(a-1000) & licktime_c_Opto_delta{f4}(t2).ActionlickTimeLeft<(a-1000+bin)));
     end
     for t3 = 1:length(licktime_w_Opto_delta{f4})
     lickNum_w_Opto_r{f4}(t3,(a/bin+1)) = length(find(licktime_w_Opto_delta{f4}(t3).ActionlickTimeRight>(a-1000) & licktime_w_Opto_delta{f4}(t3).ActionlickTimeRight<(a-1000+bin)));
     lickNum_w_Opto_l{f4}(t3,(a/bin+1)) = length(find(licktime_w_Opto_delta{f4}(t3).ActionlickTimeLeft>(a-1000) & licktime_w_Opto_delta{f4}(t3).ActionlickTimeLeft<(a-1000+bin)));
     end
     nonOptoLickRate_c_r{f4}(a/bin+1) = mean(lickNum_c_nonOpto_r{f4}(:,(a/bin+1)));
     nonOptoLickRate_c_l{f4}(a/bin+1) = mean(lickNum_c_nonOpto_l{f4}(:,(a/bin+1)));
     nonOptoLickRate_w_r{f4}(a/bin+1) = mean(lickNum_w_nonOpto_r{f4}(:,(a/bin+1)));
     nonOptoLickRate_w_l{f4}(a/bin+1) = mean(lickNum_w_nonOpto_l{f4}(:,(a/bin+1)));
     OptoLickRate_c_r{f4}(a/bin+1) = mean(lickNum_c_Opto_r{f4}(:,(a/bin+1)));
     OptoLickRate_c_l{f4}(a/bin+1) = mean(lickNum_c_Opto_l{f4}(:,(a/bin+1)));
     OptoLickRate_w_r{f4}(a/bin+1) = mean(lickNum_w_Opto_r{f4}(:,(a/bin+1)));
     OptoLickRate_w_l{f4}(a/bin+1) = mean(lickNum_w_Opto_l{f4}(:,(a/bin+1)));
 end
 end




