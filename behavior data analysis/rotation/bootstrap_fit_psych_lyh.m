% %%
% k = 0; % num of sessions
% cd('G:\behavior_data\O0302\raw data')
% fn_str_1 = {'O0302_20151210', 'O0302_20151213','O0302_20151211'};
% for a = 1:length(fn_str_1)
%     k = k+1;
%     load(fn_str_1{a});
%     [fractions_withOpto{a},fractions_withoutOpto{a},freq{a}] = get_RP(SessionResults);
% end
% fractions_withOpto_av = mean(fractions_withOpto);
% fractions_withoutOpto_av = mean(fractions_withoutOpto);
% % freq = freq{1};
% errorbar(toneOct,fractions_withoutOpto_av,std(fractions_withoutOpto_av/length(fractions_withoutOpto_av)) ,'color','k','linewidth',2);
% errorbar(toneOct,fractions_withOpto_av,std(fractions_withOpto_av/length(fractions_withOpto_av)) ,'color','b','linewidth',2);

%%
file = dir('E:\我的文件\上海神经所\轮转\徐宁龙\训鼠分析\curve\*.mat');
for f = 1:length(file)
filename = file(f).name;
load(['E:\我的文件\上海神经所\轮转\徐宁龙\训鼠分析\curve\' filename])

trial_isProbe_flag = cellfun(@(x) x.Trial_isProbeTrial == 1, SessionResults);
trial_isOptoProbe_flag = cellfun(@(x) x.Trial_isOptoProbeTrial == 1, SessionResults);
trial_isOptoTraining_flag = cellfun(@(x) x.Trial_isOptoTraingTrial == 1, SessionResults);
freq = [];
for n=1:length(SessionResults);
%     if trial_isProbe_flag(n) == 1 ;
%         if ismember(SessionResults{n}.Stim_Probe_pureTone_freq,freq) == 0;
%             freq = [freq;SessionResults{n}.Stim_Probe_pureTone_freq];
%   
%         end
%     elseif trial_isProbe_flag(n) == 0 ;
%         if ismember(SessionResults{n}.Stim_toneFreq,freq)== 0;
%             freq = [freq;SessionResults{n}.Stim_toneFreq];
%         end
%     end
        freq = [freq;SessionResults{n}.Stim_toneFreq];
end


freq = sort(freq);
freq = double(freq);

% Part2 -------------------------
% Generates matrix "right_prob" with frequences,
% rightward choice fraction and number of trials for each of those frequencies.

fractions_withOpto = [];
fractions_withoutOpto = [];
count_withoutOpto = [];
count_withOpto = [];
for m = 1:length(freq);
    r_withoutOpto = 0;
    r_withOpto = 0;
    action_choice_withoutOpto{m} = [];
    action_choice_withOpto{m} = [];
    for n = 1:length(SessionResults);
        if trial_isProbe_flag(n) == 1 && trial_isOptoProbe_flag(n) == 0;
            %if SessionResults{n}.Stim_Probe_pureTone_freq == freq(m);
            if SessionResults{n}.Stim_toneFreq == freq(m);
                action_choice_withoutOpto{m} = [action_choice_withoutOpto{m};SessionResults{n}.Action_choice];
            end
        elseif trial_isProbe_flag(n) == 0  && trial_isOptoTraining_flag(n) == 0;
            if SessionResults{n}.Stim_toneFreq == freq(m);
                action_choice_withoutOpto{m} = [action_choice_withoutOpto{m};SessionResults{n}.Action_choice];
            end
        end
    end
     for n = 1:length(SessionResults);
        if trial_isOptoProbe_flag(n) == 1 && trial_isOptoProbe_flag(n) == 1;
%             if SessionResults{n}.Stim_Probe_pureTone_freq == freq(m);
              if SessionResults{n}.Stim_toneFreq == freq(m);
                action_choice_withOpto{m} = [action_choice_withOpto{m};SessionResults{n}.Action_choice];
            end
        elseif trial_isProbe_flag(n) == 0  && trial_isOptoTraining_flag(n) == 1;
            if SessionResults{n}.Stim_toneFreq == freq(m);
                action_choice_withOpto{m} = [action_choice_withOpto{m};SessionResults{n}.Action_choice];
            end
        end
    end
    
     % Optional, remove all MISS trials.
    action_choice_withoutOpto{m}(action_choice_withoutOpto{m} == 2) = [];
    
    for n = 1:length(action_choice_withoutOpto{m});
        if action_choice_withoutOpto{m}(n) == 1;
            r_withoutOpto = r_withoutOpto + 1;
        end
    end
    count_withoutOpto = [count_withoutOpto;length(action_choice_withoutOpto{m})];
    
    rightward_withoutOpto = r_withoutOpto / length(action_choice_withoutOpto{m});
    fractions_withoutOpto = [fractions_withoutOpto;rightward_withoutOpto];
    
    action_choice_withOpto{m}(action_choice_withOpto{m} == 2) = [];
    
    for n = 1:length(action_choice_withOpto{m});
        if action_choice_withOpto{m}(n) == 1;
            r_withOpto = r_withOpto + 1;
        end
    end
    count_withOpto = [count_withOpto;length(action_choice_withOpto{m})];
    
    rightward_withOpto = r_withOpto / length(action_choice_withOpto{m});
    fractions_withOpto = [fractions_withOpto;rightward_withOpto];
    
    [phat_no{m},ci_no{m}] = binofit(length(find(action_choice_withoutOpto{m} == 1)),length(action_choice_withoutOpto{m}));
    [phat_o{m},ci_o{m}] = binofit(length(find(action_choice_withOpto{m} == 1)),length(action_choice_withOpto{m}));
   
end


tic
toneOct = log2(freq/7000);
frac_choice_right_Contr_Opto = fractions_withOpto;
frac_choice_right_Contr_nonOpto = fractions_withoutOpto;
%parpool('local', 4);
%options = statset('UseParallel',true);
%[bootstat_control, bootsam] = bootstrp(500, @(x, y) fit_logistic_psych_con(x, y,0), toneOct, frac_choice_right_Contr,'Options',options);
[bootstat_control_Opto, bootsam_Opto] = bootstrp(100, @(x, y) fit_logistic_psych_con(x, y,0), toneOct, frac_choice_right_Contr_Opto);
[bootstat_control_nonOpto, bootsam_nonOpto] = bootstrp(100, @(x, y) fit_logistic_psych_con(x, y,0), toneOct, frac_choice_right_Contr_nonOpto);
toc

%
gof_nonOpto = [bootstat_control_nonOpto.gof];
rmse_nonOpto = [gof_nonOpto.rmse];
sse_nonOpto = [gof_nonOpto.sse];
rsq_nonOpto = [gof_nonOpto.rsquare];
adjrsquare_nonOpto = [gof_nonOpto.adjrsquare];

inds_use1_nonOpto = find(rmse_nonOpto < prctile(rmse_nonOpto,50));
inds_use2_nonOpto = find(rsq_nonOpto > prctile(rsq_nonOpto,50));
inds_use3_nonOpto = find(mode(rmse_nonOpto));

a1_nonOpto = median([bootstat_control_nonOpto(inds_use1_nonOpto).a]);
b1_nonOpto = median([bootstat_control_nonOpto(inds_use1_nonOpto).b]);
c1_nonOpto = median([bootstat_control_nonOpto(inds_use1_nonOpto).c]);

gof_Opto = [bootstat_control_Opto.gof];
rmse_Opto = [gof_Opto.rmse];
sse_Opto = [gof_Opto.sse];
rsq_Opto = [gof_Opto.rsquare];
adjrsquare_Opto = [gof_Opto.adjrsquare];

inds_use1_Opto = find(rmse_Opto < prctile(rmse_Opto,50));
inds_use2_Opto = find(rsq_Opto > prctile(rsq_Opto,50));
inds_use3_Opto = find(mode(rmse_Opto));

a1_Opto = median([bootstat_control_Opto(inds_use1_Opto).a]);
b1_Opto = median([bootstat_control_Opto(inds_use1_Opto).b]);
c1_Opto = median([bootstat_control_Opto(inds_use1_Opto).c]);

%
% pct = 40;
% a1 = prctile([bootstat_control.a], pct);
% b1 = prctile([bootstat_control.b], pct);
% c1 = prctile([bootstat_control.c], pct);
%
%a1 = mode([bootstat_control.a]);
% b1 = mode([bootstat_control.b]);
% c1 = mode([bootstat_control.c]);

%
x1 = linspace(min(toneOct)-0.1, max(toneOct)+0.1, 100);

% c1 = prctile([bootstat_control.c],35);
figure;
clf; hold on;
h_data_withoutOpto = plot(toneOct, fractions_withoutOpto, 'k.','markersize',15); 
h_data_withOpto = plot(toneOct, fractions_withOpto, 'b.','markersize',15); 
y1_nonOpto = a1_nonOpto./(1+exp(-(x1 - b1_nonOpto)/c1_nonOpto));
y1_Opto = a1_Opto./(1+exp(-(x1 - b1_Opto)/c1_Opto));
h_curve_nonOpto = plot(x1,y1_nonOpto,'k','linewidth',2);
h_curve_Opto = plot(x1,y1_Opto,'b','linewidth',2);
xtick = get(gca,'XTick');
ylim([0,1.1]);
%set(gca,'XTickLabel', round(2.^xtick*10)/10);
xlabel('Octave','fontWeight','bold');
ylabel('Frac R-Choice','fontWeight','bold');
% for m = 1:length(freq)
%     ci_nonOpto_1(m) = ci_no{m}(1);
%     ci_nonOpto_2(m) = ci_no{m}(2);
%     ci_Opto_1(m) = ci_o{m}(1);
%     ci_Opto_2(m) = ci_o{m}(2);
% end
% 
%     plot(toneOct,ci_nonOpto_1,'.k');
%     hold on
%     plot(toneOct,ci_nonOpto_2,'.k');

% ci_no = cell2mat(ci_no);
% ci_o = cell2mat(ci_o);
% dy1 = ci_no(1:8);
%             dy1 = dy1';
%             x = toneOct;%[0:(length(ave_sp_dff_stim{rn,tyn,fn})-1)];
%             x = x';
%             fill([x;flipud(x)],[ci_nonOpto_1';flipud(ci_nonOpto_2')],[.5 .5 .5],'linestyle','none','facealpha',0.5);
%             hold on
%             fill([x;flipud(x)],[ci_Opto_1';flipud(ci_Opto_2')],'b','linestyle','none','facealpha',0.5);
%             hold on



filename_2 = strrep(filename,'.mat','_curve');

title(filename_2,'fontWeight','bold','interprete','none');

saveas(gcf,['D:\behavior_data\temp\performance\' filename_2 '.png']);
saveas(gcf,['D:\behavior_data\temp\performance\' filename_2 '.fig']);

close all;
clear all;


file = dir('D:\behavior_data\temp\*.mat');
end


%%
fit_param.slope = 1/c1;
fit_param.bias = b1;
fit_param.lapse = a1;
fit_param.gof_bootstrap = gof;
% %%
% [bootstat_opto, bootsam] = bootstrp(1000, @(x, y) fit_logistic_psych_opto(x, y,0), toneOct, frac_choice_right_Opto);
% %%
% gof = [bootstat_control.gof];
% rmse = [gof.rmse];
% inds_use = find(rmse < prctile(rmse,50));
% 
% a2 = median([bootstat_control(inds_use).a]);
% b2 = median([bootstat_control(inds_use).b]);
% c2 = median([bootstat_control(inds_use).c]);
% 
% 
% x2 = linspace(min(toneOct)-0.1, max(toneOct)+0.1, 100);
% 
% % c1 = prctile([bootstat_control.c],35);
% 
% y2 = a1./(1+exp(-(x2 - b2)/c2));
% h_curve_2 = plot(x2, y2,'g','linewidth',2);