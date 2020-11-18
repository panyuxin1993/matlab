%%
file = dir('E:\我的文件\上海神经所\轮转\徐宁龙\训鼠分析\curve\*.mat');
for f = 1:length(file)
filename = file(f).name;
load(['E:\我的文件\上海神经所\轮转\徐宁龙\训鼠分析\curve\' filename])

[first_lick,freq_type,licktime_c_Opto_delta,licktime_w_Opto_delta,licktime_m_Opto_delta,licktime_c_nonOpto_delta,licktime_w_nonOpto_delta,licktime_m_nonOpto_delta] = get_lick_time_probe(SessionResults);
%   opto
% h= figure(2);
% set(h,'name','Lick Time(OptoTrials)');

% lick_Opto = {licktime_c_Opto_delta{1},licktime_w_Opto_delta{1},licktime_c_Opto_delta{2},licktime_w_Opto_delta{2},licktime_c_Opto_delta{3},licktime_w_Opto_delta{3},licktime_c_Opto_delta{4},licktime_w_Opto_delta{4},licktime_c_Opto_delta{5},licktime_w_Opto_delta{5},licktime_c_Opto_delta{6},licktime_w_Opto_delta{6},licktime_c_Opto_delta{7},licktime_w_Opto_delta{7},licktime_c_Opto_delta{8},licktime_w_Opto_delta{8}};
% fig_name_Opto = {'7k opto correct','7k opto wrong','10k opto correct','10k opto wrong','11k opto correct','11k opto wrong','12k opto correct','12k opto wrong','16k opto correct','16k opto wrong','18k opto correct','18k opto wrong','20k opto correct','20k opto wrong','28k opto correct','28k opto wrong'};


for g = 1:length(freq_type)
lick_Opto{1,(2*g-1)} = {licktime_c_Opto_delta{g}};
lick_Opto{1,(2*g)} = {licktime_w_Opto_delta{g}};
fig_name_Opto{(2*g-1)} = {sprintf('%d opto correct',freq_type(g))};
fig_name_Opto{(2*g)} = {sprintf('%d opto wrong',freq_type(g))};
end

%len = {lenRtime_hc_Opto,lenLtime_hc_Opto,lenRtime_hw_Opto,lenLtime_hw_Opto,lenRtime_lc_Opto,lenLtime_lc_Opto,lenRtime_lw_OptolenLtime_lw_Opto,};
for c = 1:2*(length(freq_type))
subplot(length(freq_type),2,c)
%title('28000 correct trials')
for n=1:length(lick_Opto{1,c}{1,1})
  
    for j = 1:length(lick_Opto{1,c}{1,1}(n).ActionlickTimeRight)
        line([(lick_Opto{1,c}{1,1}(n).ActionlickTimeRight(j)-0.5),(lick_Opto{1,c}{1,1}(n).ActionlickTimeRight(j)+0.5)],[(n-0.5),(n+0.5)],'color','r','linewidth',2)
        %plot(licktime_hc_Opto(timeonset_inds_hc_Opto(n)).ActionlickTimeRight(j),n,'b.','markersize',4);
        hold on
    end
       for j = 1:length(lick_Opto{1,c}{1,1}(n).ActionlickTimeLeft)
        line([(lick_Opto{1,c}{1,1}(n).ActionlickTimeLeft(j)-0.5),(lick_Opto{1,c}{1,1}(n).ActionlickTimeLeft(j)+0.5)],[(n-0.5),(n+0.5)],'color','k','linewidth',2)
        %plot(licktime_hc_Opto(timeonset_inds_hc_Opto(n)).ActionlickTimeLeft(j),n,'r.','markersize',4);
        hold on
     end
    
     patch([0,0,300,300],[0,length(lick_Opto{1,c}{1,1}),length(lick_Opto{1,c}{1,1}),0],[0.5,0.5,0.5],'linestyle','none','facealpha',0.5);
   %line([(0-0.5),(0+0.5)],[(n-0.5),(n+0.5)],'color','k','linewidth',3);
    %patch([(time_OptoOnset(n)-time_stimOnset(n)),(time_OptoOnset(n)-time_stimOnset(n)),(time_OptoOff(n)-time_stimOnset(n)),(time_OptoOff(n)-time_stimOnset(n))],[(n-0.5),(n+0.5),(n+0.5),(n-0.5)],'b','linestyle','none','facealpha',0.5);
    %patch([SessionResults{trial_num_opto_hc(n)}.Time_optoStimOnset,SessionResults{trial_num_opto_hc(n)}.Time_optoStimOnset,SessionResults{trial_num_opto_hc(n)}.Time_optoStimOffTime,SessionResults{trial_num_opto_hc(n)}.Time_optoStimOffTime],[(n-0.5),(n+0.5),(n+0.5),(n-0.5)],[0.5,0.5,0.5],'linestyle','none','facealpha',0.5);

end
xlim([-1000,6000]);
%ylim([0,length(lick_Opto{1,c}{1,1})]);
title(fig_name_Opto{c});
end


filename_2 = strrep(filename,'.mat','_Opto_lickRaster');

filename_4 = strsplit(filename_2,'_');
filename_sup = [filename_4{1} '-' filename_4{2} '-Opto-lickRaster'];
suptitle(filename_sup);

set(gcf,'PaperUnits','inches');
x_width = 7;
y_width = 18;
set(gcf,'PaperPosition',[0 0 x_width y_width]);

saveas(gcf,['D:\behavior_data\temp\lick\' filename_2 '.png']);
saveas(gcf,['D:\behavior_data\temp\lick\' filename_2 '.fig']);
close all;
% lick_nonOpto = {licktime_c_nonOpto_delta{1},licktime_w_nonOpto_delta{1},licktime_c_nonOpto_delta{2},licktime_w_nonOpto_delta{2},licktime_c_nonOpto_delta{3},licktime_w_nonOpto_delta{3},licktime_c_nonOpto_delta{4},licktime_w_nonOpto_delta{4},licktime_c_nonOpto_delta{5},licktime_w_nonOpto_delta{5},licktime_c_nonOpto_delta{6},licktime_w_nonOpto_delta{6},licktime_c_nonOpto_delta{7},licktime_w_nonOpto_delta{7},licktime_c_nonOpto_delta{8},licktime_w_nonOpto_delta{8}};
% fig_name_nonOpto = {'7k nonOpto correct','7k nonOpto wrong','10k nonOpto correct','10k nonOpto wrong','11k nonOpto correct','11k nonOpto wrong','12k nonOpto correct','12k nonOpto wrong','16k nonOpto correct','16k nonOpto wrong','18k nonOpto correct','18k nonOpto wrong','20k nonOpto correct','20k nonOpto wrong','28k nonOpto correct','28k nonOpto wrong'};
%len = {lenRtime_hc_Opto,lenLtime_hc_Opto,lenRtime_hw_Opto,lenLtime_hw_Opto,lenRtime_lc_Opto,lenLtime_lc_Opto,lenRtime_lw_OptolenLtime_lw_Opto,};
for h = 1:length(freq_type)
lick_nonOpto{1,2*h-1} = {licktime_c_nonOpto_delta{h}};
lick_nonOpto{1,2*h} = {licktime_w_nonOpto_delta{h}};
fig_name_nonOpto{(2*h-1)} = {sprintf('%d nonOpto correct',freq_type(h))};
fig_name_nonOpto{(2*h)} = {sprintf('%d nonOpto wrong',freq_type(h))};
end

for b = 1:2*(length(freq_type))
subplot((length(freq_type)),2,b)
%title('28000 correct trials')
for n=1:length(lick_nonOpto{1,b}{1,1})
  
    for j = 1:length(lick_nonOpto{1,b}{1,1}(n).ActionlickTimeRight)
        line([(lick_nonOpto{1,b}{1,1}(n).ActionlickTimeRight(j)-0.5),(lick_nonOpto{1,b}{1,1}(n).ActionlickTimeRight(j)+0.5)],[(n-0.5),(n+0.5)],'color','r','linewidth',2)
        %plot(licktime_hc_nonOpto(timeonset_inds_hc_nonOpto(n)).ActionlickTimeRight(j),n,'b.','markersize',4);
        hold on
    end
       for j = 1:length(lick_nonOpto{1,b}{1,1}(n).ActionlickTimeLeft)
        line([(lick_nonOpto{1,b}{1,1}(n).ActionlickTimeLeft(j)-0.5),(lick_nonOpto{1,b}{1,1}(n).ActionlickTimeLeft(j)+0.5)],[(n-0.5),(n+0.5)],'color','k','linewidth',2)
        %plot(licktime_hc_nonOpto(timeonset_inds_hc_nonOpto(n)).ActionlickTimeLeft(j),n,'r.','markersize',4);
        hold on
     end
    
     patch([0,0,300,300],[0,length(lick_nonOpto{1,b}{1,1}),length(lick_nonOpto{1,b}{1,1}),0],[0.5,0.5,0.5],'linestyle','none','facealpha',0.5);
   %line([(0-0.5),(0+0.5)],[(n-0.5),(n+0.5)],'color','k','linewidth',3);
    %patch([(time_OptoOnset(n)-time_stimOnset(n)),(time_OptoOnset(n)-time_stimOnset(n)),(time_OptoOff(n)-time_stimOnset(n)),(time_OptoOff(n)-time_stimOnset(n))],[(n-0.5),(n+0.5),(n+0.5),(n-0.5)],'b','linestyle','none','facealpha',0.5);
    %patch([SessionResults{trial_num_opto_hc(n)}.Time_optoStimOnset,SessionResults{trial_num_opto_hc(n)}.Time_optoStimOnset,SessionResults{trial_num_opto_hc(n)}.Time_optoStimOffTime,SessionResults{trial_num_opto_hc(n)}.Time_optoStimOffTime],[(n-0.5),(n+0.5),(n+0.5),(n-0.5)],[0.5,0.5,0.5],'linestyle','none','facealpha',0.5);

end
xlim([-1000,6000]);
%ylim([0,length(lick_nonOpto{1,b}{1,1})]);
title(fig_name_nonOpto{b});
end


filename_3 = strrep(filename,'.mat','_nonOpto_lickRaster');
filename_5 = strsplit(filename_3,'_');
filename_sup = [filename_5{1} '-' filename_5{2} '-nonOpto-lickRaster'];
suptitle(filename_sup);

set(gcf,'PaperUnits','inches');
x_width = 7;
y_width = 18;
set(gcf,'PaperPosition',[0 0 x_width y_width]);
saveas(gcf,['D:\behavior_data\temp\lick\' filename_3 '.png']);
saveas(gcf,['D:\behavior_data\temp\lick\' filename_3 '.fig']);
close all;
clear all;
file = dir('D:\behavior_data\temp\*.mat');
end