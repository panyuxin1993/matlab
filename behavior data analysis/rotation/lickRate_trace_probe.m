file = dir('F:\behavior_data\temp\*.mat');
for f = 1:length(file)
filename = file(f).name;
load(['F:\behavior_data\temp\' filename])
binSize = 100;
trialLength = 6000;
sampleNum = trialLength/binSize +1;

%stimOnset = get_stimOnset(SessionResults);
   
[freq_type,nonOptoLickRate_c_r,nonOptoLickRate_c_l,nonOptoLickRate_w_r,nonOptoLickRate_w_l,OptoLickRate_c_r,OptoLickRate_c_l,OptoLickRate_w_r,OptoLickRate_w_l]  = get_lick_rate_probe(SessionResults, binSize,trialLength,sampleNum);


for g = 1:length(freq_type)
% lick_Opto{1,(2*g-1)} = {licktime_c_Opto_delta{g}};
% lick_Opto{1,(2*g)} = {licktime_w_Opto_delta{g}};
fig_name{(2*g-1)} = {sprintf('%d correct',freq_type(g))};
fig_name{(2*g)} = {sprintf('%d wrong',freq_type(g))};
end
   


lickrate_nr = {nonOptoLickRate_c_r,nonOptoLickRate_w_r};
lickrate_or = {OptoLickRate_c_r,OptoLickRate_w_r};
lickrate_nl = {nonOptoLickRate_c_l,nonOptoLickRate_w_l};
lickrate_ol = {OptoLickRate_c_l,OptoLickRate_w_l};

for n = 1:length(freq_type)*2
subplot(length(freq_type),2,n);
plot1 = plot([-1000:binSize:(trialLength-1000)],lickrate_nr{2-rem(n,2)}{((n+rem(n,2))/2)},'color','k','linewidth',2);
hold on 
plot2 = plot([-1000:binSize:(trialLength-1000)],lickrate_or{2-rem(n,2)}{((n+rem(n,2))/2)},'color','b','linewidth',2);
hold on 
plot3 = plot([-1000:binSize:(trialLength-1000)],lickrate_nl{2-rem(n,2)}{((n+rem(n,2))/2)},'--k','linewidth',2);
hold on
plot4 = plot([-1000:binSize:(trialLength-1000)],lickrate_ol{2-rem(n,2)}{((n+rem(n,2))/2)},'--b','linewidth',2);
patch([0,0,300,300],[0,1,1,0],[0.5,0.5,0.5],'linestyle','none','facealpha',0.5);
title(fig_name{n},'fontWeight','bold');

xlabel('Time');
ylabel('Lick rate');
end

filename_2 = strrep(filename,'.mat','_lickTrace');
filename_3 = strsplit(filename_2,'_');
filename_sup = [filename_3{1} '-' filename_3{2} '-lickTrace'];

suptitle(filename_sup);
%title(strrep(filename,'.mat', ' lick'));
set(gcf,'PaperUnits','inches');
x_width = 7;
y_width = 18;
set(gcf,'PaperPosition',[0 0 x_width y_width]);

saveas(gcf,['F:\behavior_data\temp\lick\' filename_2 '.png']);
saveas(gcf,['F:\behavior_data\temp\lick\' filename_2 '.fig']);

close all;
clear all;


file = dir('F:\behavior_data\temp\*.mat');
end