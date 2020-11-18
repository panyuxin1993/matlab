%08/23/2017 %
clear all;  % for guofen's data
close all;

% 读取Excel文件信息


xlsname = 'G:\brain_slice_recording\direct_connection_to_INs\XUNL\SC_f_M2_vglut';%ORB PV neurons received RSP input


[num, txt]= xlsread(xlsname);

%需要分析的trail number
num0 = 0;
num1 = 0;
for i = 1:size(num,1)
    if num(i,1)==0
        num0 = num0+1;
    elseif num(i,1)==1
        num1 = num1+1;
    end
end
%给同一个神经元的其他trial赋值地址
for b = 2 : size(txt,1)
   if length(cell2mat(txt(b,1))) == 0
   txt(b,1) = txt(b-1,1);
   end
end

% for the final summary
num_0 = 1;
num_1 = 1;

%从每种处理的第一个开始分析
for bb = 2:size(num,1)+1
       tmp = cell2mat(txt(bb,2));


        current_folder_name = ['G:\brain_slice_recording\direct_connection_to_INs\XUNL'];


%        datapath = horzcat(current_folder_name, tmp(48:end),'\'); % Cg_VIP_f_MD_EPSP_single.xls'; %Cg VIP neurons received MD input
%        current_folder_name = ['E:\brain_slice_recording_SH\Inject_V1_Record_Cg_beads_in_V1'];
%        datapath = horzcat(current_folder_name, tmp(48:end),'\'); % Cg_VIP_f_MD_EPSP_single.xls'; %Cg VIP neurons received MD input
       datapath = horzcat(current_folder_name, '\', tmp,'\'); %Cg VIP neurons received ORB input
%        datapath = horzcat(current_folder_name, tmp(44:end),'\');
       abffilename_E_tmp = cell2mat(txt(bb,4)); 
       if abffilename_E_tmp(end-2:end) ~= 'abf';
           abffilename_E = [abffilename_E_tmp, '.abf'];  %赋值需要读取的abf文件 excitatory
           abffilename_E_train_tmp = cell2mat(txt(bb,5)); 
           abffilename_E_train = [abffilename_E_train_tmp, '.abf'];
       else
           abffilename_E = abffilename_E_tmp;
           abffilename_E_train_tmp = cell2mat(txt(bb,5)); 
           abffilename_E_train = abffilename_E_train_tmp;  %赋值需要读取的abf文件 inhibitory
       end
     %%%%%%%%%%%  E component  %%%%%%%%%%%%%%%%  
%        % Load the ABFfile
       [This_data,Sampling_Interval,Headers]=abfload([datapath abffilename_E]);
       % Store information from the ABFfile
       Sample_Rate=10^6/Sampling_Interval; % Sample Rate per second
       Sweep_Time_sec= Headers.dataPtsPerChan/Sample_Rate; % Time of each sweep
       Num_Sweeps=size(This_data,3); % Num of Sweeps
       Channels=size(This_data,2);
       % analysis
       st=1.8; % 画图和算charge的起始时间 Start of the Pulse in Sec
       et=2.4; % 画图和算charge的结束时间 End of the Pulse in Sec
       pulse_number = 1;
        
      % Determine the peak response (import averaged trace)
       for Sweeps=1:Num_Sweeps
           Voltage_Trace=This_data(st*Sample_Rate:et*Sample_Rate,1,Sweeps);
           laser_duration = find(This_data(:,2)>3);
           laser_start = laser_duration(1);
           STD_baseline = std(This_data(laser_start-3000:laser_start));  %取刺激开始前的一段计算baseline的standard deviation
    
           time_start = laser_start;
           for i = 1:pulse_number
              [Peak_E,PI_E]=max(This_data(time_start:time_start+1*Sample_Rate));  %找分析区间的最大值 从laser_start开始后1s之内
              if abs(Peak_E) < STD_baseline
                Peaks_E(i) = 0;
              else
                Peaks_E(i) = abs(Peak_E);
              end
              PIs_E(i) = PI_E;  %PI_E是出现最小值的index
           end 
           %%%% charge %%%%%%
           Charge_E = abs(sum(Voltage_Trace));           
       end

%figure
figure
subplot(1,2,1)
x = [1:length(Voltage_Trace)]/Sample_Rate*1000; %ms
plot(x, Voltage_Trace)
xlabel('Time (ms)')
ylabel('Amplitude (pA)')
% ylim([-2 26])
xlim([0 600])
title(tmp)
subplot(1,2,2)
x2 = 1:pulse_number;
plot(x2,Peaks_E,'ko-')
% ylim([-0.1 1.3])
xlim([0.5 1.5])
ylabel('Amplitude (pA)')

% saveppt('Cg_Pyr.ppt')
       
       
%%%%%%%%%%%  EPSP PPR 10Hz  %%%%%%%%%%%%%%%%  
    if abffilename_E_train(end) == 'f'; 
        % Load the ABFfile
       [This_data,Sampling_Interval,Headers]=abfload([datapath abffilename_E_train]);
       % Store information from the ABFfile
       Sample_Rate=10^6/Sampling_Interval; % Sample Rate per second
       Sweep_Time_sec= Headers.dataPtsPerChan/Sample_Rate; % Time of each sweep
       Num_Sweeps=size(This_data,3); % Num of Sweeps
       Channels=size(This_data,2);
       % analysis
       st=1.6; % Start of the Pulse in Sec
       et=3.4; % End of the Pulse in Sec
       pulse_number = 10;
       sti_frequency = 10;

      % Determine the peak response (import averaged trace)
       for Sweeps=1:Num_Sweeps
           Voltage_Trace=This_data(st*Sample_Rate:et*Sample_Rate,1,Sweeps);
           laser_duration = find(This_data(:,2)>3);
           laser_start = laser_duration(1);
           STD_baseline = std(This_data(laser_start-3000:laser_start));

           time_start = laser_start;
           for i = 1:pulse_number
              [Peak_E_train,PI_E_train]=max(This_data(time_start:time_start+1/sti_frequency*Sample_Rate));
              time_start = time_start+0.1*Sample_Rate;
              if Peak_E_train < STD_baseline
                Peaks_E_train(i) = 0;
              else
                Peaks_E_train(i) = Peak_E_train;
              end
              PIs_E_train(i) = PI_E_train;
              i = i+1;
           end 
           PPR = Peaks_E_train(2)/Peaks_E_train(1); % calculate pair pulse ratio
           %%%% charge %%%%%%
           Charge_E_train = sum(Voltage_Trace);
       end

        %figure
        figure
        subplot(1,3,1)
        x = [1:length(Voltage_Trace)]/Sample_Rate*1000; %ms
        plot(x, Voltage_Trace)
        xlabel('Time (ms)')
        ylabel('Amplitude (pA)')
        % ylim([-2 26])
        xlim([0 1800])
        title(tmp)
        subplot(1,3,2)
        x2 = 1:10;
        plot(x2,Peaks_E_train,'ko-')
        % ylim([-0.1 1.3])
        xlim([0.5 10.5])
        ylabel('Amplitude (pA)')
        subplot(1,3,3)
        x3 = 1;
        plot(x3, PPR, 'ro')
    else
      Peaks_E_train = 0;
      Charge_E = 0;
      PPR = 0;
    end

%     saveppt('Cg_Pyr.ppt')

%%%%%%% calculation finished, start data accumulation
% num_0 = 1 and num_1 = 1 were set at the beginning
%        if num(bb,1) == 0  %Cg_VIP_f_MD_EPSP_single.xls and Cg_VIP_f_ORB_EPSP_single.xls' excel file num has 'NaN' in the first row
       if num(bb-1,1) == 0
           Peak_E_all0(num_0) = Peak_E;
           Peaks_E_all0(num_0,:) = Peaks_E_train;
           Charge_E_all0(num_0) = Charge_E;
           PPR_E_all0(num_0) = PPR;
           num_0 = num_0+1;
%         elseif num(bb,1) == 1  %Cg_VIP_f_MD_EPSP_single.xls and Cg_VIP_f_ORB_EPSP_single.xls'this excel file num has 'NaN' in the first row
        elseif num(bb-1,1) == 1
           Peak_E_all1(num_1) = Peak_E;
           Peaks_E_all1(num_1,:) = Peaks_E_train;
           Charge_E_all1(num_1) = Charge_E;
           PPR_E_all1(num_1) = PPR;
           num_1 = num_1+1;       
       end
end
    
%%%%% remove the unwanted items
%find the index for items to be removed
% ii = 1;
% for i = 1:size(Peaks_I_all0)
%     if Peaks_I_all0(i) == -10000
%         remove_idx(ii) = i;
%         ii = ii+1;
%         i = i +1;
%     end
% end
% %remove
% Peaks_I_all0(remove_idx) = [];
% Peaks_I_all1(remove_idx) = [];
% Charge_I_all0(remove_idx) = [];
% Charge_I_all1(remove_idx) = [];


% %%%%% normalize by beads+ cells 
%    %%%%%%%%  charge 
% Norm_Charge_E0 = ones(size(Charge_E_all0));
% Norm_Charge_I0 = ones(size(Charge_I_all0));
% for i = 1:length(Charge_E_all0)
%    Norm_Charge_E1(i) = Charge_E_all1(i)/Charge_E_all0(i);   
% end
% 
% for i = 1:length(Charge_I_all0)
%    Norm_Charge_I1(i) = Charge_I_all1(i)/Charge_I_all0(i);
%    if Norm_Charge_I1(i) < 0
%        Norm_Charge_I1(i) = 0;
%    end
% end
% 
% figure
% x = [0.5,2];
% % 兴奋性charge
% y = [mean(Norm_Charge_E0),mean(Norm_Charge_E1)];
% plot(x,y,'ro-')
% hold on
% e1 = std(Norm_Charge_E1) / sqrt(length(Charge_E_all1));
% h1 = errorbar(2, mean(Norm_Charge_E1), e1, 'ro-', 'linewidth', 1.5);
% % 抑制性charge
% y2 = [mean(Norm_Charge_I0),mean(Norm_Charge_I1)];
% plot(x,y2,'bo-')
% hold on
% e2 = std(Norm_Charge_I1) / sqrt(length(Charge_I_all1));
% h2 = errorbar(2, mean(Norm_Charge_I1), e2, 'bo-', 'linewidth', 1.5);
% ylabel('Normalized charge')
% ylim([0 1.5])
% xlim([0 2.5])
% 
% [H1,P1] = ttest(Norm_Charge_E1,1);
% [H2,P2] = ttest(Norm_Charge_I1,1);
% 
% 
% % saveppt('Cg_Pyr.ppt')
% s = ['V1_Pyr_E_I_charge_evoked by Cg sti'];
% % s = ['Cg_Pyr_E_I_charge_evoked by V1 sti'];
% saveas(gcf,s,'pdf');
% 
% %%%%%%%% peak
% Norm_Peak_E0 = ones(size(Peaks_E_all0,1),1);
% Norm_Peak_I0 = ones(size(Peaks_I_all0,1),1);
% for i = 1:size(Peaks_E_all0,1)
%    Norm_Peak_E1(i) = Peaks_E_all1(i, 1)/Peaks_E_all0(i, 1);   
% end
% 
% for i = 1:size(Peaks_I_all0,1)
%    Norm_Peak_I1(i) = Peaks_I_all1(i, 1)/Peaks_I_all0(i, 1);   
% end
% 
% figure
% x = [0.5,2];
% % 兴奋性charge
% y = [mean(Norm_Peak_E0),mean(Norm_Peak_E1)];
% plot(x,y,'ro-')
% hold on
% e1 = std(Norm_Peak_E1) / sqrt(length(Peaks_E_all1));
% h1 = errorbar(2, mean(Norm_Peak_E1), e1, 'ro-', 'linewidth', 1.5);
% % 抑制性charge
% y2 = [mean(Norm_Peak_I0),mean(Norm_Peak_I1)];
% plot(x,y2,'bo-')
% hold on
% e2 = std(Norm_Peak_I1) / sqrt(length(Peaks_I_all1));
% h2 = errorbar(2, mean(Norm_Peak_I1), e2, 'bo-', 'linewidth', 1.5);
% ylabel('Normalized Peak')
% ylim([0 1.5])
% xlim([0 2.5])
% 
% 


% % saveppt('Cg_Pyr.ppt')
% s = ['V1_Pry_E_I_peak_evoked by Cg sti'];
% % s = ['Cg_Pry_E_I_peak_evoked by V1 sti'];
% saveas(gcf,s,'pdf');

% plot PPR
[H1,P1] = ttest2(PPR_E_all0, PPR_E_all1);
[H2,P2] = ttest(PPR_E_all0, 1);
[H3,P3] = ttest(PPR_E_all1, 1);

figure
x = [0.5,2.5];
% PPR
y = [mean(PPR_E_all0),mean(PPR_E_all1)];
plot(x,y,'bo')
hold on
e1 = std(PPR_E_all0) / sqrt(length(PPR_E_all0));
h1 = errorbar(0.5, mean(PPR_E_all0), e1, 'bo-', 'linewidth', 1.5);
hold on
e2 = std(PPR_E_all1) / sqrt(length(PPR_E_all1));
h2 = errorbar(2.5, mean(PPR_E_all1), e2, 'bo-', 'linewidth', 1.5);
hold on
plot(1, PPR_E_all0, 'bo')
hold on
plot(2, PPR_E_all1, 'bo')
ylabel('Pair pulse ratio')
ylim([0 4])
xlim([0 3])
s = ['RSP VIP and Pyr PPR-V1 input'];
% s = ['Cg_Pyr_E_I_peak_evoked by V1 sti_all points'];
saveas(gcf,s,'pdf');

% % from excel Cg VIP and pyr from MD
% Normalized_Peak_E_0 = [1.137791586, 0.86220877, 0.348052356, 1.255944557, 1.609740038, 0.786263358, ...
% 1.076028706, 0.923971595, 0.874352982, 1.125647532];
% Normalized_Peak_E_1 = [0.474491791, 2.568631883, 0.792300718, 0.198388394, 0.268355449,...
% 0.775155187, 0.40035659,0.247434707, 0.530823644, 0.298899206, 0.184273203, 0.673090064,...
% 1.325209545];
% 
% [H2, P2] = ttest2(Normalized_Peak_E_0, Normalized_Peak_E_1);

% % from excel Cg VIP and pyr from ORB
% Normalized_Peak_E_0 = [1, 0.3388, 1.6611,1];
% Normalized_Peak_E_1 = [1.371712458, 1.178530597, 1.21530357, 2.888328811, 3.353105235,1.108831052];
% [H2, P2] = ttest2(Normalized_Peak_E_0, Normalized_Peak_E_1);

%存储mat文件
filepath = 'G:\brain_slice_recording\direct_connection_to_INs\VIP\summary\final_results\';
savematfile = 'RSP VIP and Pyr PPR-V1 input.mat'
% savematfile = 'Beads_Cg_Pyr_E_I_evoked_by_V1_sti.mat'
save([filepath savematfile], 'Peak_E_all0', 'Peak_E_all1', 'Peaks_E_all0', 'Peaks_E_all1','PPR_E_all0', 'PPR_E_all1')

Peak_E_all0 = Peak_E_all0';
Peak_E_all1 = Peak_E_all1';


