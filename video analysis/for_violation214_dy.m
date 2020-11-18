%根据led on和trial start之间的差值找到trial start的点
dbstop if error
close all
clear

path='D:\xulab\behavior\video tracking\example_tracking_licking_without_lickport';
datapath='F:\FP\pyx214_20190802';
led_file=[path,filesep,'pyx214_led.csv'];
file_trace=[datapath,filesep,'video\2019-08-02-133157-FP-pyx214DeepCut_resnet50_tracking-lickport-tongueAug24shuffle1_100000.csv'];
behdata=[datapath,filesep,'2019_08_02_pyx214-FP'];
%ITI时间差为120，最终预测trial_start_OLED为647

% led_file=[path,'\233\cropped_30scut.csv'];%2019-07-25-163638vio-pyx233.avi源视频
% file_trace=[path,'\233\2019-07-25-163638vio-pyx233-cut30s-deepcut.csv'];
% behdata=[path,'\233\pyx233_20190725'];
%需要修改的参数：ITI时间差为60，最终预测trial_start_OLED为743（30scut）

%导入数据
T_trace=readtable(file_trace,'HeaderLines',2);
T=readtable(led_file);
load(behdata);
fr=24;
%根据荧光亮度变化获得视频中trial开始
% if exist([path,filesep,'trialStartVideoByOLEDOn.mat'],'file')
%     load([path,filesep,'trialStartVideoByOLEDOn.mat']);
% else
    [ trial_start_OLED ] = fTrialStartByOLED_dy( T, fr ,path);
% end
temp_trial_start_OLED=[trial_start_OLED;length(T.X)/fr];
    
%behdata处理：将检测到的lick按照trial画出来，类似实验过程中绘出来的图
% figure_lick=figure;
% subplot(3,2,1);
% set(gcf,'Position',[100,100,600,800]);
% height=1;
%% plot licking using arduino data
licktime_arduino=cell(3,length(SessionResults));
for i=1:length(SessionResults)
    left=strsplit(SessionResults{i}.Action_lickTimeLeft,'|');
    left(:,1)=[];
    lickleft=str2double(left)/1000;
%     for d=1:length(lickleft)
%         dot_left=plot([lickleft(d),lickleft(d)],[i*height,(i+1)*height],'b','LineWidth',2);
%         hold on
%     end
    right=strsplit(SessionResults{i}.Action_lickTimeRight,'|');
    right(:,1)=[];
    lickright=str2double(right)/1000;
%     for d=1:length(lickright)
%         dot_right=plot([lickright(d),lickright(d)],[i*height,(i+1)*height],'r','LineWidth',2);
%         hold on
%     end
    licktime_arduino{1,i}=lickright;
    licktime_arduino{2,i}=lickleft;
    licktime_arduino{3,i}=[lickleft,lickright];
end
% legend([dot_left,dot_right],{'left lick','right lick'});
% xlabel('time(s)');
% title('lick detected by arduino');
% set(gca,'Xlim',[0,12]);
% set(gca,'FontSize',14);
%%
%舌头到水杆距离预测lick
likelihood_portR=min(T_trace.likelihood_3);
likelihood_portL=min(T_trace.likelihood_4);

% 左右port可靠点占比
rows_right=T_trace.likelihood_3>0.95;
q=sum(rows_right);
reliable_r=q/length(rows_right);

rows_left=T_trace.likelihood_4>0.95;
p=sum(rows_left);
reliable_l=p/length(rows_left);

%由于port处于固定位置，故选择预测位置分布最多的点作为port的坐标
figure_reliableport=figure;
scatter(T_trace.x_3,T_trace.y_3,'b');
hold on
x0_r=fmode(T_trace.x_3);
y0_r=fmode(T_trace.y_3);
scatter(x0_r,y0_r,'r');

scatter(T_trace.x_4,T_trace.y_4,'y');
hold on
x0_l=fmode(T_trace.x_4);
y0_l=fmode(T_trace.y_4);
scatter(x0_l,y0_l,'g');
title('location of lickports');
set(gca,'FontSize',14);

%提取捕捉到的tongue的可靠坐标
figure_tongue_to_lickports=figure;
T_trace.distance_r=sqrt((T_trace.x-x0_r).^2+(T_trace.y-y0_r).^2);
T_trace.distance_l=sqrt((T_trace.x-x0_l).^2+(T_trace.y-y0_l).^2);
subplot(3,1,1);
plot(1:length(T_trace.x),T_trace.distance_r,'r');
hold on
plot(1:length(T_trace.x),T_trace.distance_l+max(T_trace.distance_r),'b');
title('distance of tongue to lickports');
set(gca,'FontSize',14);

times_as_threshold=2;

subplot(3,1,2);
histogram(T_trace.distance_r);
threshold_distance_r=mean(T_trace.distance_r(rows_right))-std(T_trace.distance_r(rows_right))*times_as_threshold;%y轴方向划分threshold（超出范围才认为移动）
ylim=get(gca,'Ylim');
hold on
plot([threshold_distance_r,threshold_distance_r],[ylim(1),ylim(end)]);%画出threshold范围
text(threshold_distance_r,mean(ylim),num2str(threshold_distance_r));%列出数值
xlabel('distance of tongue to right lickport');
set(gca,'FontSize',14);


subplot(3,1,3);
histogram(T_trace.distance_l);
threshold_distance_l=mean(T_trace.distance_l(rows_left))-std(T_trace.distance_l(rows_left))*times_as_threshold;%y轴方向划分threshold（超出范围才认为移动）
ylim_l=get(gca,'Ylim');
hold on
plot([threshold_distance_l,threshold_distance_l],[ylim_l(1),ylim_l(end)]);%画出threshold范围
text(threshold_distance_l,mean(ylim_l),num2str(threshold_distance_l));%列出数值
xlabel('distance of tongue to left lickport');
set(gca,'FontSize',14);
%计算舌头到水嘴距离预测lick的阈值
%% predict lick by distance from lickport to tongue
%distance小于阈值的,且在此之前distance大于阈值的可能为lick
T_trace.distance_r(~rows_right)=max(T_trace.distance_r);
right_distance_lick=(T_trace.distance_r<threshold_distance_r);
right_lick_distance=(diff(right_distance_lick)>0);
right_lick_distance=[0;right_lick_distance];
fr_right_lick_distance=find(right_lick_distance);
%把不可靠的点滤除，满足条件的点（距离，连续变化）——逻辑判断，差分，find；
T_trace.distance_l(~rows_left)=max(T_trace.distance_l);
left_distance_lick=(T_trace.distance_l<threshold_distance_l);
left_lick_distance=(diff(left_distance_lick)>0);
left_lick_distance=[0;left_lick_distance];
fr_left_lick_distance=find(left_lick_distance);

time_right_lick_distance=fr_right_lick_distance/fr;
time_left_lick_distance=fr_left_lick_distance/fr;
    
% figure(figure_lick)
% hold on
% subplot(3,2,2);
% if ~isempty(time_right_lick_distance)&&~isempty(time_left_lick_distance)
    licktime_by_distance=cell(3,length(trial_start_OLED));
    %转为时间；cell储存
    for g=1:length(trial_start_OLED)
        right_cycle=time_right_lick_distance(logical((time_right_lick_distance>temp_trial_start_OLED(g)).*(time_right_lick_distance<temp_trial_start_OLED(g+1))));
        left_cycle=time_left_lick_distance(logical((time_left_lick_distance>temp_trial_start_OLED(g)).*(time_left_lick_distance<temp_trial_start_OLED(g+1))));
        licktime_by_distance{1,g}=right_cycle-trial_start_OLED(g);
        licktime_by_distance{2,g}=left_cycle-trial_start_OLED(g);
        licktime_by_distance{3,g}=[right_cycle-trial_start_OLED(g);left_cycle-trial_start_OLED(g)];
%         %将每个预测lick对应到每个trial（时间信息），cell储存相对trial的起始时间
%         for ri=1:length(right_cycle)
%             do_ri=plot([right_cycle(ri)-temp_trial_start_OLED(g),right_cycle(ri)-temp_trial_start_OLED(g)],[g*height,(g+1)*height],'r','LineWidth',2);
%             hold on
%         end
%         %画点
%         
%         for le=1:length(left_cycle)
%             do_le=plot([left_cycle(le)-temp_trial_start_OLED(g),left_cycle(le)-temp_trial_start_OLED(g)],[g*height,(g+1)*height],'b','LineWidth',2);
%             hold on
%         end
    end
%     legend([do_le,do_ri],{'left lick','right lick'});
%     xlabel('time(s)');
%     title('lick predicted by distance of tongue to lickports');
%     set(gca,'FontSize',14);
% else
%     title('no lick predicted by distance of tongue to lickports');
% end
%%  predict lick by y   
figure_y=figure;
rows=T_trace.likelihood>0.95;
histogram(T_trace.y(rows));hold on;%hold on即保留已有图片
xlabel('Coordinates');
title('y distribution');
set(gca,'FontSize',14);
thresholdy=mean(T_trace.y(rows))+std(T_trace.y(rows))*times_as_threshold;%y轴方向划分threshold（超出范围才认为移动）
ylim=get(gca,'Ylim');%？？？？？？
plot([thresholdy,thresholdy],[ylim(1),ylim(end)]);%画出threshold范围
text(thresholdy,mean(ylim),num2str(thresholdy));%列出数值
set(gca,'FontSize',14);

y0=fmode(T_trace.y);

T_trace.y(~rows)=y0;%filter the points that likelihood<0.95
t_lick=(T_trace.y>thresholdy);%using histogram to find about 3.3 may be threshold for licking
ind_lick3=(diff(t_lick)>0);
ind_lick3=[0;ind_lick3];
p3=find(ind_lick3);
licktime_by_y=p3/fr;
% figure(figure_lick); 
% hold on
% subplot(3,2,3);

licktime_by_y_trial=cell(1,length(trial_start_OLED));
for g=1:length(trial_start_OLED)
    y_cycle=licktime_by_y(logical((licktime_by_y>temp_trial_start_OLED(g)).*(licktime_by_y<temp_trial_start_OLED(g+1))));
    licktime_by_y_trial{1,g}=y_cycle-trial_start_OLED(g);    
    
%      for le=1:length(y_cycle)
%         dot_y=plot([y_cycle(le)-temp_trial_start_OLED(g),y_cycle(le)-temp_trial_start_OLED(g)],[g*height,(g+1)*height],'k','LineWidth',2);
%         hold on
%     end
end
% xlabel('time(s)');
% title('lick predicted by y of tongue');
% set(gca,'FontSize',14);

%% predict lick by distance of tongue to its origin
x0=fmode(T_trace.x);%对T的x列进行计算，确定目标的平衡位点（分布最多）？？？？？？？？？
y0=fmode(T_trace.y);%use mode as origin
T_trace.distance=sqrt((T_trace.x-x0).^2+(T_trace.y-y0).^2);%计算位移距离

figure_distance=figure;
histogram(T_trace.distance(rows));hold on;
xlabel('Distance');
title('distance distribution');
threshold=mean(T_trace.distance(rows))+std(T_trace.distance(rows));
ylim=get(gca,'Ylim');
plot([threshold,threshold],[ylim(1),ylim(end)]);%画出一条线
text(threshold,mean(ylim),num2str(threshold));%定义text的坐标及数值
set(gca,'FontSize',14);

T_trace.distance(~rows)=0;%filter the points that likelihood<0.95
t_lick=(T_trace.distance>threshold);%using histogram to find about 3.3 may be threshold for licking
ind_lick2=(diff(t_lick)>0);
ind_lick2=[0;ind_lick2];
p2=find(ind_lick2);

licktime_by_distance_tongue=p2/fr;
% figure(figure_lick);
% hold on
% subplot(3,2,4);

licktime_by_distance_trial=cell(1,length(trial_start_OLED));
for g=1:length(trial_start_OLED)
    tongue_cycle=licktime_by_distance_tongue(logical((licktime_by_distance_tongue>temp_trial_start_OLED(g)).*(licktime_by_distance_tongue<temp_trial_start_OLED(g+1))));
    licktime_by_distance_trial{1,g}=tongue_cycle-trial_start_OLED(g);    
    
%      for le=1:length(tongue_cycle)
%         dot_tongue=plot([tongue_cycle(le)-temp_trial_start_OLED(g),tongue_cycle(le)-temp_trial_start_OLED(g)],[g*height,(g+1)*height],'k','LineWidth',2);
%         hold on
%     end
end
% xlabel('time(s)');
% title('lick predicted by shift distance of tongue');
% set(gca,'FontSize',14);
%% predict lick by filtered y
%傅里叶变换滤除舌头的y噪音点，根据lick rate
figfft=figure;
S=T_trace.y;
Sf=T_trace.y(rows);
n = length(Sf);
X = fft(Sf);
f = (0:n-1)*(fr/n);     %frequency range
power = abs(X).^2/n;    %power
% plot(f,power);
Y = fftshift(X);
fshift = (-n/2:n/2-1)*(fr/n); % zero-centered frequency range
powershift = abs(Y).^2/n;     % zero-centered power
plot(fshift,powershift);
hold on;
set(gca,'Ylim',[0,10000]);
%normal fitting of licking at around 7Hz(using 6~8Hz data)
% fun = fittype('y0+(a/(w*sqrt(pi/2)))*exp(-2*((x-xc)/w).^2)', 'independent', 'x', 'dependent', 'y');
% fitrange=(fshift'>6&fshift'<8);
% [fitresult, gof] = fit(fshift(fitrange)',powershift(fitrange),fun,'Start',[1.1,1.1,1.1,1.1]);
% figure(figfft);
% powerfit=fitresult(fshift);
% plot(fitresult,'-k');
%according to previous section, I just filter the out signals not related
%to licking
wp = [6 8 ] / (fr/2);%通带截止频率
ws = [5 9 ] / (fr/2);  %阻带截止频率
alpha_p = 3; %通带允许最大衰减为  db
alpha_s = 10;%阻带允许最小衰减为  db
%获取阶数和截止频率
[ N3, wn ] = buttord( wp , ws , alpha_p , alpha_s);
%获得转移函数系数
[ b, a ] = butter(N3,wn,'bandpass');
filter_bp_sf = filter(b,a,Sf);
filter_bp_s=zeros(length(rows),1);
j=1;
for i=1:length(rows)
    if rows(i)==1
        filter_bp_s(i)=filter_bp_sf(j);
        j=j+1;
    end
end
figpose=figure;
plot((1:size(T_trace,1)),filter_bp_s);
%set a threshold for licking- method-I
std_fbp_s=std(filter_bp_s);
times=2;
aa=(filter_bp_s>times*std_fbp_s);
bb=(diff(aa)>0);%diff results will have less element,so add back
bb=[0;bb];
ind_lick=bb;
% ind_lick=(filter_bp_s==3*std_fbp_s).*(diff(filter_bp_s)>0);%choose 3 std and increase time point as licking time point
pp=find(ind_lick);

licktime_by_y_filter=pp/fr;

% figure(figure_lick);
% hold on
% subplot(3,2,5);

licktime_by_y_filter_trial=cell(1,length(trial_start_OLED));
for g=1:length(trial_start_OLED)
    y_filter_cycle=licktime_by_y_filter(logical((licktime_by_y_filter>temp_trial_start_OLED(g)).*(licktime_by_y_filter<temp_trial_start_OLED(g+1))));
    licktime_by_y_filter_trial{1,g}=y_filter_cycle-trial_start_OLED(g);    
    
%      for le=1:length(y_filter_cycle)
%         dot_y_filter=plot([y_filter_cycle(le)-temp_trial_start_OLED(g),y_filter_cycle(le)-temp_trial_start_OLED(g)],[g*height,(g+1)*height],'k','LineWidth',2);
%         hold on
%     end
end
% xlabel('time(s)');
% title('lick predicted by y filtered of tongue');
% set(gca,'FontSize',14);
%% 评价同种对齐方法，不同指标对licking的不同预测效果
%evaluate prediction based on video
% 1)compare recall and precision of different methods 在时间窗范围内寻找hit的lick
%   a)recall
setwin_ind=2;%the time window for each method to compare prediction accuracy
figEvaluate=figure;
subplot(2,2,1);
win=0.1:0.1:1;

recall_distance_trial=zeros(1,length(win));
recall_y=zeros(1,length(win));
recall_y_filter=zeros(1,length(win));
for i=1:length(win)
    window=[win(i),win(i)];
%     recall_distance_r(i)=frecall(licktime_by_distance(1,:),licktime_arduino(1,4:size(licktime_arduino,2)),window);
%     recall_distance_l(i)=frecall(licktime_by_distance(2,:),licktime_arduino(2,4:size(licktime_arduino,2)),window);
    recall_y_filter(i)=frecall(licktime_by_y_filter_trial(1,1:size(licktime_arduino,2)),licktime_arduino(3,:),window);
    recall_y(i)=frecall(licktime_by_y_trial(1,1:size(licktime_arduino,2)),licktime_arduino(3,:),window);
    recall_distance_trial(i)=frecall(licktime_by_distance_trial(1,1:size(licktime_arduino,2)),licktime_arduino(3,:),window);
end
% plot(recall_distance_r,'b-');hold on;
% plot(recall_distance_l,'r-');
plot(recall_y_filter,'k-');hold on;
plot(recall_y,'m-');
plot(recall_distance_trial,'g-');

ylim=get(gca,'Ylim');
plot([setwin_ind,setwin_ind],[ylim(1),ylim(end)],'k-');
legend('prediction by y filtered of tongue','prediction by y of tongue','prediction by location shift of tongue');
set(gca,'XTick',1:length(win),'XTickLabel',num2str(win'));
ylabel('recall');
xlabel('time window(s)');
set(gca,'FontSize',14);

subplot(2,2,2);
% bar(1,recall_distance_r(setwin_ind),'FaceColor','b');hold on;
% bar(2,recall_distance_l(setwin_ind),'FaceColor','r');
bar(1,recall_y_filter(setwin_ind),'FaceColor','k');hold on;
bar(2,recall_y(setwin_ind),'FaceColor','m');
bar(3,recall_distance_trial(setwin_ind),'FaceColor','g');

set(gca,'XTick',[]);
set(gca,'FontSize',14);

%   b)precision
subplot(2,2,3);
win=0.1:0.1:1;

for i=1:length(win)
    window=[win(i),win(i)];
%     recall_distance_r(i)=frecall(licktime_arduino(1,4:size(licktime_arduino,2)),licktime_by_distance(1,:),window);
%     recall_distance_l(i)=frecall(licktime_arduino(2,4:size(licktime_arduino,2)),licktime_by_distance(2,:),window);

    recall_y_filter(i)=frecall(licktime_arduino(3,:),licktime_by_y_filter_trial(1,1:size(licktime_arduino,2)),window);
    recall_y(i)=frecall(licktime_arduino(3,:),licktime_by_y_trial(1,1:size(licktime_arduino,2)),window);
    recall_distance_trial(i)=frecall(licktime_arduino(3,:),licktime_by_distance_trial(1,1:size(licktime_arduino,2)),window);
end
% plot(recall_distance_r,'b-');hold on;
% plot(recall_distance_l,'r-');
plot(recall_y_filter,'k-');hold on;
plot(recall_y,'m-');
plot(recall_distance_trial,'g-');

ylim=get(gca,'Ylim');
plot([setwin_ind,setwin_ind],[ylim(1),ylim(end)],'k-');
legend('prediction by y filtered of tongue','prediction by y of tongue','prediction by location shift of tongue');
set(gca,'XTick',1:length(win),'XTickLabel',num2str(win'));
ylabel('precision');
xlabel('time window(s)');
set(gca,'FontSize',14);

subplot(2,2,4);
% bar(1,recall_distance_r(setwin_ind),'FaceColor','b');hold on;
% bar(2,recall_distance_l(setwin_ind),'FaceColor','r');
bar(1,recall_y_filter(setwin_ind),'FaceColor','k');hold on;
bar(2,recall_y(setwin_ind),'FaceColor','m');
bar(3,recall_distance_trial(setwin_ind),'FaceColor','g');
set(gca,'XTick',[]);
set(gca,'FontSize',14);
%% plot lick rate around real licks(see if time shifts for prediction)
fig_histo_lickrate=figure;
timeMin=1000;
binStep=50;
trialLength=3000;
binSize=500;

[ ~,lickRate_distance_tongue ] = fLickRateAlignedCell(licktime_by_distance_trial(1,1:size(licktime_arduino,2)),licktime_arduino(3,:),binSize,trialLength,timeMin, binStep);
% [ ~,lickRate_distance_port ] = fLickRateAligned( licktime_by_distance(3,:), licktime_arduino(3,4:end),binSize,trialLength,timeMin, binStep);
[ ~,lickRate_y ] = fLickRateAlignedCell( licktime_by_y_trial(1,1:size(licktime_arduino,2)),licktime_arduino(3,:),binSize,trialLength,timeMin, binStep);
[ ~,lickRate_y_filter ] = fLickRateAlignedCell( licktime_by_y_filter_trial(1,1:size(licktime_arduino,2)),licktime_arduino(3,:),binSize,trialLength,timeMin, binStep);
ts=-timeMin:binStep:trialLength-timeMin-binSize;
plot(ts',lickRate_y);hold on;
plot(ts',lickRate_y_filter);
plot(ts',lickRate_distance_tongue);
% plot(ts',lickRate_distance_port);

xlabel('Time from true lick(ms)');
ylabel('Lick rate(/s)');
legend('y of tongue','y filtered of tongue','distance of tongue');
set(gca,'FontSize',14,'FontName','Arial');
%% evaluate prediction based on video
% 2) trial type prediction

n_trial=size(licktime_arduino,2);

time_stim_onset=cellfun(@(x) double(x.Time_stimOnset)/1000,SessionResults);
time_delay_off=cellfun(@(x) (double(x.Time_stimOffset)+double(x.Delay_duration))/1000,SessionResults);

vio_y=fvio(licktime_by_y_trial,time_stim_onset,time_delay_off);
vio_distance_tongue=fvio(licktime_by_distance_trial,time_stim_onset,time_delay_off);
vio_distance_port=fvio(licktime_by_distance,time_stim_onset,time_delay_off);
vio_y_filter=fvio(licktime_by_y_filter_trial,time_stim_onset,time_delay_off);
vio_arduino_lick=fvio(licktime_arduino(3,:),time_stim_onset,time_delay_off);%violation trial predicted by 
vio_arduino=cellfun(@(x) x.Action_choice==3,SessionResults);
% vio_arduino(1:3)=[];

nh=1:n_trial;
figVio=figure;
subplot(2,2,1);
% dot_vio_arduino=scatter(nh(vio_arduino),ones(length(nh(vio_arduino)),1));
% hold on

legendstr={'arduino','lick by arduino','y of tongue','y filtered of tongue','distance of tongue','distance from tongue to lickport'};
dot_vio=fviodot(legendstr,vio_arduino,vio_arduino_lick,vio_y,vio_y_filter,vio_distance_tongue,vio_distance_port);%五种预测方法的vio
vio_predict=cell(1,5);
vio_predict{1}=vio_arduino_lick;
vio_predict{2}=vio_y;
vio_predict{3}=vio_y_filter;
vio_predict{4}=vio_distance_tongue;
vio_predict{5}=vio_distance_port;
[vio_true{1:5}] =deal(vio_arduino);

[vio_predict_TP,vio_predict_TN]=cellfun(@fTPTN, vio_true,vio_predict);

figure(figVio);
subplot(2,2,2);
for i=1:length(vio_predict_TP)
    bar(i,vio_predict_TP(1));hold on;
end
legend(legendstr(2:end));
set(gca,'XTick',[]);
title('violation predicted correctly');
set(gca,'FontSize',14);

subplot(2,2,3);
for i=1:length(vio_predict_TN)
    bar(i,vio_predict_TN(1));hold on;
end
title('no violation predicted correctly');
set(gca,'XTick',[]);
set(gca,'FontSize',14);

subplot(2,2,4);
scatter(1,sum(vio_arduino)/n_trial,'b');
hold on
for i=1:length(vio_predict)
    bar(i+1,sum(vio_predict{i})/n_trial);
end
legend(legendstr);
title('violation rate');
set(gca,'XTick',[]);
set(gca,'FontSize',14);
%%
%评价同种预测的对齐效果

%%
function [f0]=fmode(f)%找到最大的分布值区间？
[N,edges,bin] = histcounts(f,100);%将f分成100个大小均一的bin
f0 = edges(N == max(N));%内含最多元素的bin的edge
if length(f0)==2
    f0=mean(f0);
end
end

function [vio]=fvio(licktime,stimulus_onset,delay_off)
%based on licking time in a trial to decide trial type (violation or not) 
%Input-
%   licktime- 1-by-n cell array,each cell represent a trial
%   stimulus_onset,delay_off- vector, each element a trial
n_trial=min([length(licktime),length(stimulus_onset),length(delay_off)]);%use the smallest size
vio=false(n_trial,1);
for i=1:n_trial
    temp=logical(licktime{i}>=stimulus_onset(i)).*(licktime{i}<=delay_off(i));
    if ~isempty(temp)
        vio(i)=(sum(temp)>0);
    end
end
end

function [dot_vio]=fviodot(legendstr,varargin)
for i=1:length(varargin)
    vio=varargin{i};
    nh=1:length(vio);
    dot_vio(i)=scatter(nh(vio),ones(length(nh(vio)),1)*i,50);
    hold on
end
legend(dot_vio,legendstr);
end

function [TP,TN]=fTPTN(groundtrue,data)
n_trial=min(length(groundtrue),length(data));
groundtrue=reshape(groundtrue(1:n_trial),1,[]);
data=reshape(data(1:n_trial),1,[]);
TP=sum(groundtrue.*(data==1))/sum(groundtrue);
TN=sum((groundtrue+data)==0)/sum(groundtrue);
end