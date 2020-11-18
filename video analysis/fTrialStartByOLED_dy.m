function [ trial_start_final ] = fTrialStartByOLED_dy( T, framerate ,path)
%FTRIALSTARTBYOLED_DY using OLED brightness to predict trial start in the
%video, revised from Yun Du. 
%Input-
%   T- table read from the csv file containing brightness change extracted
%   from the .tiff file that transformed from raw .webm file
%   framerate- frame rate of the video, usually 24.
%   path- file path saving the result
%Output-
%   trial_start_final- n-by-1 vector, indicating the start time(s) of 
%   the n trials from session start. 

T.ts=T.X/framerate;%time stamp, using this to plot easier to be understand
figure_start=figure;%各种预测画在一张图上，人眼进行质量检测
figure(figure_start);
original=plot(T.ts,T.Y,'b');%原始LED亮度图
%led on预测
led_on=(diff(T.Y)>6);
led_on=[false;led_on];
no_led_on=sum(led_on);
x_on=T.X(led_on);%根据亮度剧变预测led on时间
y_on=T.Y(led_on);
ts_on=T.ts(led_on);
figure(figure_start);
hold on
start_by_delta=scatter(ts_on,y_on,'b');%作图,根据亮度剧变预测led on时间
%限制一：删除连续的预测点
led_continue=(diff(x_on)==1);%连续的预测点，需要被删除
led_continue=[led_continue;false];
x_continue=x_on(led_continue);
figure(figure_start);
hold on
exclusion=scatter(ts_on(led_continue),ones(length(x_continue),1)*(max(T.Y)+5),'k');%作图：去除的点（连续的预测点） 
ts_on=ts_on(~led_continue);%从预测点中去除连续点 %ts_on=setdiff(ts_on,ts_on(led_continue));
y_on=y_on(~led_continue);
x_on=x_on(~led_continue);
%对于pyx214_20190802的led至此根据delta>6及不连续出现的原则获得初步预测的668个Led亮点
led_after_continuity_exclusion=false(1,length(T.X));
for w=1:length(x_on)
    led_after_continuity_exclusion(T.X==x_on(w))=true;
end
% no_led=sum(led_after_continuity_exclusion);
% led_after_continuity_exclusion2=setdiff(x_on,x_continue);
% temp=T.X(led_after_continuity_exclusion)-led_after_continuity_exclusion2;
figure(figure_start);
hold on
on_predicted=scatter(T.ts(led_after_continuity_exclusion),T.Y(led_after_continuity_exclusion),'g');%在亮度变化线上直接画出预测点

%限制二：ITI时间相对恒定，不满足这个条件的点也要被删除
% led_interval=diff(ts_on);%led on的interval，过短的可能就是异常点，可以作为计算ITI阈值的参考
% figure_interval=figure;
% histogram(led_interval);
% xlabel('time(s)');
% title('led on interval');
%寻找灭灯时刻
delta_off=fdiff(T.Y,5);%由于led灭比亮更缓慢，故使用5个Y的差值预测off点
led_off=delta_off<-8;%this criterion is mannually set
add=false(5,1);
led_off=[add',led_off];%对预测off的统计——主要看下预测的是否准确（如果有异常点可以进一步调整条件（如dleta_off<-8)后再进行下一步）
% figure_off=figure;
% histogram(delta_off);hold on;
% title('5 Y delta');
figure(figure_start);%将off预测点和led on预测点画在一张图：便于计算ITI（每个视频存在差别）
curve_led_off=scatter(T.ts(led_off),T.Y(led_off),'c');

off_original=T.ts(led_off);
%计算预测led on在ITI范围内是否存在off——不存在则说明预测错误
hit=zeros(1,length(ts_on));
for f=1:length(ts_on)
    hit(f)=sum((off_original>ts_on(f)-4.5).*(off_original<ts_on(f)));%(off_original>ts_on(f)-4.5)减4.5需要手动调整此处ITI为4.5
end
hit(1)=sum(min(off_original)>min(x_on));
rows_start_f=(hit==0);
rows_start_final=(hit>0);
trial_start_final=ts_on(rows_start_final);
%根据deltaY巨变，不连续及ITI预测的最终led on的点
%video 233：切掉30s后余下743个trial

figure(figure_start);
start_exclused_by_off=scatter(ts_on(rows_start_f),ones(1,length(x_on(rows_start_f)))*(max(T.Y)+3),'k');
curve_trial_start_final=scatter(ts_on(rows_start_final),y_on(rows_start_final),'r');
%根据ITI滤除的点
%需要根据figure(figure_start)看一下led on和off之间的间距，ITI需要根据实际情况调整（看ITI滤除点是否正确）
legend([start_by_delta,exclusion,original,on_predicted,curve_led_off,start_exclused_by_off,curve_trial_start_final],{'start time predicted by delta Y','exclused time by sequence','brightness detected by ImageJ','predicted start time','off time predicted by 5 delta Y','led on exclused by ITI','led on final'});

title('led on predicted by movie analysis');
xlabel('time(s)');
set(gca,'FontSize',14);

%将led on换算为trial start（减去shift常数：根据水杆伸出计算而来）
shift_led=53.209300982091460;
%单位为ms，实际在画面中只相差1-2帧
shift_lickpredict=0;
trial_start_final=trial_start_final-shift_led/1000+shift_lickpredict;%对齐trial起始点
trial_start_OLED=trial_start_final;
save([path,filesep,'trialStartVideoByOLEDOn.mat'],'trial_start_OLED');
close(figure_start);
end

function [delta]=fdiff(X,n)
delta=[];
for i=1:(length(X)-n)
    delta_m=X(i+n)-X(i);
    delta=[delta,delta_m];
end
end
