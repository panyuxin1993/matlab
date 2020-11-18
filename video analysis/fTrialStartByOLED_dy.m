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
figure_start=figure;%����Ԥ�⻭��һ��ͼ�ϣ����۽����������
figure(figure_start);
original=plot(T.ts,T.Y,'b');%ԭʼLED����ͼ
%led onԤ��
led_on=(diff(T.Y)>6);
led_on=[false;led_on];
no_led_on=sum(led_on);
x_on=T.X(led_on);%�������Ⱦ��Ԥ��led onʱ��
y_on=T.Y(led_on);
ts_on=T.ts(led_on);
figure(figure_start);
hold on
start_by_delta=scatter(ts_on,y_on,'b');%��ͼ,�������Ⱦ��Ԥ��led onʱ��
%����һ��ɾ��������Ԥ���
led_continue=(diff(x_on)==1);%������Ԥ��㣬��Ҫ��ɾ��
led_continue=[led_continue;false];
x_continue=x_on(led_continue);
figure(figure_start);
hold on
exclusion=scatter(ts_on(led_continue),ones(length(x_continue),1)*(max(T.Y)+5),'k');%��ͼ��ȥ���ĵ㣨������Ԥ��㣩 
ts_on=ts_on(~led_continue);%��Ԥ�����ȥ�������� %ts_on=setdiff(ts_on,ts_on(led_continue));
y_on=y_on(~led_continue);
x_on=x_on(~led_continue);
%����pyx214_20190802��led���˸���delta>6�����������ֵ�ԭ���ó���Ԥ���668��Led����
led_after_continuity_exclusion=false(1,length(T.X));
for w=1:length(x_on)
    led_after_continuity_exclusion(T.X==x_on(w))=true;
end
% no_led=sum(led_after_continuity_exclusion);
% led_after_continuity_exclusion2=setdiff(x_on,x_continue);
% temp=T.X(led_after_continuity_exclusion)-led_after_continuity_exclusion2;
figure(figure_start);
hold on
on_predicted=scatter(T.ts(led_after_continuity_exclusion),T.Y(led_after_continuity_exclusion),'g');%�����ȱ仯����ֱ�ӻ���Ԥ���

%���ƶ���ITIʱ����Ժ㶨����������������ĵ�ҲҪ��ɾ��
% led_interval=diff(ts_on);%led on��interval�����̵Ŀ��ܾ����쳣�㣬������Ϊ����ITI��ֵ�Ĳο�
% figure_interval=figure;
% histogram(led_interval);
% xlabel('time(s)');
% title('led on interval');
%Ѱ�����ʱ��
delta_off=fdiff(T.Y,5);%����led���������������ʹ��5��Y�Ĳ�ֵԤ��off��
led_off=delta_off<-8;%this criterion is mannually set
add=false(5,1);
led_off=[add',led_off];%��Ԥ��off��ͳ�ơ�����Ҫ����Ԥ����Ƿ�׼ȷ��������쳣����Խ�һ��������������dleta_off<-8)���ٽ�����һ����
% figure_off=figure;
% histogram(delta_off);hold on;
% title('5 Y delta');
figure(figure_start);%��offԤ����led onԤ��㻭��һ��ͼ�����ڼ���ITI��ÿ����Ƶ���ڲ��
curve_led_off=scatter(T.ts(led_off),T.Y(led_off),'c');

off_original=T.ts(led_off);
%����Ԥ��led on��ITI��Χ���Ƿ����off������������˵��Ԥ�����
hit=zeros(1,length(ts_on));
for f=1:length(ts_on)
    hit(f)=sum((off_original>ts_on(f)-4.5).*(off_original<ts_on(f)));%(off_original>ts_on(f)-4.5)��4.5��Ҫ�ֶ������˴�ITIΪ4.5
end
hit(1)=sum(min(off_original)>min(x_on));
rows_start_f=(hit==0);
rows_start_final=(hit>0);
trial_start_final=ts_on(rows_start_final);
%����deltaY�ޱ䣬��������ITIԤ�������led on�ĵ�
%video 233���е�30s������743��trial

figure(figure_start);
start_exclused_by_off=scatter(ts_on(rows_start_f),ones(1,length(x_on(rows_start_f)))*(max(T.Y)+3),'k');
curve_trial_start_final=scatter(ts_on(rows_start_final),y_on(rows_start_final),'r');
%����ITI�˳��ĵ�
%��Ҫ����figure(figure_start)��һ��led on��off֮��ļ�࣬ITI��Ҫ����ʵ�������������ITI�˳����Ƿ���ȷ��
legend([start_by_delta,exclusion,original,on_predicted,curve_led_off,start_exclused_by_off,curve_trial_start_final],{'start time predicted by delta Y','exclused time by sequence','brightness detected by ImageJ','predicted start time','off time predicted by 5 delta Y','led on exclused by ITI','led on final'});

title('led on predicted by movie analysis');
xlabel('time(s)');
set(gca,'FontSize',14);

%��led on����Ϊtrial start����ȥshift����������ˮ��������������
shift_led=53.209300982091460;
%��λΪms��ʵ���ڻ�����ֻ���1-2֡
shift_lickpredict=0;
trial_start_final=trial_start_final-shift_led/1000+shift_lickpredict;%����trial��ʼ��
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
