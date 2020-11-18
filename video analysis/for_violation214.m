%����led on��trial start֮��Ĳ�ֵ�ҵ�trial start�ĵ�
dbstop if error
close all
clear

path='C:\Users\dy\Desktop\2019 ion\ʵ����ѧϰ\xulab\���ݷ���\��Ƶ����';
led_file=[path,'\214\pyx214_led.csv'];
file_trace=[path,'\214\2019-08-02-133157-FP-pyx214DeepCut_resnet50_tracking-lickport-tongueAug24shuffle1_100000.csv'];
behdata=[path,'\214\2019_08_02_pyx214-FP'];
%ITIʱ���Ϊ120������Ԥ��trial_start_finalΪ647

% led_file=[path,'\233\cropped_30scut.csv'];%2019-07-25-163638vio-pyx233.aviԴ��Ƶ
% file_trace=[path,'\233\2019-07-25-163638vio-pyx233-cut30s-deepcut.csv'];
% behdata=[path,'\233\pyx233_20190725'];
%��Ҫ�޸ĵĲ�����ITIʱ���Ϊ60������Ԥ��trial_start_finalΪ743��30scut��



%��������
T=readtable(led_file);
T_trace=readtable(file_trace,'HeaderLines',2);
load(behdata);



%behdata��������⵽��lick����trial������������ʵ������л������ͼ
figure_lick=figure;
subplot(3,2,1);
height=1;

licktime_arduino=cell(3,length(SessionResults));
for i=1:length(SessionResults)
    left=strsplit(SessionResults{i}.Action_lickTimeLeft,'|');
    left(:,1)=[];
    lickleft=str2double(left)/1000;
    for d=1:length(lickleft)
        dot_left=plot([lickleft(d),lickleft(d)],[i*height,(i+1)*height],'b','LineWidth',2);
        hold on
    end
    right=strsplit(SessionResults{i}.Action_lickTimeRight,'|');
    right(:,1)=[];
    lickright=str2double(right)/1000;
    for d=1:length(lickright)
        dot_right=plot([lickright(d),lickright(d)],[i*height,(i+1)*height],'r','LineWidth',2);
        hold on
    end
    licktime_arduino{1,i}=lickright;
    licktime_arduino{2,i}=lickleft;
    licktime_arduino{3,i}=[lickleft,lickright];
end
legend([dot_left,dot_right],{'left lick','right lick'});
xlabel('time(s)');
title('lick detected by arduino');
set(gca,'Xlim',[0,12]);
set(gca,'FontSize',14);
%%
%led onԤ��
figure_start=figure;
figure(figure_start);
original=plot(T.X,T.Y,'b');
%ԭʼLED����ͼ

led_delta=(diff(T.Y)>6);
led_delta=[false;led_delta];
no_led_delta=sum(led_delta);
x_delta=T.X(led_delta);
y_delta=T.Y(led_delta);
%�������Ⱦ��Ԥ��led onʱ��
figure(figure_start);
hold on
start_by_delta=scatter(x_delta,y_delta,'r');
%��ͼ

led_excluse=(diff(x_delta)==1);
led_excluse=[led_excluse;false];
x_excluse=x_delta(led_excluse);
%������Ԥ���
figure(figure_start);
hold on
exclusion=scatter(x_excluse,ones(length(x_excluse),1)*(max(T.Y)+5),'k');
%��ͼ��ȥ���ĵ�

for i=1:length(x_excluse)
x_delta(x_delta==x_excluse(i))=[];
end
%��Ԥ�����ȥ��������
%���˸���delta>6�����������ֵ�ԭ���ó���Ԥ���668��Led����

led_after_sequence_exclusion=false(1,length(T.X));
for w=1:length(x_delta)
    led_after_sequence_exclusion(T.X==x_delta(w))=true;
end
no_led=sum(led_after_sequence_exclusion);
figure(figure_start);
hold on
predicted=scatter(T.X(led_after_sequence_exclusion),T.Y(led_after_sequence_exclusion),'g');
%�����ȱ仯����ֱ�ӻ���Ԥ���
legend([start_by_delta,exclusion,original,predicted],{'start time predicted by delta Y','exclused time by sequence','brightness detected by ImageJ','predicted start time'});


delta_off=fdiff(T.Y,5);
%����led���������������ʹ��5��Y�Ĳ�ֵԤ��off��
figure_off=figure;
subplot(2,1,1);
histogram(delta_off)
title('5 Y delta');
hold on
rows_off=delta_off<-8;
subplot(2,1,2);
plot(T.X,T.Y,'b');
hold on
add=false(5,1);
rows_off=[add',rows_off];
off_predict=scatter(T.X(rows_off),T.Y(rows_off),'r');
title('off time predicted by 5 Y delta');
%��Ԥ��off��ͳ�ơ�����Ҫ����Ԥ����Ƿ�׼ȷ��������쳣����Խ�һ��������������dleta_off<-8)���ٽ�����һ����
figure(figure_start)
hold on
led_off=scatter(T.X(rows_off),T.Y(rows_off),'m');
legend([led_off],['led off predicted by 5 Y delta']);
%��offԤ����led onԤ��㻭��һ��ͼ�����ڼ���ITI��ÿ����Ƶ���ڲ��
off_original=T.X(rows_off);
led_interval=diff(x_delta);
figure_interval=figure;
histogram(led_interval);
title('led on interval');
%led on��interval�����̵Ŀ��ܾ����쳣�㣬������Ϊ����ITI��ֵ�Ĳο�

hit=zeros(1,length(x_delta));
for f=1:length(x_delta)
%     hit(f)=sum((off_original>x_delta(f)-60).*(off_original<x_delta(f)));
    hit(f)=sum((off_original>x_delta(f)-120).*(off_original<x_delta(f)));
end
%����Ԥ��led on��ITI��Χ���Ƿ����off������������˵��Ԥ�����
%(off_original>x_delta(f)-120)��120��Ҫ�ֶ�����

hit(1)=sum(min(off_original)>min(x_delta));
rows_start_f=(hit==0);
figure(figure_start);
hold on
start_exclused_by_off=scatter(x_delta(rows_start_f),ones(1,length(x_delta(rows_start_f)))*(max(T.Y)+3),'b');
%����ITI�˳��ĵ�
%��Ҫ����figure(figure_start)��һ��led on��off֮��ļ�࣬ITI��Ҫ����ʵ�������������ITI�˳����Ƿ���ȷ��
legend([start_by_delta,exclusion,original,off_predict,start_exclused_by_off],{'led on predicted by delta Y','exclused by sequence','brightness detected by ImageJ','off time predicted by 5 delta Y','led on exclused by ITI'});
title('led on predicted by movie analysis');
set(gca,'FontSize',14);

rows_start_final=(hit>0);
trial_start_final=x_delta(rows_start_final);
%����deltaY�ޱ䣬��������ITIԤ�������led on�ĵ�
%video 233���е�30s������743��trial


%��led on����Ϊtrial start����ȥshift����������ˮ��������������
shift_led=53.209300982091460;
%��λΪms��ʵ���ڻ�����ֻ���1-2֡
fr=24;
trial_start_final=trial_start_final/fr;
shift_lickpredict=0;
%��led on��֡������Ϊ��Ƶ�ڵ�ʱ����Ϣ����λΪs
trial_start_final=trial_start_final-shift_led/1000+shift_lickpredict;
%����trial��ʼ��
%%
%��ͷ��ˮ�˾���Ԥ��lick
% likelihood_portR=min(T_trace.likelihood_3);
% likelihood_portL=min(T_trace.likelihood_4);
% 
% % ����port�ɿ���ռ��
% rows_right=T_trace.likelihood_3>0.95;
% q=sum(rows_right);
% reliable_r=q/length(rows_right);
% 
% rows_left=T_trace.likelihood_4>0.95;
% p=sum(rows_left);
% reliable_r=p/length(rows_left);
% 
% %����port���ڹ̶�λ�ã���ѡ��Ԥ��λ�÷ֲ����ĵ���Ϊport������
% figure_reliableport=figure;
% scatter(T_trace.x_3,T_trace.y_3,'b');
% hold on
% x0_r=fmode(T_trace.x_3);
% y0_r=fmode(T_trace.y_3);
% scatter(x0_r,y0_r,'r');

% scatter(T_trace.x_4,T_trace.y_4,'y');
% hold on
% x0_l=fmode(T_trace.x_4);
% y0_l=fmode(T_trace.y_4);
% scatter(x0_l,y0_l,'g');
% title('location of lickports');
% set(gca,'FontSize',14);

%��ȡ��׽����tongue�Ŀɿ�����
% figure_tongue_to_lickports=figure;
% T_trace.distance_r=sqrt((T_trace.x-x0_r).^2+(T_trace.y-y0_r).^2);
% T_trace.distance_l=sqrt((T_trace.x-x0_l).^2+(T_trace.y-y0_l).^2);
% subplot(3,1,1);
% plot(1:length(T_trace.x),T_trace.distance_r,'r');
% hold on
% plot(1:length(T_trace.x),T_trace.distance_l+max(T_trace.distance_r),'b');
% title('distance of tongue to lickports');
% set(gca,'FontSize',14);

times_as_threshold=2;

% subplot(3,1,2);
% histogram(T_trace.distance_r);
% threshold_distance_r=mean(T_trace.distance_r(rows_right))-std(T_trace.distance_r(rows_right))*times_as_threshold;%y�᷽�򻮷�threshold��������Χ����Ϊ�ƶ���
% ylim=get(gca,'Ylim');
% hold on
% plot([threshold_distance_r,threshold_distance_r],[ylim(1),ylim(end)]);%����threshold��Χ
% text(threshold_distance_r,mean(ylim),num2str(threshold_distance_r));%�г���ֵ
% xlabel('distance of tongue to right lickport');
% set(gca,'FontSize',14);
% 
% 
% subplot(3,1,3);
% histogram(T_trace.distance_l);
% threshold_distance_l=mean(T_trace.distance_l(rows_left))-std(T_trace.distance_l(rows_left))*times_as_threshold;%y�᷽�򻮷�threshold��������Χ����Ϊ�ƶ���
% ylim_l=get(gca,'Ylim');
% hold on
% plot([threshold_distance_l,threshold_distance_l],[ylim_l(1),ylim_l(end)]);%����threshold��Χ
% text(threshold_distance_l,mean(ylim_l),num2str(threshold_distance_l));%�г���ֵ
% xlabel('distance of tongue to left lickport');
% set(gca,'FontSize',14);
%������ͷ��ˮ�����Ԥ��lick����ֵ
%%
% %distanceС����ֵ��,���ڴ�֮ǰdistance������ֵ�Ŀ���Ϊlick
% T_trace.distance_r(~rows_right)=max(T_trace.distance_r);
% right_distance_lick=(T_trace.distance_r<threshold_distance_r);
% right_lick_distance=(diff(right_distance_lick)>0);
% right_lick_distance=[0;right_lick_distance];
% fr_right_lick_distance=find(right_lick_distance);
% %�Ѳ��ɿ��ĵ��˳������������ĵ㣨���룬�����仯�������߼��жϣ���֣�find��
% T_trace.distance_l(~rows_left)=max(T_trace.distance_l);
% left_distance_lick=(T_trace.distance_l<threshold_distance_l);
% left_lick_distance=(diff(left_distance_lick)>0);
% left_lick_distance=[0;left_lick_distance];
% fr_left_lick_distance=find(left_lick_distance);
% 
% 
% time_right_lick_distance=fr_right_lick_distance/fr;
% time_left_lick_distance=fr_left_lick_distance/fr;
% figure(figure_lick)
% hold on
% subplot(3,2,2);
% licktime_by_distance=cell(3,length(trial_start_final));
temp_trial_start=[trial_start_final;length(T.X)/fr];
% %תΪʱ�䣻cell����
% for g=1:length(trial_start_final)
%     right_cycle=time_right_lick_distance(logical((time_right_lick_distance>temp_trial_start(g)).*(time_right_lick_distance<temp_trial_start(g+1))));
%     left_cycle=time_left_lick_distance(logical((time_left_lick_distance>temp_trial_start(g)).*(time_left_lick_distance<temp_trial_start(g+1))));
%     licktime_by_distance{1,g}=right_cycle-trial_start_final(g);
%     licktime_by_distance{2,g}=left_cycle-trial_start_final(g);
%     licktime_by_distance{3,g}=[right_cycle-trial_start_final(g);left_cycle-trial_start_final(g)];
%     %��ÿ��Ԥ��lick��Ӧ��ÿ��trial��ʱ����Ϣ����cell�������trial����ʼʱ��
%     for ri=1:length(right_cycle)
%         do_ri=plot([right_cycle(ri)-temp_trial_start(g),right_cycle(ri)-temp_trial_start(g)],[g*height,(g+1)*height],'r','LineWidth',2);
%         hold on
%     end
%   %����
%     
%      for le=1:length(left_cycle)
%         do_le=plot([left_cycle(le)-temp_trial_start(g),left_cycle(le)-temp_trial_start(g)],[g*height,(g+1)*height],'b','LineWidth',2);
%         hold on
%     end
% end
% legend([do_le,do_ri],{'left lick','right lick'});
% xlabel('time(s)');
% title('lick predicted by distance of tongue to lickports');
% set(gca,'FontSize',14);
%%        
figure_y=figure;
rows=T_trace.likelihood>0.95;
histogram(T_trace.y(rows));hold on;%hold on����������ͼƬ
xlabel('Coordinates');
title('y distribution');
set(gca,'FontSize',14);
thresholdy=mean(T_trace.y(rows))+std(T_trace.y(rows))*times_as_threshold;%y�᷽�򻮷�threshold��������Χ����Ϊ�ƶ���
ylim=get(gca,'Ylim');%������������
plot([thresholdy,thresholdy],[ylim(1),ylim(end)]);%����threshold��Χ
text(thresholdy,mean(ylim),num2str(thresholdy));%�г���ֵ
set(gca,'FontSize',14);

y0=fmode(T_trace.y);

T_trace.y(~rows)=y0;%filter the points that likelihood<0.95
t_lick=(T_trace.y>thresholdy);%using histogram to find about 3.3 may be threshold for licking
ind_lick3=(diff(t_lick)>0);
ind_lick3=[0;ind_lick3];
p3=find(ind_lick3);
licktime_by_y=p3/fr;
figure(figure_lick); 
hold on
subplot(3,2,3);

licktime_by_y_trial=cell(1,length(trial_start_final));
for g=1:length(trial_start_final)
    y_cycle=licktime_by_y(logical((licktime_by_y>temp_trial_start(g)).*(licktime_by_y<temp_trial_start(g+1))));
    licktime_by_y_trial{1,g}=y_cycle-trial_start_final(g);    
    
     for le=1:length(y_cycle)
        dot_y=plot([y_cycle(le)-temp_trial_start(g),y_cycle(le)-temp_trial_start(g)],[g*height,(g+1)*height],'k','LineWidth',2);
        hold on
    end
end
xlabel('time(s)');
title('lick predicted by y of tongue');
set(gca,'FontSize',14);

%%
x0=fmode(T_trace.x);%��T��x�н��м��㣬ȷ��Ŀ���ƽ��λ�㣨�ֲ���ࣩ������������������
y0=fmode(T_trace.y);%use mode as origin
T_trace.distance=sqrt((T_trace.x-x0).^2+(T_trace.y-y0).^2);%����λ�ƾ���

figure_distance=figure;
histogram(T_trace.distance(rows));hold on;
xlabel('Distance');
title('distance distribution');
threshold=mean(T_trace.distance(rows))+std(T_trace.distance(rows));
ylim=get(gca,'Ylim');
plot([threshold,threshold],[ylim(1),ylim(end)]);%����һ����
text(threshold,mean(ylim),num2str(threshold));%����text�����꼰��ֵ
set(gca,'FontSize',14);

T_trace.distance(~rows)=0;%filter the points that likelihood<0.95
t_lick=(T_trace.distance>threshold);%using histogram to find about 3.3 may be threshold for licking
ind_lick2=(diff(t_lick)>0);
ind_lick2=[0;ind_lick2];
p2=find(ind_lick2);

licktime_by_distance_tongue=p2/fr;
figure(figure_lick);
hold on
subplot(3,2,4);

licktime_by_distance_trial=cell(1,length(trial_start_final));
for g=1:length(trial_start_final)
    tongue_cycle=licktime_by_distance_tongue(logical((licktime_by_distance_tongue>temp_trial_start(g)).*(licktime_by_distance_tongue<temp_trial_start(g+1))));
    licktime_by_distance_trial{1,g}=tongue_cycle-trial_start_final(g);    
    
     for le=1:length(tongue_cycle)
        dot_tongue=plot([tongue_cycle(le)-temp_trial_start(g),tongue_cycle(le)-temp_trial_start(g)],[g*height,(g+1)*height],'k','LineWidth',2);
        hold on
    end
end
xlabel('time(s)');
title('lick predicted by shift distance of tongue');
set(gca,'FontSize',14);
%%
%����Ҷ�任�˳���ͷ��y�����㣬����lick rate
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
wp = [6 8 ] / (fr/2);%ͨ����ֹƵ��
ws = [5 9 ] / (fr/2);  %�����ֹƵ��
alpha_p = 3; %ͨ���������˥��Ϊ  db
alpha_s = 10;%���������С˥��Ϊ  db
%��ȡ�����ͽ�ֹƵ��
[ N3, wn ] = buttord( wp , ws , alpha_p , alpha_s);
%���ת�ƺ���ϵ��
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


figure(figure_lick);
hold on
subplot(3,2,5);

licktime_by_y_filter_trial=cell(1,length(trial_start_final));
for g=1:length(trial_start_final)
    y_filter_cycle=licktime_by_y_filter(logical((licktime_by_y_filter>temp_trial_start(g)).*(licktime_by_y_filter<temp_trial_start(g+1))));
    licktime_by_y_filter_trial{1,g}=y_filter_cycle-trial_start_final(g);    
    
     for le=1:length(y_filter_cycle)
        dot_y_filter=plot([y_filter_cycle(le)-temp_trial_start(g),y_filter_cycle(le)-temp_trial_start(g)],[g*height,(g+1)*height],'k','LineWidth',2);
        hold on
    end
end
xlabel('time(s)');
title('lick predicted by y filtered of tongue');
set(gca,'FontSize',14);
%%
%����ͬ�ֶ���Ĳ�ͬԤ��Ч��

%evaluate prediction based on video
% 1)compare recall and precision of different methods ��ʱ�䴰��Χ��Ѱ��hit��lick
setwin_ind=2;
figEvaluate=figure;
subplot(2,2,1);
win=0.1:0.1:1;

% % recall_distance=zeros(1,length(win));
% recall_y=zeros(1,length(win));
% recall_filter=zeros(1,length(win));
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

%%
subplot(2,2,2);
% bar(1,recall_distance_r(setwin_ind),'FaceColor','b');hold on;
% bar(2,recall_distance_l(setwin_ind),'FaceColor','r');
bar(1,recall_y_filter(setwin_ind),'FaceColor','k');hold on;
bar(2,recall_y(setwin_ind),'FaceColor','m');
bar(3,recall_distance_trial(setwin_ind),'FaceColor','g');

set(gca,'XTick',[]);
set(gca,'FontSize',14);

%%
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


%%
subplot(2,2,4);
% bar(1,recall_distance_r(setwin_ind),'FaceColor','b');hold on;
% bar(2,recall_distance_l(setwin_ind),'FaceColor','r');
bar(1,recall_y_filter(setwin_ind),'FaceColor','k');hold on;
bar(2,recall_y(setwin_ind),'FaceColor','m');
bar(3,recall_distance_trial(setwin_ind),'FaceColor','g');
set(gca,'XTick',[]);
set(gca,'FontSize',14);
%%
%??????
dbstop if error
fig_histo_lickrate=figure;
timeMin=1000;
binStep=50;
trialLength=3000;
binSize=500;


[ ~,lickRate_distance_tongue ] = fLickRateAligned(licktime_by_distance_trial(1,1:size(licktime_arduino,2)),licktime_arduino(3,:),binSize,trialLength,timeMin, binStep);
% [ ~,lickRate_distance_port ] = fLickRateAligned( licktime_by_distance(3,:), licktime_arduino(3,4:end),binSize,trialLength,timeMin, binStep);
[ ~,lickRate_y ] = fLickRateAligned( licktime_by_y_trial(1,1:size(licktime_arduino,2)),licktime_arduino(3,:),binSize,trialLength,timeMin, binStep);
[ ~,lickRate_y_filter ] = fLickRateAligned( licktime_by_y_filter_trial(1,1:size(licktime_arduino,2)),licktime_arduino(3,:),binSize,trialLength,timeMin, binStep);
ts=-timeMin:binStep:trialLength-timeMin-binSize;
plot(ts',lickRate_y);
hold on;
plot(ts',lickRate_y_filter);
plot(ts',lickRate_distance_tongue);
% plot(ts',lickRate_distance_port);

xlabel('Time from true lick(ms)');
ylabel('Lick rate(/s)');
legend('y of tongue','y filtered of tongue','distance of tongue');
set(gca,'FontSize',14,'FontName','Arial');
%%
%evaluate prediction based on video
% 2) trial type prediction

n_trial=size(licktime_arduino,2);


time_delay_off=zeros(n_trial,1);
time_stim_onset=zeros(n_trial,1);

time_stim_onset=cellfun(@(x) double(x.Time_stimOnset)/1000,SessionResults);
% time_stim_onset(1:3)=[];
time_delay_off=cellfun(@(x) (double(x.Time_stimOffset)+double(x.Delay_duration))/1000,SessionResults);
% time_delay_off(1:3)=[];


vio_y=fvio(licktime_by_y_trial,time_stim_onset,time_delay_off,n_trial,1);
vio_distance_tongue=fvio(licktime_by_distance_trial,time_stim_onset,time_delay_off,n_trial,1);
% vio_distance_port=fvio(licktime_by_distance,time_stim_onset,time_delay_off,n_trial,3);
vio_y_filter=fvio(licktime_by_y_filter_trial,time_stim_onset,time_delay_off,n_trial,1);


vio_arduino=cellfun(@(x) x.Action_choice==3,SessionResults);
% vio_arduino(1:3)=[];

nh=1:n_trial;
figVio=figure;
subplot(2,2,1);
% dot_vio_arduino=scatter(nh(vio_arduino),ones(length(nh(vio_arduino)),1));
% hold on

legendstr={'arduino','y of tongue','y filtered of tongue','distance of tongue'};

dot_vio=fviodot(legendstr,vio_arduino,vio_y,vio_y_filter,vio_distance_tongue);
%����Ԥ�ⷽ����vio

% evio_y_t=sum(vio_arduino'.*vio_y==1)/sum(vio_arduino);%������
% evio_y_f=sum((vio_arduino'+vio_y)==0)/sum(vio_arduino==0);%������
% evio_y_filter_t=sum(vio_arduino'.*vio_y_filter==1)/sum(vio_arduino);
% evio_y_filter_f=sum((vio_arduino'+vio_y_filter)==0)/sum(vio_arduino==0);
% evio_v_distance_tongue_t=sum(vio_arduino'.*vio_distance_tongue==1)/sum(vio_arduino);
% evio_v_distance_tongue_f=sum((vio_arduino'+vio_distance_tongue)==0)/sum(vio_arduino==0);
% figure(figVio);
% subplot(2,2,2);

% bar(1,evio_y_t,'FaceColor','r');hold on;
% bar(2,evio_y_filter_t,'FaceColor','k');
% bar(3,evio_v_distance_tongue_t,'FaceColor','m');
% bar(4,evio_v_distance_port_t,'FaceColor','g');
% legend('y of tongue','y filtered of tongue','distance of tongue','distance to port of tongue');
% title('violation predicted correctly');
% set(gca,'XTick',[]);
% set(gca,'FontSize',14);
% 
% subplot(2,2,3);
% bar(1,evio_y_f,'FaceColor','r');hold on;
% bar(2,evio_y_filter_f,'FaceColor','k');
% bar(3,evio_v_distance_tongue_f,'FaceColor','m');
% bar(4,evio_v_distance_port_f,'FaceColor','g');
% legend('y of tongue','y filtered of tongue','distance of tongue','distance to port of tongue');
% title('no violation predicted correctly');
% set(gca,'XTick',[]);
% set(gca,'FontSize',14);
% 
% subplot(2,2,4);
% scatter(1,sum(vio_arduino)/n_trial,'b');
% hold on
% scatter(2,sum(vio_y)/n_trial,'r');
% scatter(3,sum(vio_y_filter)/n_trial,'k');
% scatter(4,sum(vio_distance_tongue)/n_trial,'m');
% scatter(5,sum(vio_distance_port)/n_trial,'g');
% legend('arduino','y of tongue','y filtered of tongue','distance of tongue','distance to port of tongue');
% title('violation rate');
% set(gca,'XTick',[]);
% set(gca,'FontSize',14);
%%
%����ͬ��Ԥ��Ķ���Ч��






%%
function [f0]=fmode(f)%�ҵ����ķֲ�ֵ���䣿
[N,edges,bin] = histcounts(f,100);%��f�ֳ�100����С��һ��bin
f0 = edges(N == max(N));%�ں����Ԫ�ص�bin��edge
if length(f0)==2
    f0=mean(f0);
end
end


function [delta]=fdiff(X,n)
delta=[];
for i=1:(length(X)-n)
    delta_m=X(i+n)-X(i);
    delta=[delta,delta_m];
end
end



function [ratio,hit]= frecall(data,template,window)
%input: data and template are both vectors,use the second as
%baseline,window means the time window that tolerate time shifts,[a,b]
%output: ratio of elements in template, for which data recalled
% upgrade 2020.1.12 to compatible with data that organized by trials(each
% trial a cell)
if (~iscell(data))&&(~iscell(template))
    hit=zeros(1,length(template));
    for i=1:length(template)
        hit(i)=sum((data<template(i)+window(2)).*(data>template(i)-window(1)));
    end
    ratio=1-sum(hit==0)/length(hit);
elseif iscell(data)&&iscell(template)&&(length(data)==length(template))%here data and template should both be cells with same length
    hit=cell(1,length(data)); 
    for n_trial=1:length(data)
        if isnan(template{1,n_trial})
            hit{1,n_trial}=[];
        else
            hit{1,n_trial}=zeros(1,length(template{1,n_trial}));
            for i=1:length(template{1,n_trial})
                hit{1,n_trial}(i)=sum((data{1,n_trial}<template{1,n_trial}(i)+window(2)).*(data{1,n_trial}>template{1,n_trial}(i)-window(1)));
            end
        end
    end
    n_template=cellfun(@length , hit);
    n_template=sum(n_template);
    n_hit=cellfun(@(x) sum(x~=0), hit);
    n_hit=sum(n_hit);
    ratio=n_hit/n_template;
else
    warning('Two input variables are not same type');
end
end

% function [ lickMat,lickRate ] = fLickRateAligned( lickTime, base,binSize,trialLength,timeMin, binStep)
% sampleNum=ceil((trialLength-binSize+binStep)/binStep)+1;%��λms
% lickMat=zeros(length(base),trialLength);%��һάtrials,�ڶ�ά����trial����
% lickRate=zeros(sampleNum,1);%��һά����trial����/��������
% for i=1:length(base)
%     lickT=(lickTime-base(i))*1000;%��λs
%     for j=1:length(lickT)
%         if ~isnan(lickT(j)) && lickT(j)>-timeMin+1 && lickT(j)<trialLength-timeMin-1
%             lickMat(i,round(lickT(j))+timeMin)=1;
%         end
%     end
% end
% for j=1:sampleNum-ceil(binSize/binStep)
%     lickRate(j+ceil(binSize/binStep))=nansum(nansum(lickMat(:,(binStep*(j+ceil(binSize/binStep)-1)-binSize+1):(binStep*(j+ceil(binSize/binStep)-1)))))*1000/binSize/length(base);
% end
% 
% end

function [vio]=fvio(licktime,stimulus_onset,delay_off,n_trial,cell_lines)
    vio=false(n_trial,1);
    for i=1:n_trial
        temp=logical(licktime{cell_lines,i}>=stimulus_onset(i)).*(licktime{cell_lines,i}<=delay_off(i));
        if ~isempty(temp)
            vio(i)=(sum(temp)>0);
        end
    end
end


function [dot_vio]=fviodot(legendstr,varargin)
    for i=1:length(varargin)
        vio=varargin{i};
        nh=1:length(vio);
        dot_vio(i)=scatter(nh(vio),ones(length(nh(vio)),1)*i,500);
        hold on
    end
    legend(dot_vio,legendstr);
end


