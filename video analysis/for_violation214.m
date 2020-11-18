%根据led on和trial start之间的差值找到trial start的点
dbstop if error
close all
clear

path='C:\Users\dy\Desktop\2019 ion\实验室学习\xulab\数据分析\视频对齐';
led_file=[path,'\214\pyx214_led.csv'];
file_trace=[path,'\214\2019-08-02-133157-FP-pyx214DeepCut_resnet50_tracking-lickport-tongueAug24shuffle1_100000.csv'];
behdata=[path,'\214\2019_08_02_pyx214-FP'];
%ITI时间差为120，最终预测trial_start_final为647

% led_file=[path,'\233\cropped_30scut.csv'];%2019-07-25-163638vio-pyx233.avi源视频
% file_trace=[path,'\233\2019-07-25-163638vio-pyx233-cut30s-deepcut.csv'];
% behdata=[path,'\233\pyx233_20190725'];
%需要修改的参数：ITI时间差为60，最终预测trial_start_final为743（30scut）



%导入数据
T=readtable(led_file);
T_trace=readtable(file_trace,'HeaderLines',2);
load(behdata);



%behdata处理：将检测到的lick按照trial画出来，类似实验过程中绘出来的图
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
%led on预测
figure_start=figure;
figure(figure_start);
original=plot(T.X,T.Y,'b');
%原始LED亮度图

led_delta=(diff(T.Y)>6);
led_delta=[false;led_delta];
no_led_delta=sum(led_delta);
x_delta=T.X(led_delta);
y_delta=T.Y(led_delta);
%根据亮度剧变预测led on时间
figure(figure_start);
hold on
start_by_delta=scatter(x_delta,y_delta,'r');
%作图

led_excluse=(diff(x_delta)==1);
led_excluse=[led_excluse;false];
x_excluse=x_delta(led_excluse);
%连续的预测点
figure(figure_start);
hold on
exclusion=scatter(x_excluse,ones(length(x_excluse),1)*(max(T.Y)+5),'k');
%作图：去除的点

for i=1:length(x_excluse)
x_delta(x_delta==x_excluse(i))=[];
end
%从预测点中去除连续点
%至此根据delta>6及不连续出现的原则获得初步预测的668个Led亮点

led_after_sequence_exclusion=false(1,length(T.X));
for w=1:length(x_delta)
    led_after_sequence_exclusion(T.X==x_delta(w))=true;
end
no_led=sum(led_after_sequence_exclusion);
figure(figure_start);
hold on
predicted=scatter(T.X(led_after_sequence_exclusion),T.Y(led_after_sequence_exclusion),'g');
%在亮度变化线上直接画出预测点
legend([start_by_delta,exclusion,original,predicted],{'start time predicted by delta Y','exclused time by sequence','brightness detected by ImageJ','predicted start time'});


delta_off=fdiff(T.Y,5);
%由于led灭比亮更缓慢，故使用5个Y的差值预测off点
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
%对预测off的统计——主要看下预测的是否准确（如果有异常点可以进一步调整条件（如dleta_off<-8)后再进行下一步）
figure(figure_start)
hold on
led_off=scatter(T.X(rows_off),T.Y(rows_off),'m');
legend([led_off],['led off predicted by 5 Y delta']);
%将off预测点和led on预测点画在一张图：便于计算ITI（每个视频存在差别）
off_original=T.X(rows_off);
led_interval=diff(x_delta);
figure_interval=figure;
histogram(led_interval);
title('led on interval');
%led on的interval，过短的可能就是异常点，可以作为计算ITI阈值的参考

hit=zeros(1,length(x_delta));
for f=1:length(x_delta)
%     hit(f)=sum((off_original>x_delta(f)-60).*(off_original<x_delta(f)));
    hit(f)=sum((off_original>x_delta(f)-120).*(off_original<x_delta(f)));
end
%计算预测led on在ITI范围内是否存在off——不存在则说明预测错误
%(off_original>x_delta(f)-120)减120需要手动调整

hit(1)=sum(min(off_original)>min(x_delta));
rows_start_f=(hit==0);
figure(figure_start);
hold on
start_exclused_by_off=scatter(x_delta(rows_start_f),ones(1,length(x_delta(rows_start_f)))*(max(T.Y)+3),'b');
%根据ITI滤除的点
%需要根据figure(figure_start)看一下led on和off之间的间距，ITI需要根据实际情况调整（看ITI滤除点是否正确）
legend([start_by_delta,exclusion,original,off_predict,start_exclused_by_off],{'led on predicted by delta Y','exclused by sequence','brightness detected by ImageJ','off time predicted by 5 delta Y','led on exclused by ITI'});
title('led on predicted by movie analysis');
set(gca,'FontSize',14);

rows_start_final=(hit>0);
trial_start_final=x_delta(rows_start_final);
%根据deltaY巨变，不连续及ITI预测的最终led on的点
%video 233：切掉30s后余下743个trial


%将led on换算为trial start（减去shift常数：根据水杆伸出计算而来）
shift_led=53.209300982091460;
%单位为ms，实际在画面中只相差1-2帧
fr=24;
trial_start_final=trial_start_final/fr;
shift_lickpredict=0;
%将led on的帧数换算为视频内的时间信息，单位为s
trial_start_final=trial_start_final-shift_led/1000+shift_lickpredict;
%对齐trial起始点
%%
%舌头到水杆距离预测lick
% likelihood_portR=min(T_trace.likelihood_3);
% likelihood_portL=min(T_trace.likelihood_4);
% 
% % 左右port可靠点占比
% rows_right=T_trace.likelihood_3>0.95;
% q=sum(rows_right);
% reliable_r=q/length(rows_right);
% 
% rows_left=T_trace.likelihood_4>0.95;
% p=sum(rows_left);
% reliable_r=p/length(rows_left);
% 
% %由于port处于固定位置，故选择预测位置分布最多的点作为port的坐标
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

%提取捕捉到的tongue的可靠坐标
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
% threshold_distance_r=mean(T_trace.distance_r(rows_right))-std(T_trace.distance_r(rows_right))*times_as_threshold;%y轴方向划分threshold（超出范围才认为移动）
% ylim=get(gca,'Ylim');
% hold on
% plot([threshold_distance_r,threshold_distance_r],[ylim(1),ylim(end)]);%画出threshold范围
% text(threshold_distance_r,mean(ylim),num2str(threshold_distance_r));%列出数值
% xlabel('distance of tongue to right lickport');
% set(gca,'FontSize',14);
% 
% 
% subplot(3,1,3);
% histogram(T_trace.distance_l);
% threshold_distance_l=mean(T_trace.distance_l(rows_left))-std(T_trace.distance_l(rows_left))*times_as_threshold;%y轴方向划分threshold（超出范围才认为移动）
% ylim_l=get(gca,'Ylim');
% hold on
% plot([threshold_distance_l,threshold_distance_l],[ylim_l(1),ylim_l(end)]);%画出threshold范围
% text(threshold_distance_l,mean(ylim_l),num2str(threshold_distance_l));%列出数值
% xlabel('distance of tongue to left lickport');
% set(gca,'FontSize',14);
%计算舌头到水嘴距离预测lick的阈值
%%
% %distance小于阈值的,且在此之前distance大于阈值的可能为lick
% T_trace.distance_r(~rows_right)=max(T_trace.distance_r);
% right_distance_lick=(T_trace.distance_r<threshold_distance_r);
% right_lick_distance=(diff(right_distance_lick)>0);
% right_lick_distance=[0;right_lick_distance];
% fr_right_lick_distance=find(right_lick_distance);
% %把不可靠的点滤除，满足条件的点（距离，连续变化）——逻辑判断，差分，find；
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
% %转为时间；cell储存
% for g=1:length(trial_start_final)
%     right_cycle=time_right_lick_distance(logical((time_right_lick_distance>temp_trial_start(g)).*(time_right_lick_distance<temp_trial_start(g+1))));
%     left_cycle=time_left_lick_distance(logical((time_left_lick_distance>temp_trial_start(g)).*(time_left_lick_distance<temp_trial_start(g+1))));
%     licktime_by_distance{1,g}=right_cycle-trial_start_final(g);
%     licktime_by_distance{2,g}=left_cycle-trial_start_final(g);
%     licktime_by_distance{3,g}=[right_cycle-trial_start_final(g);left_cycle-trial_start_final(g)];
%     %将每个预测lick对应到每个trial（时间信息），cell储存相对trial的起始时间
%     for ri=1:length(right_cycle)
%         do_ri=plot([right_cycle(ri)-temp_trial_start(g),right_cycle(ri)-temp_trial_start(g)],[g*height,(g+1)*height],'r','LineWidth',2);
%         hold on
%     end
%   %画点
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
%评价同种对齐的不同预测效果

%evaluate prediction based on video
% 1)compare recall and precision of different methods 在时间窗范围内寻找hit的lick
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
%五种预测方法的vio

% evio_y_t=sum(vio_arduino'.*vio_y==1)/sum(vio_arduino);%真阳性
% evio_y_f=sum((vio_arduino'+vio_y)==0)/sum(vio_arduino==0);%真阴性
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
%评价同种预测的对齐效果






%%
function [f0]=fmode(f)%找到最大的分布值区间？
[N,edges,bin] = histcounts(f,100);%将f分成100个大小均一的bin
f0 = edges(N == max(N));%内含最多元素的bin的edge
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
% sampleNum=ceil((trialLength-binSize+binStep)/binStep)+1;%单位ms
% lickMat=zeros(length(base),trialLength);%第一维trials,第二维单个trial长度
% lickRate=zeros(sampleNum,1);%第一维代表trial长度/采样数量
% for i=1:length(base)
%     lickT=(lickTime-base(i))*1000;%单位s
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


