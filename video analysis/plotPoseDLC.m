dbstop if error
path='D:\xulab\behavior\video tracking\tracking-lickport-tongue-pyx-2019-06-20\videos\';
file_filter=[path,'2019-06-18-111206-retract-day1DeepCut_resnet50_tracking-lickport-tongueJun20shuffle1_210000filtered.csv'];
file_raw=[path,'2019-06-18-111206-retract-day1DeepCut_resnet50_tracking-lickport-tongueJun20shuffle1_210000.csv'];
file_test='D:\xulab\behavior\pyx233\2019-07-25-163638vio-pyx233DeepCut_resnet50_tracking-lickport-tongueJun20shuffle1_210000.csv';
file_predict=[path,'15s\trancatedDeepCut_resnet50_tracking-lickport-tongueJun20shuffle1_210000.csv'];
file_predict_compress=[path,'15s\smallDeepCut_resnet50_tracking-lickport-tongueJun20shuffle1_210000.csv'];
% behdata=[path,'pyx226_20190618.mat'];
behdata='D:\xulab\behavior\pyx233\pyx233_20190725.mat';
%data=csvread(file_raw,3,0);
% opts = detectImportOptions(file_raw)
%T.Properties.VariableNames = {'Gender' 'Age' 'Height' 'Weight'}

% T=readtable(file_raw,'HeaderLines',2);%skip two rows
T=readtable(file_test,'HeaderLines',2);
Tf=readtable(file_filter,'HeaderLines',2);%skip two rows
Tt=readtable(file_test,'HeaderLines',2);
Tp=readtable(file_predict,'HeaderLines',2);
Tpc=readtable(file_predict_compress,'HeaderLines',2);
%use ground true licking
load(behdata);
t0=strsplit(SessionResults{1}.Time_trialStart,'_');
t0=str2double(t0);
time_trial_start0=t0(4)*3600+t0(5)*60+t0(6)+t0(7)/1000;
time_shift=2;%video record 2s before 1st trial start
licktimeL=[];
licktimeR=[];
k=5811/5770;
time_trial_start=zeros(length(SessionResults),1);% time(s)
time_delay_off=zeros(length(SessionResults),1);
time_stim_onset=zeros(length(SessionResults),1);
for i=1:length(SessionResults)
    t=strsplit(SessionResults{i}.Time_trialStart,'_');
    t=str2double(t);
    time_trial_start(i)=t(4)*3600+t(5)*60+t(6)+t(7)/1000-time_trial_start0+time_shift;
    time_stim_onset(i)=time_trial_start(i)+SessionResults{i}.Time_stimOnset/1000;
    time_delay_off(i)=time_trial_start(i)+SessionResults{i}.Time_stimOffset/1000+SessionResults{i}.Delay_duration/1000;
    left=strsplit(SessionResults{i}.Action_lickTimeLeft,'|');
    left(:,1)=[];
    lickleft=str2double(left)/1000+time_trial_start(i);
    licktimeL=[licktimeL,lickleft];
    right=strsplit(SessionResults{i}.Action_lickTimeRight,'|');
    right(:,1)=[];
    lickright=str2double(right)/1000+time_trial_start(i);
    licktimeR=[licktimeR, lickright];
end
%calibarate time of lick port detected licking to video
licktimeL=k*licktimeL;
licktimeR=k*licktimeR;
time_trial_start=k*time_trial_start;
time_stim_onset=k*time_stim_onset;
time_delay_off=k*time_delay_off;

fr=24;%frame rate is 24Hz
times_as_threshold=2;
fr_t=24;
ts=(1:size(T,1))/fr;
ts_t=(1:size(Tt,1))/fr_t;
x0=fmode(T.x);
y0=fmode(T.y);%use mode as origin
T.distance=sqrt((T.x-x0).^2+(T.y-y0).^2);
rows=T.likelihood>0.95;%1:length(T.likelihood);%
rows2=1:size(T,1);
%plot histgram of tongue position and likelihood
fig_hist_position=figure;
subplot(2,2,1);
histogram(T.x(rows));
xlabel('Coordinates');
title('x distribution');
set(gca,'FontSize',14);
subplot(2,2,2);
histogram(T.y(rows));hold on;
xlabel('Coordinates');
title('y distribution');
set(gca,'FontSize',14);
thresholdy=mean(T.y(rows))+std(T.y(rows))*times_as_threshold;
ylim=get(gca,'Ylim');
plot([thresholdy,thresholdy],[ylim(1),ylim(end)]);
text(thresholdy,mean(ylim),num2str(thresholdy));
set(gca,'FontSize',14);
subplot(2,2,3);
histogram(T.likelihood(rows));
xlabel('likelihood');
set(gca,'FontSize',14);
subplot(2,2,4);
histogram(T.distance(rows));hold on;
xlabel('Distance');
title('distance distribution');
threshold=mean(T.distance(rows))+std(T.distance(rows))*times_as_threshold;
ylim=get(gca,'Ylim');
plot([threshold,threshold],[ylim(1),ylim(end)]);
text(threshold,mean(ylim),num2str(threshold));
set(gca,'FontSize',14);

figtrajectory=figure;
scatter(T.x(rows),T.y(rows),'k');
hold on;
scatter(T.x_3(rows),T.y_3(rows),'b');
scatter(T.x_4(rows),T.y_4(rows),'g');
scatter(x0,y0,'r');
legend('tongue','lickport-right','lickport-left','origin-tongue');
xlabel('x-tongue');
ylabel('y-tongue');
set(gca,'FontSize',14);

figpose=figure;
% %perform smooth, then set the threshold for licking
% T.distance(rows)=T.distance(rows)-smooth(T.distance(rows),fr+1)+mean(T.distance(rows));
% T.x(rows)=T.x(rows)-smooth(T.x(rows),fr+1)+mean(T.x(rows));
% T.y(rows)=T.y(rows)-smooth(T.y(rows),fr+1)+mean(T.y(rows));
% set(gcf,'Position',[100,100,500,400]);
% x_lim=[520,580];
% subplot(2,1,1);
y_adjust_d=0;%(mean(T.x(rows))+mean(T.y(rows)))/2;%平移distance曲线至x,y之间合适位置
plot(ts(rows),T.distance(rows)+y_adjust_d);hold on;
y_adjust_x=max(T.distance(rows))-min(T.distance(rows))+y_adjust_d;
plot(ts(rows),T.x(rows)-min(T.x(rows))+y_adjust_x);
y_adjust_y=max(T.x(rows))-min(T.x(rows))+y_adjust_x;
plot(ts(rows),T.y(rows)-min(T.y(rows))+y_adjust_y);
y_adjust_filter=max(T.y(rows))-min(T.y(rows))+y_adjust_y;
% plot(ts(rows),T.y_4(rows));
% plot(ts(rows),T.y_3(rows));
% h=legend('y-right-lick-port','y-left-lick-port','x-tongue','y-tongue','distance-tongue','AutoUpdate','off');
h=legend('x-tongue','y-tongue','distance-tongue','AutoUpdate','off');
set(h,'Box','off');
ylabel('Coordinates');
xlabel('Time(s)');
% set(gca,'XTickLabel',[]);
% set(gca,'Xlim',x_lim);
set(gca,'FontSize',14);
box off;
%%simultaneously see likelihood
% subplot(2,1,2);
% plot(ts,T.likelihood_4);hold on;
% plot(ts,T.likelihood_3);
% plot(ts,T.likelihood);
% % plot(ts_t,Tt.likelihood);
% % legend('likelihood','likelihood_t');
% ylabel('likelihood')
% xlabel('Time(s)');
% h=legend('right-lick-port','left-lick-port','tongue','AutoUpdate','off');
% set(h,'Box','off');
% % set(gca,'Xlim',x_lim);
% set(gca,'FontSize',14);
% box off;

%fft, try to see frequency distribution
figfft=figure;
S=T.y;
Sf=T.y(rows);
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
figure(figpose);
plot(ts,filter_bp_s+y_adjust_filter);
%set a threshold for licking- method-I
std_fbp_s=std(filter_bp_s);
times=2;
a=(filter_bp_s>times*std_fbp_s);
b=(diff(a)>0);%diff results will have less element,so add back
b=[0;b];
ind_lick=b;
% ind_lick=(filter_bp_s==3*std_fbp_s).*(diff(filter_bp_s)>0);%choose 3 std and increase time point as licking time point
p=find(ind_lick);
figure(figpose);
y_adjust_lick=y_adjust_filter+prctile(filter_bp_s,99)-prctile(filter_bp_s,1);
for i=1:length(p)
    curve_lick_predict=plot([p(i)/fr,p(i)/fr],[y_adjust_lick+10,y_adjust_lick+15],'k-','LineWidth',2);
    hold on;
end
for i=1:length(licktimeL)
    curve_lickL=plot([licktimeL(i),licktimeL(i)],[y_adjust_lick+5,y_adjust_lick+10],'b-','LineWidth',2);
    hold on;
end
for i=1:length(licktimeR)
    curve_lickR=plot([licktimeR(i),licktimeR(i)],[y_adjust_lick,y_adjust_lick+5],'r-','LineWidth',2);
    hold on;
end
%find licking using distance,method 2
T.distance(~rows)=0;%filter the points that likelihood<0.95
t_lick=(T.distance>threshold);%using histogram to find about 3.3 may be threshold for licking
ind_lick2=(diff(t_lick)>0);
ind_lick2=[0;ind_lick2];
p2=find(ind_lick2);
figure(figpose);
for i=1:length(p2)
    curve_lick_predict2=plot([p2(i)/fr,p2(i)/fr],[y_adjust_lick+15,y_adjust_lick+20],'-','Color',[0.5,0.5,0.5],'LineWidth',2);
    hold on;
end
plot([1,ts(end)],[threshold+y_adjust_d,threshold+y_adjust_d]);
%find licking using y,method 3
T.y(~rows)=y0;%filter the points that likelihood<0.95
t_lick=(T.y>thresholdy);%using histogram to find about 3.3 may be threshold for licking
ind_lick3=(diff(t_lick)>0);
ind_lick3=[0;ind_lick3];
p3=find(ind_lick3);
figure(figpose);
for i=1:length(p3)
    curve_lick_predict3=plot([p3(i)/fr,p3(i)/fr],[y_adjust_lick+20,y_adjust_lick+25],'-','Color',[1,0.5,0.5],'LineWidth',2);
    hold on;
end
plot([1,ts(end)],[thresholdy-min(T.y(rows))+y_adjust_y,thresholdy-min(T.y(rows))+y_adjust_y]);
legend([curve_lickL,curve_lickR,curve_lick_predict,curve_lick_predict2,curve_lick_predict3],{'lick-left','lick-right','lick-prediction-filter','lick-prediction-distance','lick-prediction-y'},'location','bestoutside');
set(gca,'FontSize',14);

%evaluate prediction based on video
% 1)compare recall and precision of different methods
lickTime=[licktimeL,licktimeR]';%ground true
prediction_filter=p/fr;
prediction_distance=p2/fr;
prediction_y=p3/fr;
Cprediction=cell(1,4);%each cell store a vector containing predicted/true licking
Cprediction{1}=prediction_filter;
Cprediction{2}=prediction_distance;
Cprediction{3}=prediction_y;
Cprediction{4}=lickTime;
setwin_ind=2;
figEvaluate=figure;
subplot(2,2,1);
win=0.1:0.1:1;
recall_filter=zeros(1,length(win));
recall_distance=zeros(1,length(win));
recall_y=zeros(1,length(win));
for i=1:length(win)
    window=[win(i),win(i)];
    recall_filter(i)=frecall(prediction_filter,lickTime,window);
    recall_distance(i)=frecall(prediction_distance,lickTime,window);
    recall_y(i)=frecall(prediction_y,lickTime,window);
end
plot(recall_filter,'r-');hold on;
plot(recall_distance,'b-');
plot(recall_y,'k-');
ylim=get(gca,'Ylim');
plot([setwin_ind,setwin_ind],[ylim(1),ylim(end)],'k-');
legend('prediction by filter','prediction by distance','prediction by y');
set(gca,'XTick',1:length(win),'XTickLabel',num2str(win'));
ylabel('recall');
xlabel('time window(s)');
set(gca,'FontSize',14);
subplot(2,2,2);
% bar([1,2,3],[recall_filter(setwin_ind),recall_distance(setwin_ind),recall_y(setwin_ind)]);
% set(gca,'XTick',[1,2,3],'XTickLabel',{'prediction by filter','prediction by distance','prediction by y'});
bar(1,recall_filter(setwin_ind),'FaceColor','r');hold on;
bar(2,recall_distance(setwin_ind),'FaceColor','b');
bar(3,recall_y(setwin_ind),'FaceColor','k');
set(gca,'XTick',[]);
set(gca,'FontSize',14);

subplot(2,2,3);
win=0.1:0.1:1;
recall_filter=zeros(1,length(win));
recall_distance=zeros(1,length(win));
recall_y=zeros(1,length(win));
for i=1:length(win)
    window=[win(i),win(i)];
    recall_filter(i)=frecall(lickTime,prediction_filter,window);
    recall_distance(i)=frecall(lickTime,prediction_distance,window);
    recall_y(i)=frecall(lickTime,prediction_y,window);
end
plot(recall_filter,'r-');hold on;
plot(recall_distance,'b-');
plot(recall_y,'k-');
ylim=get(gca,'Ylim');
plot([setwin_ind,setwin_ind],[ylim(1),ylim(end)],'k-');
legend('prediction by filter','prediction by distance','prediction by y');
set(gca,'XTick',1:length(win),'XTickLabel',num2str(win'));
ylabel('precision');
xlabel('time window(s)');
set(gca,'FontSize',14);
subplot(2,2,4);
bar(1,recall_filter(setwin_ind),'FaceColor','r');hold on;
bar(2,recall_distance(setwin_ind),'FaceColor','b');
bar(3,recall_y(setwin_ind),'FaceColor','k');
set(gca,'XTick',[]);
set(gca,'FontSize',14);

fig_histo_lickrate=figure;
timeMin=1000;
binStep=50;
trialLength=3000;
binSize=500;
[ ~,lickRate_predict_filter ] = fLickRateAligned( prediction_filter, lickTime,binSize,trialLength,timeMin, binStep);
[ ~,lickRate_predict_distance ] = fLickRateAligned( prediction_distance, lickTime,binSize,trialLength,timeMin, binStep);
[ ~,lickRate_predict_y ] = fLickRateAligned( prediction_y, lickTime,binSize,trialLength,timeMin, binStep);
plot(-timeMin:binStep:trialLength-timeMin-binSize+binStep,lickRate_predict_filter);
hold on;
plot(-timeMin:binStep:trialLength-timeMin-binSize+binStep,lickRate_predict_distance);
plot(-timeMin:binStep:trialLength-timeMin-binSize+binStep,lickRate_predict_y);
xlabel('Time from true lick(ms)');
ylabel('Lick rate(/s)');
legend('filtered','distance','y');
set(gca,'FontSize',14,'FontName','Arial');
% plot([0,1000],[times*std_fbp_s+400,times*std_fbp_s+400]);

%evaluate prediction based on video
% 2) trial type prediction
n_trial=length(time_trial_start);
arrayPredictIndTrialVio=zeros(n_trial,4);
indVio=cellfun(@(x) x.Action_choice==3,SessionResults);
for i=1:n_trial
    for ind_col=1:4
        temp=(Cprediction{ind_col}>=time_stim_onset(i) ).*(Cprediction{ind_col}<=time_delay_off(i));
        if  sum(temp)>0
            arrayPredictIndTrialVio(i,ind_col)=1;
        end
    end
end
TpredictIndTrialVio=array2table(arrayPredictIndTrialVio,'VariableNames',{'filter','distance','y','lick'});%table(filter,distance,y,lick);
figVio=figure;
subplot(1,2,1);
plot(1:n_trial,indVio,'ko');
hold on;
for i=1:4
    plot(1:n_trial,TpredictIndTrialVio{:,i});%here T{} change to vector, where T() is table
end
legend('vio','filter','distance','y','lick');
subplot(1,2,2);
scatter(5,sum(indVio)/n_trial);
hold on;
for i=1:4
    scatter(i,sum(TpredictIndTrialVio{:,i})/n_trial);
end
set(gca,'XTick',[1,2,3,4,5],'XTickLabel',{'filter','distance','y','lick','Vio'});
ylabel('Violation rate');
figure(figpose);%label those violated trial, from trial start to delay off
for i=1:n_trial
    if indVio(i)==0
        plot([time_trial_start(i), time_delay_off(i)],[y_adjust_lick+26,y_adjust_lick+26],'Color',[0.5,0.5,0.5],'LineWidth',2);%label non-vio trial with gray lines
    else 
        plot([time_trial_start(i), time_delay_off(i)],[y_adjust_lick+26,y_adjust_lick+26],'Color',[1,0.5,0],'LineWidth',2);%label vio trial with orange line
    end
    for ind_col=1:4
        if TpredictIndTrialVio{i,ind_col}==0
            plot([time_stim_onset(i), time_delay_off(i)],[y_adjust_lick+26+3*ind_col,y_adjust_lick+26+3*ind_col],'Color',[0.5,0.5,0.5],'LineWidth',2);
        else
            plot([time_stim_onset(i), time_delay_off(i)],[y_adjust_lick+26+3*ind_col,y_adjust_lick+26+3*ind_col],'Color',[1,0.5,0],'LineWidth',2);
        end
    end
    text(time_trial_start(i),y_adjust_lick+45,['trial',num2str(i)]);
end
%{
figure;%compare filter
plot(ts,T.x);hold on;
plot(ts,T.y);
plot(ts,Tf.x);
plot(ts,Tf.y);
legend('x','y','x_filter','y_filter');
%}
%{
figure;%compare diff model
% rows1=Tt.likelihood>0.02;
% rows2=Tpc.likelihood>0.02;
rows1=1:length(Tt.likelihood);
rows2=1:length(Tpc.likelihood);
subplot(2,1,1);
plot(ts_t(rows1),Tt.x(rows1));
hold on;
plot(ts_t(rows1),Tt.y(rows1));
plot(ts_t(rows1),Tt.y_3(rows1));
plot(ts(rows2),Tp.x(rows2));
plot(ts(rows2),Tp.y(rows2));
plot(ts(rows2),Tp.y_3(rows2));
legend('tongue_x-t','tongue_y-t','lickport_y-t','tongue_x','tongue_y','lickport_y');
set(gca,'Xlim',[0,15]);
set(gca,'FontSize',14);
subplot(2,1,2);
plot(ts_t,Tt.likelihood);hold on;
plot(ts_t,Tt.likelihood_3);

plot(ts,Tp.likelihood);
plot(ts,Tp.likelihood_3);
set(gca,'Xlim',[0,15]);
legend('tongue-t','lickport-t','tongue','lickport');
set(gca,'FontSize',14);
%}
function [f0]=fmode(f)
[N,edges,bin] = histcounts(f,100);
f0 = edges(N == max(N));
if length(f0)==2
    f0=mean(f0);
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
function [ lickMat,lickRate ] = fLickRateAligned( lickTime, base,binSize,trialLength,timeMin, binStep)
sampleNum=floor((trialLength-binSize)/binStep)+1;%单位ms
lickMat=zeros(length(base),trialLength);%第一维trials,第二维单个trial长度
lickRate=zeros(sampleNum,1);%第一维代表trial长度/采样数量
for i=1:length(base)
    lickT=(lickTime-base(i))*1000;%单位s
    for j=1:length(lickT)
        if ~isnan(lickT(j)) && lickT(j)>-timeMin+1 && lickT(j)<trialLength-timeMin-1
            lickMat(i,round(lickT(j))+timeMin)=1;
        end
    end
end
for j=1:sampleNum
    lickRate(j)=nansum(nansum(lickMat(:,(binStep*j-binStep+1):(binStep*j-binStep+binSize))))*1000/binSize/length(base);
end

end