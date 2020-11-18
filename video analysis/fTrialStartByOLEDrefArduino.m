function [ trial_start_frame_final,trial_start_time_final ] = fTrialStartByOLEDrefArduino( T, framerate ,trialStartBeh,varargin)
%FTRIALSTARTBYOLED_DY using OLED brightness to predict trial start in the
%video, using trial start time in beh file for reference. Acknowledge for
%primary attempts by Yun Du. Upgraded from fTrialStartByOLED_dy to analyze
%sessions with pause where ITI is not a constant and thus can not be used
%to verify prediction results. By aligning predicted trial starts with beh
%file records, the accuracy of prediction may be improved.
%Input-
%   T- table read from the csv file containing brightness change extracted
%   from the .tiff file that transformed from raw .webm file
%   framerate- frame rate of the video, usually 24.
%   path- file path saving the result
%   trialStartBeh- vector, time(s) of trial start, as a template for
%   prediction to align
%   method- 'meandiff'(default)|'diff'|'abs', using absolute value or difference
%   between frames to decide when trial start,'meandiff' combine mean and
%   diff to calculate
%Output-
%   trial_start_time_final- n-by-1 vector, indicating the start time(s) of
%   the n trials from session start.
%   trial_start_frame_final- n-by-1 vector, indicating the start frame of
%   the n trials from session start.
if isempty(varargin)
    method='meandiff';
else
    method=varargin{1};
end

T.ts=T.X/framerate;%time stamp, using this to plot easier to be understand
figure_start=figure;%各种预测画在一张图上，人眼进行质量检测
figure(figure_start);
original=plot(T.ts,T.Y,'k');%原始LED亮度图
%led on预测
if strcmp(method,'meandiff')
    binstep=1;
    binsize=floor(1000/framerate);%1s as a bin
    paramd=zeros(1,size(T,1)-binsize*2+1);
    for i=0:length(paramd)-1
        paramd(i+1)=mean(T.Y(i+binsize+1:i+binsize*2))-mean(T.Y(i+1:i+binsize));
    end
    figmd=figure;
    plot(paramd);
    threshold=2;%based on figmd and mannually set
    paramdflag=paramd>threshold;
    [segparamdflag, segparamd]=fSegment2Cell(paramdflag,paramd);
    value_max=cellfun(@max,segparamd,'UniformOutput',false);
    indLedOnFlagCell=cellfun(@(x,y) (x(1)==paramd).*y,value_max,segparamdflag,'UniformOutput',false);
    indLedOnFlagMat=cell2mat(indLedOnFlagCell);
    indLedOnFlagVector=logical(sum(indLedOnFlagMat,1));
    indLedOnFlagVector=[false(1,binsize),indLedOnFlagVector,false(1,binsize-1)];
    led_on=indLedOnFlagVector;
    %visualize
    x_on=T.X(led_on);%根据亮度剧变预测led on时间
    y_on=T.Y(led_on);
    ts_on=T.ts(led_on);
    figure(figure_start);
    hold on
    start_by_meandelta=scatter(ts_on,y_on,'b','filled');%作图,根据亮度剧变预测led on时间
    %using trial interval from .beh file to refined the results
    [trial_start_frame_final,trial_start_time_final,y_OLED_final ]= fRefineTrialStartRefBeh( x_on, framerate ,trialStartBeh,y_on );
    %plot refined data to check
    figure(figure_start);
    start_refined_by_behIT=scatter(trial_start_frame_final/framerate,y_OLED_final,'r','filled');
    for i=1:length(trial_start_frame_final)
        text(trial_start_frame_final(i)/framerate,y_OLED_final(i)+1,strcat('trial-',num2str(i)));
    end
    legend([start_by_meandelta,start_refined_by_behIT],'start found by mean diff','start refined by trial interval of .beh');
    box off;
%     close(figure_start);
    
else
    if strcmp(method,'diff')
        led_on=(diff(T.Y)>6);
        led_on=[false;led_on];
        %     no_led_on=sum(led_on);
    elseif strcmp(method,'abs')
        fighistoRaw=figure;
        histogram(T.Y);
        threshold=11.5;
        led_on_time=T.Y>threshold;
        led_on=(diff(led_on_time)>0);
        
    end
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
    legend([start_by_delta,exclusion,original,on_predicted,curve_trial_start_final],{'start time predicted by delta Y','exclused time by sequence','brightness detected by ImageJ','predicted start time','led on final'});
end
end

