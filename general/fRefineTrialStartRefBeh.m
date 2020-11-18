function  [ trial_start_frame_final,trial_start_time_final,y_OLED_final ] = fRefineTrialStartRefBeh( indFrameTrialStart, framerate ,trialStartBeh,varargin )
%FREFINETRIALSTARTREFBEH %using trial interval from .beh file to refined the results
%Input-
%   indFrameTrialStart- vector of frame number of trial start, length means
%   total trial number
%   framerate- frame rate of the video, usually 24.
%   path- file path saving the result
%   trialStartBeh- vector, time(s) of trial start, as a template for
%   prediction to align
%   
%Output-
%   trial_start_time_final- n-by-1 vector, indicating the start time(s) of
%   the n trials from session start.
%   trial_start_frame_final- n-by-1 vector, indicating the start frame of
%   the n trials from session start.
    if ~isempty(varargin)
        y_on=varargin{1};
    else
        y_on= indFrameTrialStart;
    end
    x_on= indFrameTrialStart;
    ts_on = indFrameTrialStart/framerate;
    figTI=figure;
    TIBeh=diff(trialStartBeh);
    ntrialBeh=length(trialStartBeh);
    TIOLED=diff(ts_on);
    plot(TIBeh,'k-');hold on;
    indexesRefInData = fBlastAnalog( TIBeh,TIOLED,[0.99,1.01] );%using homemade blast method to match trials from .beh and OLED
    ts_TIOLED=1:length(TIBeh);
    indTIOLED=~isnan(indexesRefInData);
    tempind=indexesRefInData(indTIOLED);
    if issorted(tempind)
        figure(figTI);
        plot(ts_TIOLED(indTIOLED),TIOLED(tempind),'b-');
        legend('trial interval from .beh file','trial interval from OLED');
    else
        warning('output is not sorted, may have error. Please check.');
    end
    trial_start_time_final=zeros(ntrialBeh,1);
    trial_start_frame_final=zeros(ntrialBeh,1);
    y_OLED_final=zeros(ntrialBeh,1);
    for i=1:ntrialBeh-1
        if ~isnan(indexesRefInData(i))
            trial_start_frame_final(i)=x_on(indexesRefInData(i));
            y_OLED_final(i)=y_on(indexesRefInData(i));
        elseif isnan(indexesRefInData(i)) && i==1
            trial_start_frame_final(i)=nan;
            y_OLED_final(i)=nan;
        elseif isnan(indexesRefInData(i)) && ~isnan(indexesRefInData(i-1))
            trial_start_frame_final(i)=x_on(indexesRefInData(i-1)+1);
            y_OLED_final(i)=y_on(indexesRefInData(i-1)+1);
        elseif isnan(indexesRefInData(i)) && isnan(indexesRefInData(i-1))
            trial_start_frame_final(i)=nan;
            y_OLED_final(i)=nan;
        end 
    end
    if  ~isnan(indexesRefInData(end))
        trial_start_frame_final(end)=x_on(indexesRefInData(end)+1);
        y_OLED_final(end)=y_on(indexesRefInData(end)+1);
    else
        trial_start_frame_final(end)=nan;
        y_OLED_final(i)=nan;
    end

%     close(figTI);
end

