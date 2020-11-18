function [onset,n_change] = fOnsetPhaseChangeExamAUCCase(series,AUC,ts,threshold,celltypestr,varargin)
%FONSETPHASECHANGE find onset when a dominant stage changes(I call it phase
%change). The idea is to find the longest continuous duration of one stage
%and define the onset of that
%Input-
%   series- 0-1-vector to know the onset, to find the start of '1' series
%   AUC- AUC of corresponding case, just used for seeing how it look
%   like,you can also use p-value instead
%   celltypestr- used for name of saved figures
%   ts- timestamp used for plot
%   leastduration- 5(default), a threshold for least continuous duration
%Output-
%   onset-a number of when start the longest continuous duration
%   n_change-number indicating how many times phase change happen, given
%   the criteria that duration is longer than least duration
if isempty(varargin)
    leastduration=5;
else
    leastduration=varargin{1};
end
diffSeries=diff(series);
diffSeries=reshape(diffSeries,1,[]);%to a row vector
diffSeries=[0,diffSeries];
tchange=find(diffSeries~=0);
if isempty(tchange)
    onset=[];
    n_change=0;
else
    if diffSeries(tchange(1))<0
        tchange=[1,tchange];
    end
    if diffSeries(tchange(end))>0
        tchange=[tchange,length(series)];
    end
end
%now, n-pairs, each up-down
ind_up=1:2:length(tchange);
ind_down=2:2:length(tchange);
duration=tchange(ind_down)-tchange(ind_up);
n_change=sum(duration>leastduration);
ind_max=find((max(duration)==duration).*(duration>=leastduration));%this can result more than 1 vector
ind_1st=find(duration>=leastduration);
goframe=find(ts==2);%2s from stimuli onset as longest frame
if isempty(ind_max)
    onset=goframe;
else
    onsets=tchange(2*ind_max-1);
%     onsets=tchange(2*ind_1st-1);
    onsets=onsets(onsets>threshold);
    if isempty(onsets)
        onset=goframe;
    else
        onset=onsets(1);
    end
end

%show period of changed phase as shadows
if ~isempty(AUC)
    for i=1:length(ind_up)
        x=[tchange(ind_up(i)),tchange(ind_down(i)),tchange(ind_down(i)),tchange(ind_up(i))];
        y=[0,0,1,1];
        if duration(i)==max(duration) && duration(i)>=leastduration
            patch(ts(x),y,[1,0.5,0.5]);
        elseif duration(i)>=leastduration
            patch(ts(x),y,[0.5,0.5,1]);
        else
            patch(ts(x),y,[0.5,0.5,0.5]);
        end
        hold on;
        if ~isempty(onset)
            plot([ts(onset),ts(onset)],[0,1],'Color',[1,0,1],'LineWidth',2);
        end
    end
    plot(ts,AUC,'-k');hold on;
    x_lim=get(gca,'Xlim');
    set(gca,'Ylim',[0,1]);
    plot([0,0],[0,1],'k--');
    plot([-0.5,-0.5],[0,1],'k--');
    title(celltypestr);
    %xlabel('time(s) from stimuli onset');
%     saveas(tempfig,['F:\FP\summary\AUC-cases2\',celltypestr,'.png'],'png')
%     set(tempfig,'PaperPosition',[1 1 2 1.5 ]);
%     set(gca,'Xlim',x_lim);
%     saveas(tempfig,['F:\FP\summary\AUC-cases2\',celltypestr,'.pdf'],'pdf')
    %close(tempfig);
end
end

