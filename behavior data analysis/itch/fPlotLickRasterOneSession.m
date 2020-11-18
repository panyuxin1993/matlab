function [strtitle]=fPlotLickRasterOneSession(lickTime,trialType,resultType,sortType,filename,animal,maxdelay,nresample,rasterSize)
%trialTypeѡ��trialΪ��0��1,��ѡ������[0 1];ismember(lickTime.trialType,trialType)
%resultTypeѡ��result�ǡ�withVio'����trial����noVio'����violation
%sortTypeѡ���Ƿ�sort��'sort'��ʾsort��'raw'��ʾ��sort
if strcmp(nresample,'resample')
    nsample=unidrnd(length(lickTime),1,200);%trial ����̫�࣬��ȡ����һ����
    nsample=sort(nsample);
    lickTime=lickTime(nsample);
else
    nsample=1:length(lickTime);
end
if strcmp(sortType,'sort')
    [~,index]=sort([lickTime.delay]);
elseif strcmp(sortType,'raw')
    index=1:1:length(lickTime);%when sort don't need
end
yi=1;%����lick raster�ڵڼ��У����ֻ���������Ҳ��trial�� ���԰��м�������ȥ�������¶ѵ�
for i=1:length(nsample) %i��trial
    if ismember(lickTime(index(i)).choice,trialType)
%         line([lickTime(index(i)).stimOn,lickTime(index(i)).stimOff],[i-0.5,i+0.5],'color',[0.9,0.9,0.9],'linewidth',2);%stim, indicated with grey
%         line([lickTime(index(i)).stimOn,lickTime(index(i)).stimOn],[yi-0.5,yi+0.5],'color','k','linewidth',2);%stimOn
%         hold on;
%         line([lickTime(index(i)).stimOff,lickTime(index(i)).stimOff],[yi-0.5,yi+0.5],'color','k','linewidth',2);%stimOff
%         hold on;
        line([lickTime(index(i)).go,lickTime(index(i)).go],[yi-0.5,yi+0.5],'color','k','linewidth',2);%go cue
        hold on;
        if strcmp(resultType,'noVio') && (lickTime(index(i)).choice==3)%����vio��ͼ
           continue;
        else
            for jl=1:length(lickTime(index(i)).lickleft)
                line([lickTime(index(i)).lickleft(jl),lickTime(index(i)).lickleft(jl)],[yi-0.5,yi+0.5],'color','b','linewidth',rasterSize);%left lick
                hold on;
            end
            for jr=1:length(lickTime(index(i)).lickright)
                line([lickTime(index(i)).lickright(jr),lickTime(index(i)).lickright(jr)],[yi-0.5,yi+0.5],'color','r','linewidth',rasterSize);%right lick
                hold on;
            end
        end
        yi=yi+1;
    end
end
%�����ļ��������date,animal��Ϣ
pdate=strfind(filename,'_');
pend=strfind(filename,'.');
date=filename(pdate(1)+1:pend(end)-1);
strtitle=strcat(animal,'-',date);

set(gca, 'YLim',[0 yi]);
set(gca, 'XLim',[-maxdelay-500 maxdelay+2500]);
% xlabel('time(ms) from stimulus onset','FontName','Arial','FontSize',14);
xlabel('time(ms) from go cue','FontName','Arial','FontSize',14);
ylabel('trials','FontName','Arial','FontSize',14);
title(strtitle,'FontName','Arial','FontSize',14);
plot([0,0],[0,1000],'k','linewidth',1);
hold on;
set(gca,'FontName','Arial','FontSize',14);
end

