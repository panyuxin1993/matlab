function [strtitle]=fPlotLickRasterOneSession(lickTime,trialType,resultType,sortType,filename,animal,maxdelay,nresample,rasterSize)
%trialType选择trial为左0右1,均选择则用[0 1];ismember(lickTime.trialType,trialType)
%resultType选择result是‘withVio'所有trial，’noVio'不画violation
%sortType选择是否sort，'sort'表示sort，'raw'表示不sort
if strcmp(nresample,'resample')
    nsample=unidrnd(length(lickTime),1,200);%trial 数量太多，抽取其中一部分
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
yi=1;%控制lick raster在第几行，如果只画左侧或者右侧的trial， 可以把中间空余的行去掉，向下堆叠
for i=1:length(nsample) %i个trial
    if ismember(lickTime(index(i)).choice,trialType)
%         line([lickTime(index(i)).stimOn,lickTime(index(i)).stimOff],[i-0.5,i+0.5],'color',[0.9,0.9,0.9],'linewidth',2);%stim, indicated with grey
%         line([lickTime(index(i)).stimOn,lickTime(index(i)).stimOn],[yi-0.5,yi+0.5],'color','k','linewidth',2);%stimOn
%         hold on;
%         line([lickTime(index(i)).stimOff,lickTime(index(i)).stimOff],[yi-0.5,yi+0.5],'color','k','linewidth',2);%stimOff
%         hold on;
        line([lickTime(index(i)).go,lickTime(index(i)).go],[yi-0.5,yi+0.5],'color','k','linewidth',2);%go cue
        hold on;
        if strcmp(resultType,'noVio') && (lickTime(index(i)).choice==3)%不画vio的图
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
%分析文件名，获得date,animal信息
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

