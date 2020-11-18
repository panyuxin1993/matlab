animal='pyx046';
fileFolder=strcat('D:\xulab\behavior\',animal);
dirmat=strcat(fileFolder,'\*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
sort(filenames);%将文件按照训练先后顺序排列起来
firstLick=cell(length(filenames),1);%每个cell为一个列向量，存放每个trial的firstlick，miss trial定义为nan
%读数据
for i=1:length(filenames)
    load(strcat(fileFolder,'\',filenames{i}));
    firstLick{i}=zeros(length(SessionResults),2);%第一列是相对于delayOnset，第二列是delay长度
    for j=1:length(SessionResults)
        lickLTp=strfind(SessionResults{1,j}.Action_lickTimeLeft,'|');
        lickRTp=strfind(SessionResults{1,j}.Action_lickTimeRight,'|');
        if length(lickLTp)<=2
            firstLickLT=nan;
        else
            firstLickLT=str2double(SessionResults{1,j}.Action_lickTimeLeft(lickLTp(2)+1:lickLTp(3)-1));
        end
        if length(lickRTp)<=2
            firstLickRT=nan;
        else
            firstLickRT=str2double(SessionResults{1,j}.Action_lickTimeRight(lickRTp(2)+1:lickRTp(3)-1));
        end
        if firstLickLT>firstLickRT || isnan(firstLickLT)
            firstLickTime=firstLickRT;
        else
            firstLickTime=firstLickLT;
        end            
        firstLick{i}(j,1) = firstLickTime-SessionResults{1,j}.Time_delayOnset;
        if firstLick{i}(j,1)<-500
            i;
        end
        firstLick{i}(j,2) = SessionResults{1,j}.Delay_duration;
    end
end
%plot
A=figure(1);
hold on;
subplot(1,2,1);%每个trial first lick相对于delayOnset
for i=1:length(filenames)
    x=i*ones(size(firstLick{i},1),1);
    scatter(x,firstLick{i}(:,1));
    hold on;
    plot([0,length(filenames)],[0,0]);
    hold on;
    plot([0,length(filenames)],[-500,-500]);%stimOnset
end
set(gca, 'YLim',[-2000 4500]);
%boxplot(firstLick{:}(:,1));
subplot(1,2,2);%每个trial first lick相对于delayOffset
for i=1:length(filenames)
    x=i*ones(size(firstLick{i},1),1);
    scatter(x,firstLick{i}(:,1)-firstLick{i}(:,2));
    hold on;
    %boxplot(firstLick{i}(:,2));
    %hold on;
    plot([0,length(filenames)],[0,0]);
    hold on;
    %plot([0,length(filenames)],[-500,-500]);
end
set(gca, 'YLim',[-2000 4500]);