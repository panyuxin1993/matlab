%plot opto effect in 2 ways:
% 1)compare baseline/opto groups 2)compare adjacent blocks

%load data
clear;
fileFolder='D:\xulab\behavior\practice\LTY\LTY_ChR2_01\opto_effect';
dirmat=strcat(fileFolder,'\*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
sort(filenames);
nfiles=length(filenames);
animal=cellfun(@(x) x(1:9),filenames,'UniformOutput',false);
%define dataset
dataGroup=zeros(nfiles,3); %store data for different group(left opto, baseline, right opto),data is p(answer right) or p(licking number right)
dataBlock=zeros(nfiles,5); %store data for different block(baseline, 2nd block, 3rd block, etc),data is p(choice of right) if the first opto side is right, transform others to this kind
datatype='licking';%'licking'%define which type of parameters to calculate and strore above
micegroup=zeros(nfiles,1);
for n=1:nfiles
    load(strcat(fileFolder,'\',filenames{n}));
    ntrials=length(SessionResults);
    nblock=cellfun(@(x) x.Trial_Num==1,SessionResults);
    blockInd=find(nblock);
    blockInd=[blockInd ntrials+1];
    nblocks=sum(nblock);
    optoTypeIndex=zeros(ntrials,3);
    optoTypeIndex(1:200,2)=1;%1:200trials are baseline
    optoIndex=cellfun(@(x) x.Time_optoStimOnset~=0,SessionResults);
    choice=cellfun(@(x) x.Action_choice,SessionResults);
    choice=double(choice);
    j=1;
    for i=1:length(blockInd)-1
        if i>1%i=1 first block of a session, is baseline
            temp=sum(optoIndex(blockInd(i):blockInd(i+1)-1).*choice(blockInd(i):blockInd(i+1)-1));
            if temp==0 %left opto
                optoTypeIndex(blockInd(i):blockInd(i+1)-1,1)=1;
            else % left opto
                optoTypeIndex(blockInd(i):blockInd(i+1)-1,3)=1;
            end
        end
        %calculate dataBlock
        if strcmp(datatype,'answer')%calculate answer
            dataBlock(n,j)=sum(choice(blockInd(i):blockInd(i+1)-1)==1)/sum(choice(blockInd(i):blockInd(i+1)-1)<2);
        else%calculate licking number
            leftLick=cellfun(@(x) x.Action_numLickLeft,SessionResults);
            rightLick=cellfun(@(x) x.Action_numLickRight,SessionResults);
            dataBlock(n,j)=sum(rightLick(blockInd(i):blockInd(i+1)-1))/sum(rightLick(blockInd(i):blockInd(i+1)-1)+leftLick(blockInd(i):blockInd(i+1)-1));
        end
        j=j+1;
    end
    if optoTypeIndex(201,1)==1%if first opto side is left, transform
        %dataBlock(n,:)=1-dataBlock(n,:);
        micegroup(n,1)=1; %1-->first optoblock is left; 0--> first opto block is right
    end
    %calculate dataGroup
    if strcmp(datatype,'answer')%calculate answer
        for i=1:3
            dataGroup(n,i)=sum((choice==1)'.*optoTypeIndex(:,i))/sum((choice<2)'.*optoTypeIndex(:,i));
        end
        ylabelstr='P(right choice)';
    else%calculate licking number
        leftLick=cellfun(@(x) x.Action_numLickLeft,SessionResults);
        rightLick=cellfun(@(x) x.Action_numLickRight,SessionResults);
        leftLick=double(leftLick);
        rightLick=double(rightLick);
        for i=1:3
            dataGroup(n,i)=sum(rightLick*optoTypeIndex(:,i))/sum((rightLick+leftLick)*optoTypeIndex(:,i));
        end
        ylabelstr='P(right licking)';
    end
end
%plot figures
A=figure;%plot by groups
B=figure;%plot by blocks
figure(A);
for i=1:size(dataGroup,1)
    plot(1:size(dataGroup,2),dataGroup(i,:),'k-','linewidth',1);
    hold on;
    scatter(1:size(dataGroup,2),dataGroup(i,:),10,'k');
end
ylabel(ylabelstr,'FontName','Arial','FontSize',14);
box off; %取消图标外边框
set(gca, 'XLim',[0 size(dataGroup,2)+1]);
set(gca,'XTick',1:size(dataGroup,2));%x轴只标注几个边界点
set(gca,'XTickLabel',{'opto left','baseline','opto right'});%保留到百位数的精度
set(gca,'FontName','Arial','FontSize',14);
figure(B);
set(gcf,'position',[100,100,900,600]);
for i=1:size(dataBlock,1)
    if micegroup(i,1)==0 %first optoblock is right
        curve_r=plot(1:size(dataBlock,2),dataBlock(i,:),'k-','linewidth',2);
        hold on;
        scatter(1:size(dataBlock,2),dataBlock(i,:),10,'k');
    else
        curve_l=plot(1:size(dataBlock,2),dataBlock(i,:),'k--','linewidth',2);
        hold on;
        scatter(1:size(dataBlock,2),dataBlock(i,:),10,'k');
    end
    text(1,dataBlock(i,1),animal(i),'FontName','Arial','FontSize',14);
end
ylabel(ylabelstr,'FontName','Arial','FontSize',14);
legend([curve_r curve_l],'opto right first','opto left first');
box off; %取消图标外边框
set(gca, 'XLim',[0 size(dataBlock,2)+1]);
set(gca,'XTick',1:size(dataBlock,2));%x轴只标注几个边界点
set(gca,'XTickLabel',{'baseline','1 st Opto','2nd Opto','3rd Opto','4th Opto'});%保留到百位数的精度
set(gca,'FontName','Arial','FontSize',14);
