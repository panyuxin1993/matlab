% matching law- R2/R1=b*(Rf1/Rf2)^s
% or - log2(R2/R1)= log2(b)+s*log2(Rf1/Rf2)
% for each mouse, fitting matching law and plot relationship
% each dots represent one block
ffun=fittype('log2(b)+s*log2Rf','independent','log2Rf','dependent','log2R');
% readfile
fileFolder='D:\xulab\behavior\practice\LTY\matching_law\';
dirmat=strcat(fileFolder,'*.mat');%include all files in a folder named as an animal
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=sort(dircell(1,:));
nsession=length(filenames);
choice_and_probability_data=cell(nsession,1);
nblock=[];
for i=1:nsession
    load(strcat(fileFolder,filenames{i}));
    iblock=cellfun(@(x) x.Trial_Num==1,SessionResults);
    iblock=sum(iblock);
    nblock=[nblock iblock];
end

for i=1:nsession
    load(strcat(fileFolder,filenames{i}));
    ntrial=length(SessionResults);
    temp_trialType=double(cellfun(@(x) x.Trial_Type,SessionResults));
    choice=cellfun(@(x) x.Action_choice,SessionResults);
    leftLick=cellfun(@(x) x.Action_numLickLeft,SessionResults);
    rightLick=cellfun(@(x) x.Action_numLickRight,SessionResults);
    choice_and_probability_data{i}=zeros(nblock(i),3);%each column, ratio of rigt vs. left 1)trials, 2)answers, 3)lickings
    for j=1:nblock(i)
        %pleft=double(SessionSettings{1, 1}.leftProb);%note here assume that each session one fixed leftProb
        jblock=j;
        iblock=cellfun(@(x) x.Trial_Num==1,SessionResults);
        blockInd=find(iblock);
        blockInd=[blockInd ntrial+1];
%         pleft=100*sum(temp_trialType(blockInd(j):blockInd(j+1)-1)==0)/(blockInd(j+1)-blockInd(j));
%         if pleft==0
%             warning('leftProb is zero');
%             choice_and_probability_data{i}(jblock,1)=10;
%         else
%             choice_and_probability_data{i}(jblock,1)=(100-pleft)/pleft;
%         end
        choice_and_probability_data{i}(jblock,2)=sum(choice(blockInd(j):blockInd(j+1)-1)==1)/sum(choice(blockInd(j):blockInd(j+1)-1)<2);
        choice_and_probability_data{i}(jblock,3)=sum(rightLick(blockInd(j):blockInd(j+1)-1))/(sum(leftLick(blockInd(j):blockInd(j+1)-1))+sum(rightLick(blockInd(j):blockInd(j+1)-1)));
    end
end
%plot
figure;
animal=cellfun(@(x) x(1:9),filenames,'UniformOutput',false);
set(gcf,'position',[0,0,800,400]);
ystr={' answer',' licking'};
for i=1:2
    subplot(1,2,i);
    for j=1:length(choice_and_probability_data)
        scatter(1:size(choice_and_probability_data{j},1),choice_and_probability_data{j}(:,i+1),20,'k');
        hold on;
        plot(1:size(choice_and_probability_data{j},1),choice_and_probability_data{j}(:,i+1),'k-');
        hold on;
        text(1,choice_and_probability_data{j}(1,i+1),animal(j),'FontName','Arial','FontSize',14);
    end
    xlabel('block');
    ylabel(strcat('probability of right',ystr{i}));
    set(gca, 'XLim',[0 max(nblock)+1]);
    set(gca,'FontName','Arial','FontSize',14);
end

