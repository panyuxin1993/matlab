% matching law- R2/R1=b*(Rf1/Rf2)^s
% or - log2(R2/R1)= log2(b)+s*log2(Rf1/Rf2)
% for each mouse, fitting matching law and plot relationship
% each dots represent one block
ffun=fittype('log2(b)+s*log2Rf','independent','log2Rf','dependent','log2R');
% readfile
rootpath='D:\xulab\behavior\practice\LTY\';
animal='LTY04';
fileFolder=strcat(rootpath, animal,'\');
dirmat=strcat(fileFolder,'*.mat');%include all files in a folder named as an animal
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=sort(dircell(1,:));
nsession=length(filenames);
nblock=[];
for i=1:nsession
    load(strcat(fileFolder,filenames{i}));
    iblock=cellfun(@(x) x.Trial_Num==1,SessionResults);
    nblock=[nblock iblock];
end
matching_law_data=zeros(sum(nblock),3);%each column, ratio of rigt vs. left 1)trials, 2)answers, 3)lickings
for i=1:nsession
    load(strcat(fileFolder,filenames{i}));
    for j=1:nblock(i)
        pleft=double(SessionSettings{1, 1}.leftProb);%note here assume that each session one fixed leftProb
        if pleft==0
            warning('leftProb is zero');
        end
        matching_law_data(i,1)=(100-pleft)/pleft;
        nanswer=fGetAnsNum(SessionResults);
        nlicking=fGetLickNum(SessionResults);
        matching_law_data(i,2)=nanswer(2)/nanswer(1);
        matching_law_data(i,3)=nlicking(2)/nlicking(1);
    end
end
%fitting
f_ans=fit(matching_law_data(:,1),matching_law_data(:,2),ffun,'Startpoint',[1,1]);
f_lick=fit(matching_law_data(:,1),matching_law_data(:,3),ffun,'Startpoint',[1,1]);
x=min(matching_law_data(:,1)):0.01:max(matching_law_data(:,1));
y=zeros(length(x),2);
y(:,1)=f_ans(x);
y(:,2)=f_lick(x);
%plot
figure;
set(gcf,'position',[0,0,800,400]);
ystr={' answer',' licking'};
for i=1:2
    subplot(1,2,i);
    plot(x,y(:,i));
    hold on;
    scatter(matching_law_data(:,1),matching_law_data(:,i+1));
    title('Matching Law relationship');
    xlabel('Ratio of two side reward probability(log2)');
    ylabel(strcat('Ratio of two side',ystr{i},' probability(log2)'));
    set(gca,'FontName','Arial','FontSize',14);
end

