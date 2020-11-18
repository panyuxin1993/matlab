%plot performance within session,i.e. 1-1000trials,every 20 trials as a bin
dir='D:\xulab\behavior\CNO experiment\infusion\unilateral\';
% dir='D:\xulab\behavior\CNO experiment\infusion\bilateral\';
% dir='D:\xulab\behavior\CNO experiment\ip\';
titleuni={'all','ipsi','contra'};
titlebi={'all','left','right'};
group={'saline', 'CNO'};
colorall={'k','r'};
colorindi={[0.7,0.7,0.7],[1,0.7,0.7]};
performance=cell(length(group),1);
correctrate=cell(length(group),3);%第一维表示saline/CNO,第二维表示总，ipsi tirlas,contra trials,
bin=20;%每20个取平均
A=figure;

for n=1:length(group)
    [animal_name,dataChoiceResult]=fFindChoiceMat(strcat(dir,group{n}));
    nsample=size(dataChoiceResult,1)-1;
    performance{n}=zeros(nsample,1000);%每行一个session
    side=cell(nsample,1);%用于记录side
    for i=1:nsample
        for j=1:size(dataChoiceResult{i},1)-bin%超过这个范围就是不到bin个trial的平均了
            performance{n}(i,j)=dataChoiceResult{i}(j,2);
        end
        if size(dataChoiceResult{i},1)==1000
            continue;
        else
            for j=size(dataChoiceResult{i},1)-bin+1:size(performance{n},2)
                performance{n}(i,j)=nan;
            end
        end
        %以下针对uni情况进行分类ipsi/contra；否则就是ipsi代表左边，contra代表右边
        dataChoiceResult{i}(:,1)=dataChoiceResult{i}(:,3);%左右调换的话就用第六列的声音对应的左右替换第一列频率
        protate=logical((dataChoiceResult{i}(:,6)==2)+(dataChoiceResult{i}(:,6)==3));%这些trials需要旋转，也就是上下，左右各颠倒一次；对于bilateral/ip则是0
        dataChoiceResult{i}(protate,1)=7-dataChoiceResult{i}(protate,1);
        ipsi=(dataChoiceResult{i}(:,1)<=3);
        contra=(dataChoiceResult{i}(:,1)>=3);
        all=logical(ipsi+contra);
        side{i}=[all, ipsi, contra];
    end
    for niorc=1:3
        correctrate{n,niorc}=zeros(nsample+1,1000-bin+1);%每行一个session,最后一个为mean
        for i=1:nsample
            for j=1:size(side{i},1)-bin+1 %此值小于size(correctrate{n,niorc},2)
                cor=(performance{n}(i,j:j+bin-1)==1);
                err=(performance{n}(i,j:j+bin-1)==2);
                sidenow=side{i,1}(j:j+bin-1,niorc);
                correctrate{n,niorc}(i,j)=sum(cor'.*sidenow)/sum((err+cor)'.*sidenow);
            end
            if size(side{i},1)==size(correctrate{n,niorc},2)
                continue;
            else
                for j=size(side{i},1)-bin:size(correctrate{n,niorc},2)
                    correctrate{n,niorc}(i,j)=nan;
                end
            end
        end
        for i=1:size(correctrate{n,niorc},2)
            correctrate{n,niorc}(end,i)=nanmean(correctrate{n,niorc}(1:end-1,i));
        end
    end
    %画individual曲线
    figure(A);
    for niorc=1:3
        subplot(1,3,niorc);
        for i=1:nsample
            plot(1:size(correctrate{n,niorc},2),correctrate{n,niorc}(i,:),'Color',colorindi{n},'LineWidth',0.5);
            hold on;
        end
        ylabel('Correct Rate','FontName','Arial','FontSize',14);
        if contains(dir,'uni')
            title(titleuni{niorc},'FontName','Arial','FontSize',14);
        else
            title(titlebi{niorc},'FontName','Arial','FontSize',14);
        end
    end
end
for niorc=1:3
    figure(A);
    subplot(1,3,niorc);
    %画总的曲线
    for n=1:length(group)
        plot(1:size(correctrate{n,niorc},2),correctrate{n,niorc}(end,:),'Color',colorall{n},'LineWidth',2);
        hold on;
    end
end