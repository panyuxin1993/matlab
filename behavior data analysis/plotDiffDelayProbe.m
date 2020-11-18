%������ͬ���ȵ�delay�Ĳ���probe���������
clear;
diff_delay_num=2;%���岻ͬ���ȵ�delay��3���Ǳ�ʾ�ֳɶ��г�����,2���Ƿֳɳ�������
A=figure(1);
hold on;
%suptitle(animal_name);
%ha = tight_subplot(1,2,[1 1],[1 1],[1 1]);
set(gcf, 'position', [0 0 1200 400]);%����fig�ߴ�
%suptitle(animal);
group={'saline', 'CNO'};
color={'k','r'};
type={'single case','population'};
click=[20 28 40 62 90 125];
curve=cell(length(group),1);%���ڴ洢ͼ����
%������ͬdelay��probe��CP
for n=1:length(group)
    [animal_name,dataChoiceResult]=fFindChoiceMat(strcat('D:\xulab\behavior\pyx014\probe',group{n}));
    nsample=size(dataChoiceResult,1)-1;
    curve{n}=cell(size(dataChoiceResult,1),1);
    choicep=fGetChoicep(dataChoiceResult,diff_delay_num,click);%choicep(nsample+1,6,diff_delay_num)
    for i=1:nsample
        curve{n}{i}=plotIndividualCurve(A,choicep(i,:,:),color{n},diff_delay_num,click);
    end
    %population curve
    curve{n}{nsample+1}=plotPopulationCurve(A,choicep(end,:,:),color{n},diff_delay_num,click);   
end

function [choicep] = fGetChoicep(dataChoice,diff_delay_num,click)
    
    choicep=zeros(length(dataChoice),6,diff_delay_num);%��һά��ʾÿ��session��ƽ�����ڶ�ά��ʾ6��Ƶ�ʵ���ȷ�ʣ�����ά��ʾ��ͬdelay����
    for i=1:length(dataChoice)
        %�ҳ���ͬdelay���ȣ�������delay���ܳɱ�ʶ��
          maxdelay=max(dataChoice{i}(:,4));
        mindelay=min(dataChoice{i}(:,4));
        part=(maxdelay-mindelay)/diff_delay_num;
        if diff_delay_num==3
            shortdelay=(dataChoice{i}(:,4)<=(mindelay+part));
            middledelay=(dataChoice{i}(:,4)>(mindelay+part)).*(dataChoice{i}(:,4)<(mindelay+2*part));
            longdelay=(dataChoice{i}(:,4)>=(mindelay+2*part));
            delay=[shortdelay,middledelay,longdelay];%���ܳɶ�ά���飬ÿһ�д���ͬ���ȵ�delay
            %delay=ones(length(shortdelay),3);%���Ժ��ܵ�plot�Ƿ��Ǻ�
        else
            shortdelay=(dataChoice{i}(:,4)<=(mindelay+part));
            longdelay=(dataChoice{i}(:,4)>=(mindelay+part));
            delay=[shortdelay,longdelay];%���ܳɶ�ά���飬ÿһ�д���ͬ���ȵ�delay
            %delay=ones(length(shortdelay),3);%���Ժ��ܵ�plot�Ƿ��Ǻ�
        end
        %�ҳ�ѡ��ĸ������ͣ������������ʾ���ʶ��  
    %     if clickright==20
    %         left=(dataChoice{i}(:,1)>50);
    %         right=(dataChoice{i}(:,1)<50);
    %     else
    %         left=(dataChoice{i}(:,1)<50);
    %         right=(dataChoice{i}(:,1)>50);
    %     end
    %     all=left+right;
        %side=[all,left,right];%���ܳɶ�ά���飬ÿһ�д���ͬ��
        correct=dataChoice{i}(:,2)==1;
        error=dataChoice{i}(:,2)==2;
        %miss=(dataChoice(:,2)==2);
        %violation=(dataChoice(:,2)==3);
        do=correct+error;
        %nomiss=do+violation;
        %�ҳ���ͬprobe
        probe=zeros(length(dataChoice{i}),6);
        for j=1:6
            probe(:,j)=(dataChoice{i}(:,1)==click(j));
        end
        %�������probe��delay��CP
        for j=1:6
            for k=1:diff_delay_num
                if j>3%ѡ���Ƶ��j>3ʱ������ȷ�ʣ���֮�Ǵ�����
                    choicep(i,j,k)=sum(double(correct.*probe(:,j).*delay(:,k)))/sum(double(do.*probe(:,j).*delay(:,k)));
                else
                    choicep(i,j,k)=sum(double(error.*probe(:,j).*delay(:,k)))/sum(double(do.*probe(:,j).*delay(:,k)));
                end
            end
        end
    end

end
function [curve]=plotPopulationCurve(A,choicep,color,diff_delay_num,click)
    if diff_delay_num==3    
        delay_name={'shortdelay','middledelay','longdelay'};
    else
        delay_name={'shortdelay','longdelay'};
    end
    %���
    toneOct = log2(click/click(1))';%ת�ó�������
    ffun=fittype('a+b./(1+exp(-k.*(x-c)))','independent','x');
    x=toneOct(1)-0.5:0.01:toneOct(end)+0.5;%�ڱ߽������ټ���һЩ��Χ
    yfit=cell(size(choicep,3));
    for j=1:size(choicep,3)
        f=fit(toneOct,choicep(1,:,j)',ffun,'Startpoint',[0,1,1,1]);%choicep(i,:,j)'ת�ó�������
        yfit{j}=f(x);
    end
    %��ͼ
    figure(A);
    for k=1:diff_delay_num
        subplot(1,diff_delay_num,k);
    %����ÿ��������Լ�ƽ������
    %     for j=1:size(yfit,1)-1
    %         scatter(toneOct,choicep(j,:,k),10,'b');
    %         hold on;
    %         plot(x,yfit{j,k},'color',[j/size(yfit,1),j/size(yfit,1),j/size(yfit,1)],'LineWidth',1);%ʵ������Ƚϣ���ɫ�𽥼���
    %         hold on;  
    %     end
    %     scatter(toneOct,choicep(end,:,k),50,'b','filled');
    %     hold on;
    %     plot(x,yfit{end,k},'k','LineWidth',3);
    %     hold on;
    %����Ctrl��exp������жԱ�
        scatter(toneOct,choicep(1,:,k),50,color,'filled');%ֻ��mean��ctrl 'k', exp 'r'
        hold on;
        curve=plot(x,yfit{k},color,'LineWidth',2);%ֻ��mean��ctrl 'k', exp 'r'
        hold on;
        xlabel(delay_name(k),'FontName','Arial','FontSize',14);
        set(gca, 'YLim',[0 1]);
        set(gca, 'XLim',[x(1) x(end)]); 
        plot([x(1),x(end)],[0.5,0.5],'--k','LineWidth',1);
        hold on;
        plot([1.3,1.3],[0,1],'--k','LineWidth',1);
        hold on;
    end
    %legend([curve_ctrl,curve_exp],'ctrl','bilateral infusion CNO','Location','Best');
end
function [curve]=plotIndividualCurve(A,choicep,color_population,diff_delay_num,click)
    if diff_delay_num==3    
        delay_name={'shortdelay','middledelay','longdelay'};
    else
        delay_name={'shortdelay','longdelay'};
    end   
    %���
    toneOct = log2(click/click(1))';%ת�ó�������
    ffun=fittype('a+b./(1+exp(-k.*(x-c)))','independent','x');
    x=toneOct(1)-0.5:0.01:toneOct(end)+0.5;%�ڱ߽������ټ���һЩ��Χ
    yfit=cell(size(choicep,3));
    for j=1:size(choicep,3)
        f=fit(toneOct,choicep(1,:,j)',ffun,'Startpoint',[0,1,1,1]);%choicep(i,:,j)'ת�ó�������
        yfit{j}=f(x);
    end
    figure(A);
    %��ͼ
    colormode={[0.7,0.7,0.7],[1,0.7,0.7]};
    colorp=strfind('kr',color_population);
    color=colormode{colorp};
    for k=1:diff_delay_num
        subplot(1,diff_delay_num,k);
    %����ÿ��������Լ�ƽ������
    %     for j=1:size(yfit,1)-1
    %         scatter(toneOct,choicep(j,:,k),10,'b');
    %         hold on;
    %         plot(x,yfit{j,k},'color',[j/size(yfit,1),j/size(yfit,1),j/size(yfit,1)],'LineWidth',1);%ʵ������Ƚϣ���ɫ�𽥼���
    %         hold on;  
    %     end
    %     scatter(toneOct,choicep(end,:,k),50,'b','filled');
    %     hold on;
    %     plot(x,yfit{end,k},'k','LineWidth',3);
    %     hold on;
    %����Ctrl��exp������жԱ�
        curve=plot(x,yfit{k},'Color',color,'LineWidth',2);%ֻ��mean��ctrl 'k', exp 'r'
        hold on;
        xlabel(delay_name(k),'FontName','Arial','FontSize',14);
        set(gca, 'YLim',[0 1]);
        set(gca, 'XLim',[x(1) x(end)]); 
        plot([x(1),x(end)],[0.5,0.5],'--k','LineWidth',1);
        hold on;
        plot([1.3,1.3],[0,1],'--k','LineWidth',1);
        hold on;
    end
end