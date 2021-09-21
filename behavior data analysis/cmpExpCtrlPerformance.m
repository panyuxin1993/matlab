dirE='E:\xulab\behavior\CNO experiment\infusion\bilateral\';
dirC='E:\xulab\behavior\CNO experiment\CNO ctrl\infusion\bilateral\';
dirCip='E:\xulab\behavior\CNO experiment\CNO ctrl\ip\';
dirEip='E:\xulab\behavior\CNO experiment\ip\';
group={'saline', 'CNO'};
color={'k','r'}; 
type=  {'single case','population'};
curve=cell(length(group),1);%用于存储图像句柄  
transform_method={};

%%二维点图，比较correct rate, violation rate, miss rate, curve slope，reaction time from go cue(for correct trials),区分长短delay(3*5 subplot)
%paraB=cell(5,3,length(group));%第一维代表五种指标cor,vio,miss,RT,cor_diff(endtrial-probe)，第二维代表总、短delay、长delay，第三维代表两种实验条件（saline,CNO)
% if strcmp(transform_method,'uni')
%     %group={'ipsi saline', 'ipsi CNO', 'contra saline', 'contra CNO'};
%     group={'saline', 'CNO', 'saline', 'CNO'};  
% end

%    定义delay组数或者自定义delay分组标准获得不同delay长度的参数paraB
ndelaygroup=2;
delaysetting=[300 900 900 1500];
ndifficultygroup=2;
[paraE,pparaE,paraEdiff,pparaEdiff,delaybylength]=fGetpara(dirE, transform_method, ndelaygroup ,group, ndifficultygroup);
[paraC,pparaC,paraCdiff,pparaCdiff,delaybylength]=fGetpara(dirC, transform_method, ndelaygroup ,group, ndifficultygroup);
Be=fPlotOverallPerformanceScatter(paraE,pparaE,delaybylength);%二维点图，比较correct rate, violation rate, miss rate, curve slope，reaction time from go cue(for correct trials),区分长短delay，或者难度(3*5 subplot)
Bc=fPlotOverallPerformanceScatter(paraC,pparaC,delaybylength);%for control
C=fCmpOverallPerformance(paraE,paraEdiff,paraC,paraCdiff,ndelaygroup);%画CNO-saline的差值，比较不同情况下的区别
figure(Be);
suptitle('hM4D');
figure(Bc);
suptitle('mCherry');

function [fig] = fCmpOverallPerformance(paraE,paraEdiff,paraC,paraCdiff,ndelaygroup)
    fig=figure;%2*4张图，分别存correct rate, violation rate, miss rate, reaction time；第一行总结果，  第二行区分delay，不只是长短，而是分成200ms一组的6组,第三行区分difficulty
    set(gcf, 'position', [0 0 1400 800]);%控制fig尺寸
    ylabelstr={'Correct rate change(%)';'Violation rate change(%)';'Miss rate change(%)';'Reaction time change(ms)'};
    titlestr={'Correct';'Violation';'Miss';'Reaction time'};
    a=suptitle('CNO-saline');
    set(a,'FontSize',14);
    %首行合并各种指标散点图于一张大图
    subplot(3,4,1:3);%首行合并
    pararaw=zeros(length(paraE{1,1,1,1}),4);%第一维表示样本数，值表示CNO-saline情况，第二维是四种指标
    pararawC=zeros(length(paraC{1,1,1,1}),4);
    for i=1:4
        pararaw(:,i)=paraE{i,1,2,1}-paraE{i,1,1,1};
        pararawC(:,i)=paraC{i,1,2,1}-paraC{i,1,1,1};
    end
    for k=1:3
        yyaxis left;
        scatter(ones(size(pararaw,1),1)*k-0.2,pararaw(:,k),30,'filled');
        hold on;
        scatter(ones(size(pararawC,1),1)*k+0.2,pararawC(:,k),30); 
    end
    ytext=get(gca,'Ylim');
    for k=1:3
        [~,p]=ttest(pararaw(:,k));
        text(k-0.2,ytext(end),plabelsymbol(p),'FontName','Arial','FontSize',14);
        [~,p]=ttest(pararawC(:,k));
        text(k+0.2,ytext(end),plabelsymbol(p),'FontName','Arial','FontSize',14);
        [~,p2]=ttest2(pararaw(:,k),pararawC(:,k));
%         [~,p2]=ranksum(pararaw(:,k),pararawC(:,k));
        text(k,ytext(end)+0.2*(ytext(end)-ytext(1)),plabelsymbol(p2),'FontName','Arial','FontSize',14);
        plot([k-0.2,k+0.2],[ytext(end)+0.1*(ytext(end)-ytext(1)),ytext(end)+0.1*(ytext(end)-ytext(1))],'k-');
    end
    ylabel('\Delta (%)');
%     set(gca,'Ylim',[-40,60]);
    yyaxis right;
    scatter(ones(size(pararaw,1),1)*4-0.2,pararaw(:,4),30,'filled');%exp group
    hold on;
    scatter(ones(size(pararawC,1),1)*4+0.2,pararawC(:,4),30);%ctrl group,hollow
%     for i=1:size(pararaw,1)%用线连接相同session的参数
%         plot(pararaw(i,:));
%         hold on;
%     end
    [~,p]=ttest(pararaw(:,4));
    ytext=get(gca,'Ylim');
    text(4-0.2,ytext(end),plabelsymbol(p),'FontName','Arial','FontSize',14);
    [~,p]=ttest(pararawC(:,4));
    text(4+0.2,ytext(end),plabelsymbol(p),'FontName','Arial','FontSize',14);
    plot([0,size(pararaw,2)+1],[0,0],'--k','LineWidth',0.5);
    [~,p2]=ttest2(pararaw(:,4),pararawC(:,4));
%     [~,p2]=ranksum(pararaw(:,4),pararawC(:,4));
    text(4,ytext(end)+0.1*(ytext(end)-ytext(1)),plabelsymbol(p2),'FontName','Arial','FontSize',14);
    plot([4-0.2,4+0.2],[ytext(end)+0.05*(ytext(end)-ytext(1)),ytext(end)+0.05*(ytext(end)-ytext(1))],'k-');
    hold on;
%     set(gca,'Ylim',[-400,600]);
    ylabel('\Delta (ms)');
    set(gca,'XTick',[1 2 3 4]);
    set(gca,'XTickLabel',titlestr);
    set(gca,'FontName','Arial','FontSize',14);
    box off;
    
    x=[0.8,1.2,1.8,2.2];%各点的x坐标
    for i=1:4 %控制四列       
       subplot(3,4,4+i);%第二行,delay
       parabydelay=zeros(length(paraEdiff{i,1,1,1}),ndelaygroup);
       parabydelayC=zeros(length(paraCdiff{i,1,1,1}),ndelaygroup);
       for j=2:size(paraEdiff,2)
          parabydelay(:,j-1)=paraEdiff{i,j,1,end}';
          parabydelayC(:,j-1)=paraCdiff{i,j,1,end}';
       end
       for j=1:size(parabydelay,1)
           plot([x(1),x(2)],parabydelay(j,:),'r','LineWidth',0.5);
           hold on;
       end
       for j=1:size(parabydelayC,1)
           plot([x(3),x(4)],parabydelayC(j,:),'k','LineWidth',0.5);
       end
       fige1=scatter(x(1)*ones(size(parabydelay,1),1),parabydelay(:,1),20,'r','filled');
       fige2=scatter(x(2)*ones(size(parabydelay,1),1),parabydelay(:,2),20,'r');  
       figc1=scatter(x(3)*ones(size(parabydelayC,1),1),parabydelayC(:,1),20,'k','filled');
       figc2=scatter(x(4)*ones(size(parabydelayC,1),1),parabydelayC(:,2),20,'k');
       plot([0,3],[0,0],'--k','LineWidth',0.5);%0处水平线
       hold on;
       set(gca, 'XLim',[0 size(parabydelay,2)+1]);
%        set(gca, 'YLim',[0 100]);
       set(gca,'XTick',[1,2]);
%        set(gca,'XTickLabel',['1';'2']);
       set(gca,'XTickLabel',{'hM4D','mCherry'});
%        xlabel('delay');
       ylabel(ylabelstr{i});
%        title(titlestr{i});
       [~,p]=ttest(parabydelay(:,1),parabydelay(:,end));
       xtext=get(gca,'Xlim');
       ytext=get(gca,'Ylim');
       text(1,ytext(end),plabelsymbol(p),'FontName','Arial','FontSize',14);
       [~,p]=ttest(parabydelayC(:,1),parabydelayC(:,end));
       text(2,ytext(end),plabelsymbol(p),'FontName','Arial','FontSize',14);
       parabydelaydiff=parabydelay(:,end)-parabydelay(:,1);
       parabydelaydiffC=parabydelayC(:,end)-parabydelayC(:,1);
       [~,p2]=ttest2(parabydelaydiff,parabydelaydiffC);
%        [~,p2]=ranksum(parabydelaydiff,parabydelaydiffC);
       text(1.5,ytext(end)+0.2*(ytext(end)-ytext(1)),plabelsymbol(p2),'FontName','Arial','FontSize',14);
       plot([1,2],[ytext(end)+0.1*(ytext(end)-ytext(1)),ytext(end)+0.1*(ytext(end)-ytext(1))],'k-');
       if i==4
       legend([fige1,fige2],'short delay','long delay');
       legend([figc1,figc2],'short delay','long delay');
       end
       set(gca,'FontName','Arial','FontSize',14);
       box off;
       
       subplot(3,4,8+i);%第三行,difficulty
       easy=paraEdiff{i,1,1,2};
       hard=paraEdiff{i,1,1,3};
       parabydifficulty=[easy' hard'];
       parabydifficultyC=[paraCdiff{i,1,1,2}' paraCdiff{i,1,1,3}'];
       for j=1:size(parabydifficulty,1)
           plot([x(1) x(2)],parabydifficulty(j,:),'r','LineWidth',0.5);
           hold on;
       end
       for j=1:size(parabydifficultyC,1)
           plot([x(3) x(4)],parabydifficultyC(j,:),'k','LineWidth',0.5);
           hold on;
       end

       fige1=scatter(ones(size(parabydifficulty,1),1)*x(1),parabydifficulty(:,1),20,'r','filled');
       fige2=scatter(ones(size(parabydifficulty,1),1)*x(2),parabydifficulty(:,2),20,'r');
       figc1=scatter(ones(size(parabydifficultyC,1),1)*x(3),parabydifficultyC(:,1),20,'k','filled');
       figc2=scatter(ones(size(parabydifficultyC,1),1)*x(4),parabydifficultyC(:,2),20,'k');
       plot([0,3],[0,0],'--k','LineWidth',0.5);
       hold on;
%        set(gca, 'XLim',[0 size(parabydifficulty,2)+1]);
%        set(gca, 'YLim',[0 100]);
       set(gca,'XTick',[1,2]);
%        set(gca,'XTickLabel',['1';'2']);
       set(gca,'XTickLabel',{'hM4D','mCherry'});
%        xlabel('difficulty');
       ylabel(ylabelstr{i});
%        title(titlestr{i});
       [~,p]=ttest(parabydifficulty(:,1),parabydifficulty(:,2));
%        xtext=get(gca,'Xlim');
       ytext=get(gca,'Ylim');
       text(1,ytext(end),plabelsymbol(p),'FontName','Arial','FontSize',14);
       [~,p]=ttest(parabydifficultyC(:,1),parabydifficultyC(:,2));
       text(2,ytext(end),plabelsymbol(p),'FontName','Arial','FontSize',14);
       parabydifficultydiff=parabydifficulty(:,2)-parabydifficulty(:,1);
       parabydifficultydiffC=parabydifficultyC(:,2)-parabydifficultyC(:,1);
       [~,p2]=ttest2(parabydifficultydiff,parabydifficultydiffC);
%        [~,p2]=ranksum(parabydifficultydiff,parabydifficultydiffC);
       text(1.5,ytext(end)+0.2*(ytext(end)-ytext(1)),plabelsymbol(p2),'FontName','Arial','FontSize',14);
       plot([1,2],[ytext(end)+0.1*(ytext(end)-ytext(1)),ytext(end)+0.1*(ytext(end)-ytext(1))],'k-');
       if i==4
       legend([fige1,fige2],'easy','hard');
       legend([figc1,figc2],'easy','hard');
       end
       set(gca,'FontName','Arial','FontSize',14);
       box off;
    end
end
function [fig] = fPlotOverallPerformanceScatter(paraE,pparaE,delaybylength)
fig=figure;
    set(gcf, 'position', [0 0 1500 810]);%控制fig尺寸
%     ylabelstr={' all trials',' short delay trials',' long delay trials'};
    ylabelstr={'CNO','CNO','CNO'};
    %xlabelstr={};
    titlestr={'Correct rate','Violation rate','Miss rate','Reaction time','Licking Consistency'};
    for i=1:3
        for j=1:4
            subplot(3,5,5*i+j-5);%五列分别画正确率，violation ,miss,RT,LCI

            if j==4%画对角线,RT
                plot([0,1500],[0,1500],'k','LineWidth',2);
                hold on;
            elseif j==5 %LCI
                plot([0,1],[0,1],'k','LineWidth',2);
                hold on;
            else
                plot([0,100],[0,100],'k','LineWidth',2);
                hold on;
            end
            ytext=get(gca,'Ylim');%获取坐标范围，在特定处显示text
            if size(paraE, 3)==2%bilateral
                if i==1%第一行画所有trial的总体情况
                    curve=scatter(paraE{j,i,1,1},paraE{j,i,2,1},50,'k','filled');
                    hold on;                    
                    textstr=plabel(pparaE(j,i,1,1));%自动根据p值选择符号
                elseif i==2%第二行画区分delay长度
%                     if j~=5%第五列改动为区分长短delay和hard一共四种情况
                    curve_short=scatter(paraE{j,2,1,end},paraE{j,2,2,end},50,'k','filled');%最短delay
                    hold on;
                    curve_long=scatter(paraE{j,end,1,end},paraE{j,end,2,end},50,'b','filled');%最长delay
                    hold on;
%                     textstr={strcat('short:',plabelsymbol(pparaE(j,2,1,1)));strcat('long: ',plabelsymbol(pparaE(j,3,1,1)));strcat('short vs long: ',plabelsymbol(pparaEdiff(i-1,j)))};
                   textstr={strcat('short:',plabelsymbol(pparaE(j,2,1,end)));strcat('long: ',plabelsymbol(pparaE(j,end,1,end)))};
%                     else
%                         curve_short_difficult=scatter(paraE{1,2,1,3},paraE{1,2,2,3},30,'k');
%                         hold on;
%                         curve_long_difficult=scatter(paraE{1,3,1,3},paraE{1,3,2,3},30,'b');
%                         hold on;
%                         curve_short_easy=scatter(paraE{1,2,1,2},paraE{1,2,2,2},30,'k','filled');
%                         hold on;
%                         curve_long_easy=scatter(paraE{1,3,1,2},paraE{1,3,2,2},30,'b','filled');
%                         hold on;
% 
%                         textstr={};%strcat('short p=',num2str(pparaE(j,2,1,1)));strcat('long p=',num2str(pparaE(j,3,1,1)));strcat('long p=',num2str(pparaE(j,3,1,1)));strcat('long p=',num2str(pparaE(j,3,1,1)))};
%                     end
                    if j==4%仅需在之后一列画legend
                        legend([curve_short curve_long],strcat('shortest delay(',num2str(delaybylength(1)),'-',num2str(delaybylength(2)),'ms)'),strcat('longest delay(',num2str(delaybylength(end-1)),'-',num2str(delaybylength(end)),'ms)'),'Location','Best');
%                           legend([curve_short_easy curve_long_easy curve_short_difficult curve_long_difficult],strcat('short delay(300-900ms) easy p=',num2str(pparaE(1,2,1,2))),strcat('long delay(900-1500ms) easy p=',num2str(pparaE(1,3,1,2))),strcat('short delay(300-900ms) hard p=',num2str(pparaE(1,2,1,3))),strcat('long delay(900-1500ms) hard p=',num2str(pparaE(1,3,1,3))),'Location','Best');
                    end
                else%第三行画区分难度end/probe
                    curve_end=scatter(paraE{j,1,1,2},paraE{j,1,2,2},50,'k','filled');%end
                    hold on;
                    curve_probe=scatter(paraE{j,1,1,3},paraE{j,1,2,3},50,'k');%probe
                    hold on;
%                     textstr={strcat('easy:',plabelsymbol(pparaE(j,1,1,2)));strcat('hard:',plabelsymbol(pparaE(j,1,1,3)));strcat('easy vs hard:',plabelsymbol(pparaEdiff(i-1,j)))};
                    textstr={strcat('easy:',plabelsymbol(pparaE(j,1,1,2)));strcat('hard:',plabelsymbol(pparaE(j,1,1,end)))};
                    if j==4%仅需在之后一列画legend
                        legend([curve_end curve_probe],'easy','hard','Location','Best');
                    end
                end
                text(0,ytext(end),textstr,'FontSize',12);
            else %size(paraE, 3)==4；unilateral data
                if i==1%第一行画所有trial的总体情况
                    curve_ipsi=scatter(paraE{j,i,1,1},paraE{j,i,2,1},50,'k','filled');
                    hold on;
                    curve_contra=scatter(paraE{j,i,3,1},paraE{j,i,4,1},50,'k','d');
                    hold on;
                    textstr={strcat('ipsi',plabelsymbol(pparaE(j,i,1,1)));strcat('contra',plabelsymbol(pparaE(j,i,2,1)))};
                    if j==4%仅需在之后一列画legend
                        legend([curve_ipsi curve_contra],'ipsi','contra','Location','Best');
                    end
                elseif i==2%第二行画区分delay长度
                    curve_short_ipsi=scatter(paraE{j,2,1,1},paraE{j,2,2,1},50,'k','filled');%短delay
                    hold on;
                    curve_long_ipsi=scatter(paraE{j,3,1,1},paraE{j,3,2,1},50,'b','filled');%长delay
                    hold on;                 
                    curve_short_contra=scatter(paraE{j,2,1,1},paraE{j,2,2,1},50,'k','filled','d');%短delay,contra
                    hold on;
                    curve_long_contra=scatter(paraE{j,3,1,1},paraE{j,3,2,1},50,'b','filled','d');%长delay,contra
                    hold on;
                    textstr={strcat('ipsi short p=',plabelsymbol(pparaE(j,2,1,1)));strcat('ipsi long p=',plabelsymbol(pparaE(j,3,1,1)));strcat('contra short p=',plabelsymbol(pparaE(j,2,2,1)));strcat('contra long p=',plabelsymbol(pparaE(j,3,2,1)))};
                    if j==4%仅需在之后一列画legend
                        legend([curve_short_ipsi curve_long_ipsi curve_short_contra curve_long_contra],'ipsi short delay','ipsi long delay','contra short delay','contra long delay','Location','Best');
                    end
                else%第三行画区分难度end/probe
                    curve_end_ipsi=scatter(paraE{j,1,1,2},paraE{j,1,2,2},50,'k','filled');%end
                    hold on;
                    curve_probe_ipsi=scatter(paraE{j,1,1,3},paraE{j,1,2,3},50,'k');%probe
                    hold on;
                    curve_end_contra=scatter(paraE{j,1,1,2},paraE{j,1,2,2},50,'k','filled','d');%end,contra
                    hold on;
                    curve_probe_contra=scatter(paraE{j,1,1,3},paraE{j,1,2,3},50,'k','d');%probe,contra
                    hold on;
                    textstr={strcat('ipsi easy p=',plabelsymbol(pparaE(j,1,1,2)));strcat('ipsi difficult p=',plabelsymbol(pparaE(j,1,1,3)));strcat('contra easy p=',plabelsymbol(pparaE(j,1,2,2)));strcat('contra difficult p=',plabelsymbol(pparaE(j,1,2,3)))};
                    if j==4%仅需在之后一列画legend
                        legend([curve_end_ipsi curve_probe_ipsi curve_end_contra curve_probe_contra],'ipsi easy','ipsi difficult','contra easy','contra difficult','Location','Best');
                    end
                end
                text(0,ytext(end),textstr,'FontSize',12);
            end

            % set(gca, 'YLim',[0 2000]);
            % set(gca,'XLim',[0,2000]);
%             if j==1
                ylabel(ylabelstr{i},'FontName','Arial','FontSize',14);
%             end
%             if i==3
                %xlabel(strcat('Saline',xlabelstr{i}),'FontName','Arial','FontSize',14);
                if j==4
                    xlabel('Saline(ms)','FontName','Arial','FontSize',14);
                else
                    xlabel('Saline(%)','FontName','Arial','FontSize',14);
                end
%             end
            if i==1
                title(titlestr{j},'FontName','Arial','FontSize',14);
            end
            hold on;
            set(gca,'FontName','Arial','FontSize',14);
            box off;%去掉边框
        end
    end
end
function [paraB,pparaB,paraBdiff, pparaBdiff,delaybylength]=fGetpara(dir, transform_method, vardelay,treatmentgroup, ndifficultygroup)
% paraB=cell(5,ndelaygroup+1,length(group),3);
%第一维代表五种指标（cor rate, violation rate, miss rate, RT, licking consistency index LCI)，
%第二维代表总、短、长delay，根据ndelaygroup确定组数;若vardelay为自定义划分标准，则默认分成两组
%第三维代表两种实验条件（saline,CNO)，
%第四维是难度（end/probe)
% pparaB=ones(size(paraB,1),size(paraB,2),size(paraB,3)/2,size(paraB,4));%保存各个paraB中的参数CNO-saline的p值
% paraBdiff=cell(5,ndelaygroup+1,length(group)/2,3);%区分CNO-saline的效果
% pparaBdiff=ones(2,size(paraB,1));%保存paraB中的参数CNO-saline的值在区分delay/difficulty时的p值,只为作图服务，索引对应fig位置,第一行不分组，不需比较，故此处第一维代表分delay，第二维代表分difficulty
% delayPart 根据ndelaygroup确定的每个部分的delay长度
if length(vardelay)==1
    ndelaygroup=vardelay;
elseif length(vardelay)==4
    ndelaygroup=2;
    delaysetting=vardelay;
end
if strcmp(transform_method,'uni')
    paraB=cell(5,ndelaygroup+1,2*length(treatmentgroup),ndifficultygroup+1);%第三维区分ipsi, contra来看
    paraBdiff=cell(5,ndelaygroup+1,length(treatmentgroup),ndifficultygroup+1);%区分CNO-saline的效果
else
    paraB=cell(5,ndelaygroup+1,length(treatmentgroup),ndifficultygroup+1);%第一维代表五种指标，第二维代表总、短、长delay，第三维代表两种实验条件（saline,CNO)，第四维是难度（end/probe)
    paraBdiff=cell(5,ndelaygroup+1,length(treatmentgroup)/2,ndifficultygroup+1);
end
for n=1:length(treatmentgroup)
    [animal_name,dataChoiceResult]=fFindChoiceMat(strcat(dir,treatmentgroup{n}));
    nsample=size(dataChoiceResult,1)-1;
    for i=1:nsample
        %delay先确定好
        if length(vardelay)==1
            [delaybylength,delay]=fDelayGroup(dataChoiceResult{i}(:,4),'ngroup',ndelaygroup);%汇总成二维数组，每一列代表不同长度的delay,% logical 相加就变成double
        elseif length(vardelay)==4
            [delaybylength,delay]=fDelayGroup(dataChoiceResult{i}(:,4),'user-defined',delaysetting);
        end
        %获取trial results
        correct=dataChoiceResult{i}(:,2)==1;%取值1-4分别代表cor,err,miss.err
        error=dataChoiceResult{i}(:,2)==2;
        miss=dataChoiceResult{i}(:,2)==3;
        violation=dataChoiceResult{i}(:,2)==4;
        novio=correct+error;
        do=novio+violation;
        %获取endtrial/probetrial
        click=unique(dataChoiceResult{i}(:,1));
        sort(click);
        pend=logical((dataChoiceResult{i}(:,1)==click(1))+(dataChoiceResult{i}(:,1)==click(end)));
        pprobe=logical(ones(size(dataChoiceResult{i},1),1)-pend);
        difficulty=[pend+pprobe,pend,pprobe];%汇总成二维数组，每一列代表不同难度
        %获取reaction time
        indRT=double((dataChoiceResult{i}(:,5)>0).*(dataChoiceResult{i}(:,2)==1));%只看正确的trial
        indLCI=double(dataChoiceResult{i}(:,2)==1);%只看正确的trial
%         indLCI=double(dataChoiceResult{i}(:,2)==2);%只看错误的trial
%         indLCI=double(dataChoiceResult{i}(:,2)~=3);%看所有trial，除去不舔水的miss
        for k=1:size(difficulty,2)%区分难度
            %获取ipsi/contra
            if strcmp(transform_method,'uni')
                dataChoiceResult{i}(:,1)=dataChoiceResult{i}(:,3);%左右调换的话就用第六列的声音对应的左右替换第一列频率
                protate=logical((dataChoiceResult{i}(:,6)==2)+(dataChoiceResult{i}(:,6)==3));%这些trials需要旋转，也就是上下，左右各颠倒一次
                dataChoiceResult{i}(protate,1)=7-dataChoiceResult{i}(protate,1);
                ipsi=(dataChoiceResult{i}(:,1)<=3);
                contra=(dataChoiceResult{i}(:,1)>=3);
                for j=1:size(delay,2)%各种长短的delay
                    paraB{1,j,n,k}(i)=100*sum(correct.*delay(:,j).*ipsi.*difficulty(:,k))/sum(novio.*delay(:,j).*ipsi.*difficulty(:,k));
                    paraB{2,j,n,k}(i)=100*sum(violation.*delay(:,j).*ipsi.*difficulty(:,k))/sum(do.*delay(:,j).*ipsi.*difficulty(:,k));
                    paraB{3,j,n,k}(i)=100*sum(miss.*delay(:,j).*ipsi.*difficulty(:,k))/sum(delay(:,j).*ipsi.*difficulty(:,k));
                    paraB{4,j,n,k}(i)=sum(indRT.*dataChoiceResult{i}(:,5).*delay(:,j).*ipsi.*difficulty(:,k))/sum(indRT.*delay(:,j).*ipsi.*difficulty(:,k));
                    paraB{5,j,n,k}(i)=sum(indLCI.*dataChoiceResult{i}(:,7).*delay(:,j).*ipsi.*difficulty(:,k))/sum(indLCI.*delay(:,j).*ipsi.*difficulty(:,k));
                    paraB{1,j,n+2,k}(i)=100*sum(correct.*delay(:,j).*contra.*difficulty(:,k))/sum(novio.*delay(:,j).*contra.*difficulty(:,k));
                    paraB{2,j,n+2,k}(i)=100*sum(violation.*delay(:,j).*contra.*difficulty(:,k))/sum(do.*delay(:,j).*contra.*difficulty(:,k));
                    paraB{3,j,n+2,k}(i)=100*sum(miss.*delay(:,j).*contra.*difficulty(:,k))/sum(delay(:,j).*contra.*difficulty(:,k));
                    paraB{4,j,n+2,k}(i)=sum(indRT.*dataChoiceResult{i}(:,5).*delay(:,j).*contra.*difficulty(:,k))/sum(indRT.*delay(:,j).*contra.*difficulty(:,k));
                    paraB{5,j,n+2,k}(i)=sum(indLCI.*dataChoiceResult{i}(:,7).*delay(:,j).*contra.*difficulty(:,k))/sum(indLCI.*delay(:,j).*contra.*difficulty(:,k));
                end
            else
                for j=1:size(delay,2)%各种长短的delay
                    paraB{1,j,n,k}(i)=100*sum(correct.*delay(:,j).*difficulty(:,k))/sum(novio.*delay(:,j).*difficulty(:,k));
                    paraB{2,j,n,k}(i)=100*sum(violation.*delay(:,j).*difficulty(:,k))/sum(do.*delay(:,j).*difficulty(:,k));
                    paraB{3,j,n,k}(i)=100*sum(miss.*delay(:,j).*difficulty(:,k))/sum(delay(:,j).*difficulty(:,k));
                    paraB{4,j,n,k}(i)=sum(indRT.*dataChoiceResult{i}(:,5).*delay(:,j).*difficulty(:,k))/sum(indRT.*delay(:,j).*difficulty(:,k));
                    paraB{5,j,n,k}(i)=sum(indLCI.*dataChoiceResult{i}(:,7).*delay(:,j).*difficulty(:,k))/sum(indLCI.*delay(:,j).*difficulty(:,k));
                end
            end
        end
    end
end

pparaB=ones(size(paraB,1),size(paraB,2),size(paraB,3)/2,size(paraB,4));%保存各个paraB中的参数CNO-saline的p值
for i1=1:size(pparaB,1)%parameter
    for i2=1:size(pparaB,2)%delay
        for i4=1:size(pparaB,4)%difficulty
            paraBdiff{i1,i2,1,i4}=paraB{i1,i2,2,i4}-paraB{i1,i2,1,i4};%CNO-saline
            if size(paraB, 3)==2
                [~,p]=ttest(paraB{i1,i2,1,i4},paraB{i1,i2,2,i4});
                %                     p=round(p,4);%保留小数点后4位
                pparaB(i1,i2,1,i4)=p;
            else%size(paraB, 3)==4
                [~,pipsi]=ttest(paraB{i1,i2,1,i4},paraB{i1,i2,2,i4});
                pparaB(i1,i2,1,i4)=pipsi;
                [~,pcontra]=ttest(paraB{i1,i2,1,i4},paraB{i1,i2,2,i4});
                pparaB(i1,i2,2,i4)=pcontra;
                paraBdiff{i1,i2,3,i4}=paraB{i1,i2,4,i4}-paraB{i1,i2,3,i4};
            end
        end
    end
end
%     compare difference CNO-saline
pparaBdiff=ones(2,size(paraB,1));%保存paraB中的参数CNO-saline的值在区分delay/difficulty时的p值,只为作图服务，索引对应fig位置,第一行不分组，不需比较，故此处第一维代表分delay，第二维代表分difficulty
for i=1:size(paraB,1)
    [~,pbydelay]=ttest(paraBdiff{i,2,1,1},paraBdiff{i,end,1,1});
    pparaBdiff(1,i)=pbydelay;
    [~,pbydifficulty]=ttest(paraBdiff{i,1,1,2},paraBdiff{i,1,1,end});
    pparaBdiff(2,i)=pbydifficulty;
end
end

function [delaybylength,delaygroup]= fDelayGroup(delayraw,varargin)
if nargin==3
    if strcmp(varargin{1},'ngroup')
        ngroup=varargin{2};
        %divide the delay into ngroup groups according to delay length
        mindelay=min(delayraw);
        maxdelay=max(delayraw);
        part=(maxdelay-mindelay)/ngroup;
        delaybylength=mindelay:part:maxdelay;
        delaygroup=ones(length(delayraw),1);
        for i=1:ngroup
            delaycol=(mindelay+(i-1)*part<=delayraw).*(mindelay+i*part>delayraw);
            delaygroup=[delaygroup delaycol];
        end
    elseif strcmp(varargin{1},'user-defined')
        delaybylength=varargin{2};
        if length(delaybylength)~=4
            warning('Input should be a 1*4 double vector');
        else
            delaygroup=ones(length(delayraw),1);
            for i=1:2
                delaycol=(delaybylength(2*i-1)<=delayraw).*(delaybylength(2*i)>=delayraw);
                delaygroup=[delaygroup delaycol];
            end
        end
    else
        warning('Input error for function fDelayGroup');
    end
else
    warning('Input number error for function fDelayGroup');
end
end
function [str]=plabel(pvalue)
if pvalue<0.05 && pvalue>=0.01
    str='p<0.05';
elseif pvalue<0.01 && pvalue>=0.001
    str='p<0.01';
elseif pvalue<0.001
    str='p<0.001';
else
    pvalue=round(pvalue,2);%保留两位即可
    str=strcat('p=',num2str(pvalue));
end
end
function [str]=plabelsymbol(pvalue)
if pvalue<0.05 && pvalue>=0.01
    str=' *';
elseif pvalue<0.01 && pvalue>=0.001
    str=' **';
elseif pvalue<0.001
    str=' ***';
else
    str='n.s.';
%     pvalue=round(pvalue,2);%保留两位即可
%     str=strcat('p=',num2str(pvalue));
end
end