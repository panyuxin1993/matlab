close all;
[num,txt,raw] =xlsread('D:\xulab\project\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
savepath='H:\2P\summary';
clear TAUC_combine;
trialTypeStr='cor and err';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
AUCtype='sensory';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
%% calculate AUC
%{
T=cell2table(raw(2:end,1:11));
T.Properties.VariableNames=strrep(raw(1,1:11),' ','_');%table variable name can't have ' ',so replace them
ind_session=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control');
ind_session=find(ind_session);
n_session=length(ind_session);
animal_unique=unique(T.animal(ind_session));
n_animal=length(animal_unique);

for i_session=1:n_session
    indrow=ind_session(i_session);
    TAUC_currentSession = fGetEpochAUCtableASession(T.file_path{indrow},T.animal{indrow},T.date{indrow},T.field{indrow},T.cell_type{indrow},trialTypeStr,AUCtype);
    if exist('TAUC_combine','var') && strcmp(AUCtype,'stimuli') %calculate choice probability
        if size(TAUC_currentSession.ITI,2)==2
            TAUC_combine=vertcat(TAUC_combine,TAUC_currentSession);
        elseif size(TAUC_currentSession.ITI,2)>2 %for sessions with probe sound, choose only the end sound
            TAUC_currentSession_formatted=TAUC_currentSession;
            for ivar=1:size(TAUC_currentSession_formatted,2)
                temp=table2array(TAUC_currentSession_formatted(:,ivar));
                if size(temp,2)>2
                    tempvarname=TAUC_currentSession_formatted.Properties.VariableNames{ivar};
                    TAUC_currentSession_formatted=removevars(TAUC_currentSession_formatted,TAUC_currentSession_formatted.Properties.VariableNames{ivar});
                    TAUC_currentSession_formatted = addvars(TAUC_currentSession_formatted,temp(:,[1,end]),'NewVariableNames',tempvarname,'Before',TAUC_currentSession_formatted.Properties.VariableNames{ivar});
                end
            end
            TAUC_combine=vertcat(TAUC_combine,TAUC_currentSession_formatted);
        end
    elseif exist('TAUC_combine','var') && ~strcmp(AUCtype,'stimuli')
        TAUC_combine=vertcat(TAUC_combine,TAUC_currentSession);
    else
        TAUC_combine=TAUC_currentSession;
    end
end
save([savepath,filesep,'trialType',trialTypeStr,'-',AUCtype,'TepochAUC.mat'],'TAUC_combine');
%}
%% plot histogram for choice AUC 
%
trialTypeStr='cor and err';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
AUCtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
load([savepath,filesep,'trialType',trialTypeStr,'-',AUCtype,'TepochAUC.mat']);
% load('F:\2P\summary\trialTypecorTepochAUC.mat');
celltype='syn';
manipulation='control';
%T_AUC=TAUC_combine(strcmp(TAUC_combine.celltype,'syn')|contains(TAUC_combine.celltype,'+other'),:);%syn, ..+other
T_AUC=TAUC_combine(strcmp(TAUC_combine.celltype,celltype),:);
n_animal=length(unique(T_AUC.animal));
psig=[0.025,0.975];
figAUC=figure;
set(figAUC,'Position',[100,100,800,200]);
subplot(1,5,1);
histogram(T_AUC.ITI,'FaceColor','w','BinWidth',0.05);
hold on;
histogram(T_AUC.ITI(T_AUC.pITI<psig(1)),'FaceColor','r','BinWidth',0.05);
histogram(T_AUC.ITI(T_AUC.pITI>psig(2)),'FaceColor','b','BinWidth',0.05);
ylabel('cell counts');
xlabel('AUC');
title('ITI');
set(gca,'Xlim',[0,1]);
set(gca,'FontSize',12);

box off;
subplot(1,5,2);
histogram(T_AUC.sound,'FaceColor','w','BinWidth',0.05);
hold on;
histogram(T_AUC.sound(T_AUC.psound<psig(1)),'FaceColor','r','BinWidth',0.05);
histogram(T_AUC.sound(T_AUC.psound>psig(2)),'FaceColor','b','BinWidth',0.05);
title('sound');
xlabel('AUC');
set(gca,'Xlim',[0,1]);
box off;
set(gca,'FontSize',12);
subplot(1,5,3);
histogram(T_AUC.delay,'FaceColor','w','BinWidth',0.05);
hold on;
histogram(T_AUC.delay(T_AUC.pdelay<psig(1)),'FaceColor','r','BinWidth',0.05);
histogram(T_AUC.delay(T_AUC.pdelay>psig(2)),'FaceColor','b','BinWidth',0.05);
title('delay');
xlabel('AUC');
set(gca,'Xlim',[0,1]);
box off;
set(gca,'FontSize',12);
subplot(1,5,4);
histogram(T_AUC.response,'FaceColor','w','BinWidth',0.05);
hold on;
histogram(T_AUC.response(T_AUC.presponse<psig(1)),'FaceColor','r','BinWidth',0.05);
histogram(T_AUC.response(T_AUC.presponse>psig(2)),'FaceColor','b','BinWidth',0.05);
title('response');
xlabel('AUC');
set(gca,'Xlim',[0,1]);
box off;
h=legend('n.s.','contra','ipsi','AutoUpdate','off');
set(h,'box','off');
set(gca,'FontSize',12);
subplot(1,5,5);
histogram(T_AUC.lick,'FaceColor','w','BinWidth',0.05);
hold on;
histogram(T_AUC.lick(T_AUC.plick<psig(1)),'FaceColor','r','BinWidth',0.05);
histogram(T_AUC.lick(T_AUC.plick>psig(2)),'FaceColor','b','BinWidth',0.05);
title('lick');
xlabel('AUC');
text(0.1,1,['n=',num2str(size(T_AUC,1)),'cells from ',num2str(n_animal),'animals'],'Unit','Normalized');
set(gca,'Xlim',[0,1]);
box off;
set(gca,'FontSize',12);
%indcate example as vertical line
Tsorted=sortrows(TAUC_combine,{'delay','response'},{'ascend','descend'});
indCaseRow=49;%shift from ipsi to contra
for i=1:5
    subplot(1,5,i);
    plot( table2array(Tsorted(indCaseRow,5+i))*ones(2,1),ylim,'k-');hold on;
    box off;
end
Tsorted=sortrows(TAUC_combine,{'delay','response'},{'descend','descend'});
indCaseRow=1;%shift from contra to contra
for i=1:5
    subplot(1,5,i);
    plot( table2array(Tsorted(indCaseRow,5+i))*ones(2,1),ylim,'r-');
    box off;
end
%save figure
saveas(figAUC,[savepath,filesep,celltype,'-',trialTypeStr,'-',AUCtype,'AUC at significance level of ',num2str(2*psig(1)),'.png'],'png');
set(figAUC,'PaperPosition',[0,0,8,2]);
saveas(figAUC,[savepath,filesep,celltype,'-',trialTypeStr,'-',AUCtype,'AUC at significance level of ',num2str(2*psig(1)),'.pdf'],'pdf');
%}
%% plot bar
%
n_cell= size(T_AUC,1);
psig=[0.025,0.975];
pAUCSigSC=zeros(3,4);%each row- contra, ipsi, n.s., each column- sound, delay, response, lick
pAUCSigSC(1,1)=sum(T_AUC.psound<psig(1))/n_cell;
pAUCSigSC(2,1)=sum(T_AUC.psound>psig(2))/n_cell;
pAUCSigSC(1,2)=sum(T_AUC.pdelay<psig(1))/n_cell;
pAUCSigSC(2,2)=sum(T_AUC.pdelay>psig(2))/n_cell;
pAUCSigSC(1,3)=sum(T_AUC.presponse<psig(1))/n_cell;
pAUCSigSC(2,3)=sum(T_AUC.presponse>psig(2))/n_cell;
pAUCSigSC(1,4)=sum(T_AUC.plick<psig(1))/n_cell;
pAUCSigSC(2,4)=sum(T_AUC.plick>psig(2))/n_cell;
pAUCSigSC(3,:)=1-sum(pAUCSigSC);
figBar=figure;
set(gcf,'Position',[100,100,400,160]);
subplot(1,2,2);
map=[1,0,0;0,0,1;0.5,0.5,0.5];%red-contra, ipsi-blue, n.s.-grey
colormap(map);%previous than R2017b
set(groot,'defaultAxesColorOrder',map);%after R2017b,and remember to remove the setting
bar(pAUCSigSC','stacked');

box off;
%set(gca,'XTick',1:4,'XTickLabel',{'sound','delay','response','lick'});
set(gca,'XTick',1:4,'XTickLabel',{'S','D','R','L'});
title('SC neurons');
xlim([0,5]);
set(gca,'FontSize',12);
SC_choice_AUC_005=table(n_cell*pAUCSigSC(1,:)',n_cell*pAUCSigSC(2,:)',n_cell*ones(4,1),'VariableNames',{'contra_sig','ipsi_sig','total_n'},'RowNames',{'sound';'delay';'response';'lick'});
save('H:\2P\summary\SC_choice_AUC_005.mat','SC_choice_AUC_005');

subplot(1,2,1);%for SC projecting M2 neurons
load('H:\2P\summary\M2_choice_AUC_005.mat');
pAUCSigM2=zeros(3,4);
pAUCSigM2(1,:)=M2_choice_AUC_005.contra_sig./M2_choice_AUC_005.total_n;
pAUCSigM2(2,:)=M2_choice_AUC_005.ipsi_sig./M2_choice_AUC_005.total_n;
pAUCSigM2(3,:)=1-sum(pAUCSigM2);
bar(pAUCSigM2','stacked');
box off;
%set(gca,'XTick',1:4,'XTickLabel',{'sound','delay','response','lick'});
set(gca,'XTick',1:4,'XTickLabel',{'S','D','R','L'});
title('SC-projecting M2 neurons');
set(gca,'FontSize',12);
xlim([0,5]);
h=legend('contra','ipsi','n.s.');
set(h,'box','off');
set(groot,'defaultAxesColorOrder','remove');
% set(gcf,'PaperPosition',[0,0,3,1.2]);
% saveas(figBar,[savepath,filesep,celltype,'-',trialTypeStr,'-',AUCtype,'-bar plot of AUC at significance level of ',num2str(2*psig(1)),'.pdf'],'pdf');
%}

%% choice probability
%choose cells with significant delay/sound choice AUC
trialTypeStr_datapool='cor and err';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
AUCtype_datapool='sensory';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
load([savepath,filesep,'trialType',trialTypeStr_datapool,'-',AUCtype_datapool,'TepochAUC.mat']);
celltype='syn';
T_AUCs=TAUC_combine(strcmp(TAUC_combine.celltype,celltype),:);
psig=[0,0.025;0.975,1];
TsigDataInfosensory = fGetSigDataInfoTable(T_AUCs, psig,AUCtype_datapool);

AUCtype_datapool='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
load([savepath,filesep,'trialType',trialTypeStr_datapool,'-',AUCtype_datapool,'TepochAUC.mat']);
celltype='syn';
T_AUCc=TAUC_combine(strcmp(TAUC_combine.celltype,celltype),:);
psig=[0,0.025;0.975,1];
TsigDataInfochoice = fGetSigDataInfoTable(T_AUCc, psig,AUCtype_datapool);

T_sig=join(TsigDataInfosensory,TsigDataInfochoice,'Keys',{'animal','date','nROI'});

trialTypeStr='cor';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
AUCtype='stimuli';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
load([savepath,filesep,'trialType',trialTypeStr,'-',AUCtype,'TepochAUC.mat']);
T_AUC=TAUC_combine(strcmp(TAUC_combine.celltype,celltype),:);
% T_AUC.ID=strcat(T_AUC.animal,'-',T_AUC.date,'-',int2str(T_AUC.nROI));

T_join=join(T_AUC,T_sig,'Keys',{'animal','date','nROI'});

[flag1,flag2]=fSigFlag(T_AUCs.sound,T_AUCs.psound,T_AUCc.sound,T_AUCc.psound,psig);
T_join=addvars(T_join,flag1,flag2,'NewVariableNames',{'flags1','flags2'});
[flag1,flag2]=fSigFlag(T_AUCs.delay,T_AUCs.pdelay,T_AUCc.delay,T_AUCc.pdelay,psig);
T_join=addvars(T_join,flag1,flag2,'NewVariableNames',{'flagd1','flagd2'});
[flag1,flag2]=fSigFlag(T_AUCs.response,T_AUCs.presponse,T_AUCc.response,T_AUCc.presponse,psig);
T_join=addvars(T_join,flag1,flag2,'NewVariableNames',{'flagr1','flagr2'});
[flag1,flag2]=fSigFlag(T_AUCs.lick,T_AUCs.plick,T_AUCc.lick,T_AUCc.plick,psig);
T_join=addvars(T_join,flag1,flag2,'NewVariableNames',{'flagl1','flagl2'});

figAUC=figure;
set(figAUC,'Position',[100,100,600,200]);
subplot(1,3,1);
% histogram(T_join.sound,'BinWidth',0.05);
hold on;
histogram([T_join.sound(T_join.flagCAUCsound>0,2);T_join.sound(T_join.flagCAUCsound<0,1)],'FaceColor','r','BinWidth',0.05);
histogram(T_join.sound(T_join.flagCAUCsound<0,1),'FaceColor','b','BinWidth',0.05);
% xlabel('AUC error v.s.correct for prefered stimulus')
ylabel('Cell number');
set(gca,'Xlim',[0,1]);
set(gca,'FontSize',12);
title('sound');

subplot(1,3,2);
% histogram(T_join.delay(:,2),'BinWidth',0.05);
hold on;
histogram([T_join.delay(T_join.flagCAUCdelay>0,2);T_join.delay(T_join.flagCAUCdelay<0,1)],'FaceColor','r','BinWidth',0.05);
histogram(T_join.delay(T_join.flagCAUCdelay<0,1),'FaceColor','b','BinWidth',0.05);
xlabel('AUC error v.s.correct for prefered stimulus')
ylabel('Cell number');
set(gca,'Xlim',[0,1]);
title('delay');
set(gca,'FontSize',12);

subplot(1,3,3);
% histogram(T_join.delay(:,2),'BinWidth',0.05);
hold on;
histogram([T_join.lick(T_join.flagCAUClick>0,2);T_join.lick(T_join.flagCAUClick<0,1)],'FaceColor','r','BinWidth',0.05);
histogram(T_join.lick(T_join.flagCAUClick<0,1),'FaceColor','b','BinWidth',0.05);
% xlabel('AUC error v.s.correct for prefered stimulus')
ylabel('Cell number');
set(gca,'Xlim',[0,1]);
set(gca,'FontSize',12);
% legend('n.s.','contra','ipsi');
h=legend('contra','ipsi');
set(h,'box','off');
title('lick');
set(figAUC,'PaperPosition',[0,0,6,2]);
saveas(figAUC,[savepath,filesep,celltype,'-',trialTypeStr,'-',AUCtype,'AUC error vs correct prefered stimulus.pdf'],'pdf');

figAUC2=figure;
set(figAUC2,'Position',[100,100,600,200]);
subplot(1,3,1);
% histogram(T_join.sound,'BinWidth',0.05);
hold on;
histogram([T_join.sound(T_join.flags2,2);T_join.sound(T_join.flags1,1)],'FaceColor','w','BinWidth',0.05);
histogram([T_join.sound(logical(T_join.flags2.*findSig(T_join.psound(:,2),psig).*(T_join.sound(:,2)>0.5)),2);T_join.sound(logical(T_join.flags1.*findSig(T_join.psound(:,1),psig).*(T_join.sound(:,1)>0.5)),1)],'FaceColor','r','BinWidth',0.05);
histogram([T_join.sound(logical(T_join.flags2.*findSig(T_join.psound(:,2),psig).*(T_join.sound(:,2)<0.5)),2);T_join.sound(logical(T_join.flags1.*findSig(T_join.psound(:,1),psig).*(T_join.sound(:,1)<0.5)),1)],'FaceColor','b','BinWidth',0.05);
ylabel('Cell number');
set(gca,'Xlim',[0,1]);
title('sound');
set(gca,'FontSize',12);

subplot(1,3,2);
% histogram(T_join.sound,'BinWidth',0.05);
hold on;
histogram([T_join.delay(T_join.flagd2,2);T_join.delay(T_join.flagd1,1)],'FaceColor','w','BinWidth',0.05);
histogram([T_join.delay(logical(T_join.flagd2.*findSig(T_join.pdelay(:,2),psig).*(T_join.delay(:,2)>0.5)),2);T_join.delay(logical(T_join.flagd1.*findSig(T_join.pdelay(:,1),psig).*(T_join.delay(:,1)>0.5)),1)],'FaceColor','r','BinWidth',0.05);
histogram([T_join.delay(logical(T_join.flagd2.*findSig(T_join.pdelay(:,2),psig).*(T_join.delay(:,2)<0.5)),2);T_join.delay(logical(T_join.flagd1.*findSig(T_join.pdelay(:,1),psig).*(T_join.delay(:,1)<0.5)),1)],'FaceColor','b','BinWidth',0.05);
ylabel('Cell number');
set(gca,'Xlim',[0,1]);
title('delay');
xlabel('AUC correct v.s.error for prefered stimulus');
set(gca,'FontSize',12);

subplot(1,3,3);
% histogram(T_join.sound,'BinWidth',0.05);
hold on;
histogram([T_join.lick(T_join.flagl2,2);T_join.lick(T_join.flagl1,1)],'FaceColor','w','BinWidth',0.05);
histogram([T_join.lick(logical(T_join.flagl2.*findSig(T_join.plick(:,2),psig).*(T_join.lick(:,2)>0.5)),2);T_join.lick(logical(T_join.flagl1.*findSig(T_join.plick(:,1),psig).*(T_join.lick(:,1)>0.5)),1)],'FaceColor','r','BinWidth',0.05);
histogram([T_join.lick(logical(T_join.flagl2.*findSig(T_join.plick(:,2),psig).*(T_join.lick(:,2)<0.5)),2);T_join.lick(logical(T_join.flagl1.*findSig(T_join.plick(:,1),psig).*(T_join.lick(:,1)<0.5)),1)],'FaceColor','b','BinWidth',0.05);
set(gca,'FontSize',12);

ylabel('Cell number');
set(gca,'Xlim',[0,1]);
title('lick');
h=legend('n.s.','correct','error');
set(h,'box','off');

set(figAUC2,'PaperPosition',[0,0,6,2]);
saveas(figAUC2,[savepath,filesep,celltype,'-',trialTypeStr,'-',AUCtype,'AUC cor vs err prefered stimulus significance p0.05.pdf'],'pdf');

figAUC3=figure;%M2 delay AUC
set(figAUC3,'Position',[100,100,200,200]);
load('D:\download\M2_delay_table.mat');
psig=[0,0.05];
[flag1,flag2]=fSigFlag(M2_delay_table.Sensory_auc,M2_delay_table.Sensory_p,M2_delay_table.Choice_auc,M2_delay_table.Choice_p,psig);
M2_delay_table=addvars(M2_delay_table,flag1,flag2,'NewVariableNames',{'flagd1','flagd2'});
histogram([M2_delay_table.CP_auc_contra(M2_delay_table.flagd2);M2_delay_table.CP_auc_ipsi(M2_delay_table.flagd1)],'FaceColor','w','BinWidth',0.05);
hold on;
histogram([M2_delay_table.CP_auc_contra(logical(M2_delay_table.flagd2.*findSig(M2_delay_table.CP_p_contra,psig).*(M2_delay_table.CP_auc_contra>0.5)));M2_delay_table.CP_auc_ipsi(logical(M2_delay_table.flagd1.*findSig(M2_delay_table.CP_p_ipsi,psig).*(M2_delay_table.CP_auc_ipsi>0.5)))],'FaceColor','r','BinWidth',0.05);
histogram([M2_delay_table.CP_auc_contra(logical(M2_delay_table.flagd2.*findSig(M2_delay_table.CP_p_contra,psig).*(M2_delay_table.CP_auc_contra<0.5)));M2_delay_table.CP_auc_ipsi(logical(M2_delay_table.flagd1.*findSig(M2_delay_table.CP_p_ipsi,psig).*(M2_delay_table.CP_auc_ipsi<0.5)))],'FaceColor','b','BinWidth',0.05);
ylabel('Cell number');
set(gca,'Xlim',[0,1]);
title('delay');
xlabel('AUC correct v.s.error for prefered stimulus');
set(gca,'FontSize',12);
box off;

%plot distribution of delay AUC M2 and SC together
figCmpM2SCdelay=figure;
set(figCmpM2SCdelay,'Position',[100,100,400,200]);
subplot(1,2,1);
histogram([T_join.delay(T_join.flagd2,2);T_join.delay(T_join.flagd1,1)],'FaceColor','r','BinWidth',0.05);
hold on;
histogram([M2_delay_table.CP_auc_contra(M2_delay_table.flagd2);M2_delay_table.CP_auc_ipsi(M2_delay_table.flagd1)],'FaceColor','b','BinWidth',0.05);
h=legend('SC','M2');
set(h,'box','off');
ylabel('Cell number');
set(gca,'FontSize',12);
set(gca,'Xlim',[0,1]);
box off;
subplot(1,2,2);
histogram([T_join.delay(T_join.flagd2,2);T_join.delay(T_join.flagd1,1)],'EdgeColor','r','BinWidth',0.05,'Normalization','probability','DisplayStyle','stairs');
hold on;
histogram([M2_delay_table.CP_auc_contra(M2_delay_table.flagd2);M2_delay_table.CP_auc_ipsi(M2_delay_table.flagd1)],'EdgeColor','b','BinWidth',0.05,'Normalization','probability','DisplayStyle','stairs');
plot([0.5,0.5],[0,0.3],'k--');
h=legend('SC','M2');
set(h,'box','off');
ylabel('Probabilty');
set(gca,'FontSize',12);
set(gca,'Xlim',[0,1]);
box off;
%% functions
function [flag1,flag2]=fSigFlag(SAUC,pSAUC,CAUC,pCAUC,p_range)
%p_range- the range of significance p value
if length(SAUC)~=length(CAUC) 
    warning('unequal rows of input table');
    return;
end
[indsigS,indsigC,flag1,flag2]=deal(zeros(length(SAUC),1));
for i=1:size(p_range,1) %每一行一个区间段，各行之间取并集
    indsigS=logical((pSAUC<=p_range(i,2)).*(pSAUC>=p_range(i,1))+indsigS);
    indsigC=logical((pCAUC<=p_range(i,2)).*(pCAUC>=p_range(i,1))+indsigC);
end
for irow=1:length(SAUC)
    if (indsigS(irow) && ~indsigC(irow)) || (indsigS(irow) && indsigC(irow) && abs(SAUC(irow)-0.5)>=abs(CAUC(irow)-0.5)) 
        if SAUC(irow)>0.5
            flag2(irow)=1;
        else
            flag1(irow)=1;
        end
    elseif (~indsigS(irow) && indsigC(irow))|| (indsigS(irow) && indsigC(irow) && abs(SAUC(irow)-0.5)<abs(CAUC(irow)-0.5)) 
        if CAUC(irow)>0.5
            flag2(irow)=1;
        else
            flag1(irow)=1;
        end
    end
end
flag1=logical(flag1);
flag2=logical(flag2);
end

function [TsigDataInfo] = fGetSigDataInfoTable(Tin, p_range,strAUCtype)
% summary the significance table using Tin and setted range of significant
% p value
[indsigsound,indsigdelay,indsigresponse,indsiglick]=deal(zeros(size(Tin,1),1));
for i=1:size(p_range,1) %每一行一个区间段，各行之间取并集
    indsigsound=logical((Tin.psound<=p_range(i,2)).*(Tin.psound>=p_range(i,1))+indsigsound);
    indsigdelay=logical((Tin.pdelay<=p_range(i,2)).*(Tin.pdelay>=p_range(i,1))+indsigdelay);
    indsigresponse=logical((Tin.presponse<=p_range(i,2)).*(Tin.presponse>=p_range(i,1))+indsigresponse);
    indsiglick=logical((Tin.plick<=p_range(i,2)).*(Tin.plick>=p_range(i,1))+indsiglick);
end
indsigContraSound=logical(indsigsound.*(Tin.sound>0.5));
indsigIpsiSound=logical(indsigsound.*(Tin.sound<0.5));
indNSSound= ~indsigsound;
indsigContraDelay=logical(indsigdelay.*(Tin.delay>0.5));
indsigIpsiDelay=logical(indsigdelay.*(Tin.delay<0.5));
indNSDelay= ~indsigdelay;
indsigContraResponse=logical(indsigresponse.*(Tin.response>0.5));
indsigIpsiResponse=logical(indsigresponse.*(Tin.response<0.5));
indNSResponse= ~indsigresponse;
indsigContraLick=logical(indsiglick.*(Tin.lick>0.5));
indsigIpsiLick=logical(indsiglick.*(Tin.lick<0.5));
indNSLick= ~indsiglick;
TsigDataInfo=Tin(:,{'animal','date','nROI'});
if strcmp(strAUCtype,'sensory')
    TsigDataInfo.flagSAUCsound(indsigContraSound)=1;
    TsigDataInfo.flagSAUCsound(indsigIpsiSound)=-1;
    TsigDataInfo.flagSAUCsound(indNSSound)=0;
    TsigDataInfo.flagSAUCdelay(indsigContraDelay)=1;
    TsigDataInfo.flagSAUCdelay(indsigIpsiDelay)=-1;
    TsigDataInfo.flagSAUCdelay(indNSDelay)=0;
    TsigDataInfo.flagSAUCresponse(indsigContraResponse)=1;
    TsigDataInfo.flagSAUCresponse(indsigIpsiResponse)=-1;
    TsigDataInfo.flagSAUCresponse(indNSResponse)=0;
    TsigDataInfo.flagSAUClick(indsigContraLick)=1;
    TsigDataInfo.flagSAUClick(indsigIpsiLick)=-1;
    TsigDataInfo.flagSAUClick(indNSLick)=0;
elseif strcmp(strAUCtype,'choice')
    TsigDataInfo.flagCAUCsound(indsigContraSound)=1;
    TsigDataInfo.flagCAUCsound(indsigIpsiSound)=-1;
    TsigDataInfo.flagCAUCsound(indNSSound)=0;
    TsigDataInfo.flagCAUCdelay(indsigContraDelay)=1;
    TsigDataInfo.flagCAUCdelay(indsigIpsiDelay)=-1;
    TsigDataInfo.flagCAUCdelay(indNSDelay)=0;
    TsigDataInfo.flagCAUCresponse(indsigContraResponse)=1;
    TsigDataInfo.flagCAUCresponse(indsigIpsiResponse)=-1;
    TsigDataInfo.flagCAUCresponse(indNSResponse)=0;
    TsigDataInfo.flagCAUClick(indsigContraLick)=1;
    TsigDataInfo.flagCAUClick(indsigIpsiLick)=-1;
    TsigDataInfo.flagCAUClick(indNSLick)=0;
end
end

function [ind]=findSig(var,pSig)
ind=zeros(length(var),1);
for i=1:size(pSig,1) %每一行一个区间段，各行之间取并集
    ind=logical((var<=pSig(i,2)).*(var>=pSig(i,1))+ind);
end
end