%plot correlation between activities and reaction time(RT) for preferred
%side of individual neurons (that with significant delay selectivity, AUC) 
%using only correct trials.

% close all;
clear;
[num,txt,raw] =xlsread('D:\xulab\project\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
savepath='H:\2P\summary';
trialTypeStr='cor';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
selectivityType='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
correlationMethod='cells';
displayedCellSelectivity='choice';%{'choice','sensory','sensory and choice'}
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
    T_currentSession=fGetCorrelationRT(T.file_path{indrow},T.animal{indrow},T.date{indrow},T.field{indrow},T.cell_type{indrow},trialTypeStr,selectivityType,correlationMethod);
    if exist('TcorrelationRT','var')
        TcorrelationRT=vertcat(TcorrelationRT,T_currentSession);
    else
        TcorrelationRT=T_currentSession;
    end
end
save([savepath,filesep,'trialType',trialTypeStr,'-selectivity',selectivityType,'RTcorrelation.mat'],'TcorrelationRT');
%}
load([savepath,filesep,'trialType',trialTypeStr,'-selectivity',selectivityType,'RTcorrelation.mat']);


%% choose cells with significant delay/sound choice AUC
psig=[0,0.05];
pSigAUC=[0,0.025;0.975,1];
trialTypeStr_datapool='cor and err';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
AUCtype_datapool='sensory';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
load([savepath,filesep,'trialType',trialTypeStr_datapool,'-',AUCtype_datapool,'TepochAUC.mat']);
celltype='syn';
T_AUCs=TAUC_combine(strcmp(TAUC_combine.celltype,celltype),:);
TsigDataInfosensory = fGetSigDataInfoTable(T_AUCs, pSigAUC,AUCtype_datapool);

AUCtype_datapool='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
load([savepath,filesep,'trialType',trialTypeStr_datapool,'-',AUCtype_datapool,'TepochAUC.mat']);
T_AUCc=TAUC_combine(strcmp(TAUC_combine.celltype,celltype),:);
TsigDataInfochoice = fGetSigDataInfoTable(T_AUCc, pSigAUC,AUCtype_datapool);

T_sig=join(TsigDataInfosensory,TsigDataInfochoice,'Keys',{'animal','date','nROI'});%choose neurons with either choice or sensory selectivity
% T_sig=TsigDataInfochoice;%choose neurons with only choice selectivity


TcorrRT=TcorrelationRT(strcmp(TcorrelationRT.celltype,celltype),:);
T_join=join(TcorrRT,T_sig,'Keys',{'animal','date','nROI'});

[flag1,flag2]=fSigFlag(T_AUCs.ITI,T_AUCs.pITI,T_AUCc.ITI,T_AUCc.pITI,pSigAUC);
T_join=addvars(T_join,flag1,flag2,'NewVariableNames',{'flagi1','flagi2'});
[flag1,flag2]=fSigFlag(T_AUCs.sound,T_AUCs.psound,T_AUCc.sound,T_AUCc.psound,pSigAUC);
T_join=addvars(T_join,flag1,flag2,'NewVariableNames',{'flags1','flags2'});
[flag1,flag2]=fSigFlag(T_AUCs.delay,T_AUCs.pdelay,T_AUCc.delay,T_AUCc.pdelay,pSigAUC);
T_join=addvars(T_join,flag1,flag2,'NewVariableNames',{'flagd1','flagd2'});
[flag1,flag2]=fSigFlag(T_AUCs.response,T_AUCs.presponse,T_AUCc.response,T_AUCc.presponse,pSigAUC);
T_join=addvars(T_join,flag1,flag2,'NewVariableNames',{'flagr1','flagr2'});
[flag1,flag2]=fSigFlag(T_AUCs.lick,T_AUCs.plick,T_AUCc.lick,T_AUCc.plick,pSigAUC);
T_join=addvars(T_join,flag1,flag2,'NewVariableNames',{'flagl1','flagl2'});
%% plot histogram
figRT=figure;
if strcmp(displayedCellSelectivity,'sensory and choice')
    set(figRT,'Position',[100,100,600,200]);
    subplot(1,3,1);
    hold on;
    histogram([T_join.ITI(T_join.flagi2,2);T_join.ITI(T_join.flagi1,1)],'FaceColor','w','BinWidth',0.05);
    histogram([T_join.ITI(logical(T_join.flagi2.*findSig(T_join.pITI(:,2),psig).*(T_join.ITI(:,2)>0)),2);T_join.ITI(logical(T_join.flagi1.*findSig(T_join.pITI(:,1),psig).*(T_join.ITI(:,1)>0)),1)],'FaceColor','r','BinWidth',0.05);
    histogram([T_join.ITI(logical(T_join.flagi2.*findSig(T_join.pITI(:,2),psig).*(T_join.ITI(:,2)<0)),2);T_join.ITI(logical(T_join.flagi1.*findSig(T_join.pITI(:,1),psig).*(T_join.ITI(:,1)<0)),1)],'FaceColor','b','BinWidth',0.05);
    y_lim=get(gca,'Ylim');
    plot([0,0],[0,y_lim(end)],'k-');
    x_mean=nanmean([T_join.ITI(T_join.flagi2,2);T_join.ITI(T_join.flagi1,1)]);
    plot(x_mean,y_lim(end),'k','Marker','v');
    ylabel('Cell number');
    title('ITI');
    text(0.6,0.8,['mean ',num2str(round(x_mean,3))],'Units','normalized');
    set(gca,'Xlim',[-0.5,0.5]);
    set(gca,'FontSize',12);
    %%histogram(T_join.sound,'BinWidth',0.05);
    subplot(1,3,2);
    hold on;
    histogram([T_join.sound(T_join.flags2,2);T_join.sound(T_join.flags1,1)],'FaceColor','w','BinWidth',0.05);
    histogram([T_join.sound(logical(T_join.flags2.*findSig(T_join.psound(:,2),psig).*(T_join.sound(:,2)>0)),2);T_join.sound(logical(T_join.flags1.*findSig(T_join.psound(:,1),psig).*(T_join.sound(:,1)>0)),1)],'FaceColor','r','BinWidth',0.05);
    histogram([T_join.sound(logical(T_join.flags2.*findSig(T_join.psound(:,2),psig).*(T_join.sound(:,2)<0)),2);T_join.sound(logical(T_join.flags1.*findSig(T_join.psound(:,1),psig).*(T_join.sound(:,1)<0)),1)],'FaceColor','b','BinWidth',0.05);
    y_lim=get(gca,'Ylim');
    plot([0,0],[0,y_lim(end)],'k-');
    x_mean=nanmean([T_join.sound(T_join.flags2,2);T_join.sound(T_join.flags1,1)]);
    plot(x_mean,y_lim(end),'k','Marker','v');
    
    title('sound');
    text(0.6,0.8,['mean ',num2str(round(x_mean,3))],'Units','normalized');
    set(gca,'Xlim',[-0.5,0.5]);
    set(gca,'FontSize',12);
    
    subplot(1,3,3);
    % histogram(T_join.sound,'BinWidth',0.05);
    hold on;
    histogram([T_join.delay(T_join.flagd2,2);T_join.delay(T_join.flagd1,1)],'FaceColor','w','BinWidth',0.05);
    histogram([T_join.delay(logical(T_join.flagd2.*findSig(T_join.pdelay(:,2),psig).*(T_join.delay(:,2)>0)),2);T_join.delay(logical(T_join.flagd1.*findSig(T_join.pdelay(:,1),psig).*(T_join.delay(:,1)>0)),1)],'FaceColor','r','BinWidth',0.05);
    histogram([T_join.delay(logical(T_join.flagd2.*findSig(T_join.pdelay(:,2),psig).*(T_join.delay(:,2)<0)),2);T_join.delay(logical(T_join.flagd1.*findSig(T_join.pdelay(:,1),psig).*(T_join.delay(:,1)<0)),1)],'FaceColor','b','BinWidth',0.05);
    y_lim=get(gca,'Ylim');
    plot([0,0],[0,y_lim(end)],'k-');
    x_mean=nanmean([T_join.delay(T_join.flagd2,2);T_join.delay(T_join.flagd1,1)]);
    plot(x_mean,y_lim(end),'k','Marker','v');
    title('delay');
    text(0.6,0.8,['mean ',num2str(round(x_mean,3))],'Units','normalized');
    xlabel('Correlation coefficients between RT and acitivities');
    set(gca,'Xlim',[-0.5,0.5]);
    set(gca,'FontSize',12);
    h=legend('n.s.','+','-');
    set(h,'box','off');
elseif strcmp(displayedCellSelectivity,'choice')
    set(figRT,'Position',[100,100,400,600]);
    subplot(3,2,1);
    hold on;
    histogram([T_join.sound(T_join.flagCAUCsound>0,2);T_join.sound(T_join.flagCAUCsound<0,1)],'FaceColor','w','BinWidth',0.05);
    histogram([T_join.sound(logical((T_join.flagCAUCsound>0).*findSig(T_join.psound(:,2),psig).*(T_join.sound(:,2)>0)),2);T_join.sound(logical((T_join.flagCAUCsound<0).*findSig(T_join.psound(:,1),psig).*(T_join.sound(:,1)>0)),1)],'FaceColor','r','BinWidth',0.05);
    histogram([T_join.sound(logical((T_join.flagCAUCsound>0).*findSig(T_join.psound(:,2),psig).*(T_join.sound(:,2)<0)),2);T_join.sound(logical((T_join.flagCAUCsound<0).*findSig(T_join.psound(:,1),psig).*(T_join.sound(:,1)<0)),1)],'FaceColor','b','BinWidth',0.05);
    y_lim=get(gca,'Ylim');
    plot([0,0],[0,y_lim(end)],'k-');
    x_mean=nanmean([T_join.sound(T_join.flags2,2);T_join.sound(T_join.flags1,1)]);
    plot(x_mean,y_lim(end),'k','Marker','v');
    ylabel('Cell number');
    title('sound');
    text(0.6,0.8,['mean ',num2str(round(x_mean,3))],'Units','normalized');
    set(gca,'Xlim',[-0.5,0.5]);
    set(gca,'FontSize',12);
    
    subplot(3,2,2);
    % histogram(T_join.sound,'BinWidth',0.05);
    hold on;
    histogram([T_join.delay(T_join.flagCAUCdelay>0,2);T_join.delay(T_join.flagCAUCdelay<0,1)],'FaceColor','w','BinWidth',0.05);
    histogram([T_join.delay(logical((T_join.flagCAUCdelay>0).*findSig(T_join.pdelay(:,2),psig).*(T_join.delay(:,2)>0)),2);T_join.delay(logical((T_join.flagCAUCdelay<0).*findSig(T_join.pdelay(:,1),psig).*(T_join.delay(:,1)>0)),1)],'FaceColor','r','BinWidth',0.05);
    histogram([T_join.delay(logical((T_join.flagCAUCdelay>0).*findSig(T_join.pdelay(:,2),psig).*(T_join.delay(:,2)<0)),2);T_join.delay(logical((T_join.flagCAUCdelay<0).*findSig(T_join.pdelay(:,1),psig).*(T_join.delay(:,1)<0)),1)],'FaceColor','b','BinWidth',0.05);
    y_lim=get(gca,'Ylim');
    plot([0,0],[0,y_lim(end)],'k-');
    x_mean=nanmean([T_join.delay(T_join.flagd2,2);T_join.delay(T_join.flagd1,1)]);
    plot(x_mean,y_lim(end),'k','Marker','v');
    title('delay');
    text(0.6,0.8,['mean ',num2str(round(x_mean,3))],'Units','normalized');
    xlabel('Correlation coefficients between RT and acitivities');
    set(gca,'Xlim',[-0.5,0.5]);
    set(gca,'FontSize',12);
    h=legend('n.s.','+','-');
    set(h,'box','off');
    
    subplot(3,2,3);%plot contra
    hold on;
    histogram([T_join.sound(T_join.flagCAUCsound~=0,2)],'FaceColor','w','BinWidth',0.05);
    histogram([T_join.sound(logical((T_join.flagCAUCsound~=0).*findSig(T_join.psound(:,2),psig).*(T_join.sound(:,2)>0)),2)],'FaceColor','r','BinWidth',0.05);
    histogram([T_join.sound(logical((T_join.flagCAUCsound~=0).*findSig(T_join.psound(:,2),psig).*(T_join.sound(:,2)<0)),2)],'FaceColor','b','BinWidth',0.05);
    y_lim=get(gca,'Ylim');
    plot([0,0],[0,y_lim(end)],'k-');
    x_mean=nanmean(T_join.sound(T_join.flagCAUCsound~=0,2));
    plot(x_mean,y_lim(end),'k','Marker','v');
    ylabel('Cell number');
    title('sound');
    text(0.6,0.8,['mean ',num2str(round(x_mean,3))],'Units','normalized');
    set(gca,'Xlim',[-0.5,0.5]);
    set(gca,'FontSize',12);
    
    subplot(3,2,4);%plot contra
    % histogram(T_join.sound,'BinWidth',0.05);
    hold on;
    histogram([T_join.delay(T_join.flagCAUCdelay~=0,2)],'FaceColor','w','BinWidth',0.05);
    histogram([T_join.delay(logical((T_join.flagCAUCdelay~=0).*findSig(T_join.pdelay(:,2),psig).*(T_join.delay(:,2)>0)),2)],'FaceColor','r','BinWidth',0.05);
    histogram([T_join.delay(logical((T_join.flagCAUCdelay~=0).*findSig(T_join.pdelay(:,2),psig).*(T_join.delay(:,2)<0)),2)],'FaceColor','b','BinWidth',0.05);    
    y_lim=get(gca,'Ylim');
    plot([0,0],[0,y_lim(end)],'k-');
    x_mean=nanmean(T_join.delay(T_join.flagCAUCdelay~=0,2));
    plot(x_mean,y_lim(end),'k','Marker','v');
    title('delay');
    text(0.6,0.8,['mean ',num2str(round(x_mean,3))],'Units','normalized');
    xlabel('Correlation coefficients between RT and acitivities');
    set(gca,'Xlim',[-0.5,0.5]);
    set(gca,'FontSize',12);
    
    subplot(3,2,5);%plot ipsi
    hold on;
    histogram([T_join.sound(T_join.flagCAUCsound~=0,1)],'FaceColor','w','BinWidth',0.05);
    histogram([T_join.sound(logical((T_join.flagCAUCsound~=0).*findSig(T_join.psound(:,1),psig).*(T_join.sound(:,1)>0)),1)],'FaceColor','r','BinWidth',0.05);
    histogram([T_join.sound(logical((T_join.flagCAUCsound~=0).*findSig(T_join.psound(:,1),psig).*(T_join.sound(:,1)<0)),1)],'FaceColor','b','BinWidth',0.05);
    y_lim=get(gca,'Ylim');
    plot([0,0],[0,y_lim(end)],'k-');
    x_mean=nanmean(T_join.sound(T_join.flagCAUCsound~=0,1));
    plot(x_mean,y_lim(end),'k','Marker','v');
    ylabel('Cell number');
    title('sound');
    text(0.6,0.8,['mean ',num2str(round(x_mean,3))],'Units','normalized');
    set(gca,'Xlim',[-0.5,0.5]);
    set(gca,'FontSize',12);
    
    subplot(3,2,6);%plot ipsi
    % histogram(T_join.sound,'BinWidth',0.05);
    hold on;
    histogram([T_join.delay(T_join.flagCAUCdelay~=0,1)],'FaceColor','w','BinWidth',0.05);
    histogram([T_join.delay(logical((T_join.flagCAUCdelay~=0).*findSig(T_join.pdelay(:,1),psig).*(T_join.delay(:,1)>0)),1)],'FaceColor','r','BinWidth',0.05);
    histogram([T_join.delay(logical((T_join.flagCAUCdelay~=0).*findSig(T_join.pdelay(:,1),psig).*(T_join.delay(:,1)<0)),1)],'FaceColor','b','BinWidth',0.05);    
    y_lim=get(gca,'Ylim');
    plot([0,0],[0,y_lim(end)],'k-');
    x_mean=nanmean(T_join.delay(T_join.flagCAUCdelay~=0,1));
    plot(x_mean,y_lim(end),'k','Marker','v');
    title('delay');
    text(0.6,0.8,['mean ',num2str(round(x_mean,3))],'Units','normalized');
    xlabel('Correlation coefficients between RT and acitivities');
    set(gca,'Xlim',[-0.5,0.5]);
    set(gca,'FontSize',12);
end
set(figRT,'PaperPosition',[1,1,4,2]);
saveas(figRT,[savepath,filesep,'correlation_RT_vs_activities.pdf'],'pdf');
%% function
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
    if (indsigS(irow) && ~indsigC(irow)) || (indsigS(irow) && indsigC(irow) && (abs(SAUC(irow)-0.5)>=abs(CAUC(irow)-0.5))) 
        if SAUC(irow)>0.5
            flag2(irow)=1;
        else
            flag1(irow)=1;
        end
    elseif (~indsigS(irow) && indsigC(irow))|| (indsigS(irow) && indsigC(irow) && (abs(SAUC(irow)-0.5)<abs(CAUC(irow)-0.5))) 
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