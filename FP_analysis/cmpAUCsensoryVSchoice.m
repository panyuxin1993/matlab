%combine table of sensory and choice AUC together and save combined table
%{
%for FP data
clear;
filepath={'F:\FP\summary\SC vgat-mixed-Soma-grouped bysensory-combineCorErr-binCenteredSize3-alignTo-delay onset-timeWindow1-1.5-EpochWindow stim-0.5delay-0.5-shuffle1000-n28-sites8-AUC.mat',
    'F:\FP\summary\SC vgat-mixed-Soma-grouped bychoice-combineCorErr-binCenteredSize3-alignTo-delay onset-timeWindow1-1.5-EpochWindow stim-0.5delay-0.5-shuffle1000-n28-sites8-AUC.mat',
    'F:\FP\summary\SC vglut2-mixed-Soma-grouped bysensory-combineCorErr-binCenteredSize3-alignTo-delay onset-timeWindow1-1.5-EpochWindow stim-0.5delay-0.5-shuffle1000-n14-sites6-AUC.mat',
    'F:\FP\summary\SC vglut2-mixed-Soma-grouped bychoice-combineCorErr-binCenteredSize3-alignTo-delay onset-timeWindow1-1.5-EpochWindow stim-0.5delay-0.5-shuffle1000-n14-sites6-AUC.mat'
    };
celltypestr={'vgat','vgat','vglut2','vglut2'};
AUCtypestr={'sensory','choice','sensory','choice'};
for i=1:4
load(filepath{i});
TAUCsessionTemp=TAUCepochSession(:,1:5);
TAUCsessionTemp=fExtend(TAUCsessionTemp,{celltypestr{i},AUCtypestr{i}},{'celltype','AUCtype'});
TAUCsiteTemp=TAUCepochSite(:,1:5);
TAUCsiteTemp=fExtend(TAUCsiteTemp,{celltypestr{i},AUCtypestr{i}},{'celltype','AUCtype'});
%combine AUC table
if exist('TAUCsession','var')
    TAUCsession=vertcat(TAUCsession,TAUCsessionTemp);
    TAUCsite=vertcat(TAUCsite,TAUCsiteTemp);
else
    TAUCsession=TAUCsessionTemp;
    TAUCsite=TAUCsiteTemp;
end
end

save(['F:\FP\summary\AUCsensoryVSchoiceTable.mat'],'TAUCsession','TAUCsite');
load('F:\FP\summary\AUCsensoryVSchoiceTable.mat');
%plot early AUC comparison
figAUCearly=figure;
pSig=0.05;
set(gcf,'PaperPosition',[0,0,4,4]);
subplot(2,2,1);
fCmpSensoryChoiceAUC(TAUCsession,'early',pSig,'vglut2');
title('early vglut2');
subplot(2,2,2);
fCmpSensoryChoiceAUC(TAUCsession,'late',pSig,'vglut2');
title('late vglut2');
subplot(2,2,3);
fCmpSensoryChoiceAUC(TAUCsession,'early',pSig,'vgat');
title('early vgat');
subplot(2,2,4);
fCmpSensoryChoiceAUC(TAUCsession,'late',pSig,'vgat');
title('late vgat');
saveas(figAUCearly,'F:\FP\summary\AUCsensoryVSchoice_rawFigure.pdf','pdf');
%}

%% load 2P data
%SC 2P data,need to check the assist function(especially for the pSig
%usage
%{
clear;
filepath={'H:\2P\summary\trialTypecor and err-choiceTepochAUC.mat',
    'H:\2P\summary\trialTypecor and err-sensoryTepochAUC.mat'
    };
AUCtypestr={'choice','sensory'};
for i=1:length(filepath)
    load(filepath{i});
    TAUCsessionTemp=TAUC_combine;
    TAUCsessionTemp=fExtend(TAUCsessionTemp,{AUCtypestr{i}},{'AUCtype'});
    %combine AUC table
    if exist('TAUCsession','var')
        TAUCsession=vertcat(TAUCsession,TAUCsessionTemp);
    else
        TAUCsession=TAUCsessionTemp;
    end
end
TAUCsession.celltype=categorical(TAUCsession.celltype);
save(['H:\2P\summary\AUCsensoryVSchoiceTable.mat'],'TAUCsession');
load('H:\2P\summary\AUCsensoryVSchoiceTable.mat');
%plot early AUC comparison
figAUCearly=figure;
set(gcf,'position',[100,100,1000,200]);
pSig=0.05;
epochstr={'ITI','sound','delay','response','lick'};
set(gcf,'PaperPosition',[0,0,2*length(epochstr),2]);
for i=1:length(epochstr)
    subplot(1,length(epochstr),i);
    fCmpSensoryChoiceAUC(TAUCsession,epochstr{i},pSig,'syn');
    title(epochstr{i});
end

saveas(figAUCearly,'H:\2P\summary\SC AUCsensoryVSchoice_rawFigure.pdf','pdf');
%}

%% load 2P data
%M2 2P data, need to check the assist function(especially for the pSig
%usage
%{
clear;
filepath={'D:\download\M2_delay_table.mat'};
AUCtypestr={'sensory','choice'};
load(filepath{1});
for i=1:length(AUCtypestr)
    TAUCsessionTemp= table(table2array(M2_delay_table(:,2*i-1)),table2array(M2_delay_table(:,2*i)),'VariableNames',{'delay','pdelay'});
    TAUCsessionTemp=fExtend(TAUCsessionTemp,{AUCtypestr{i}},{'AUCtype'});
    %combine AUC table
    if exist('TAUCsession','var')
        TAUCsession=vertcat(TAUCsession,TAUCsessionTemp);
    else
        TAUCsession=TAUCsessionTemp;
    end
end
[celltypcell{1:size(TAUCsession,1)}]=deal('syn');
TAUCsession.celltype=categorical(celltypcell');
save(['H:\2P\summary\M2_AUCsensoryVSchoiceTable.mat'],'TAUCsession');
load('H:\2P\summary\M2_AUCsensoryVSchoiceTable.mat');
%plot early AUC comparison
figAUCearly=figure;
set(gcf,'position',[100,100,1000,200]);
pSig=0.05;
epochstr={'ITI','sound','delay','response','lick'};
set(gcf,'PaperPosition',[0,0,2*length(epochstr),2]);
for i=3%length(epochstr)
    subplot(1,length(epochstr),i);
    fCmpSensoryChoiceAUC(TAUCsession,epochstr{i},pSig,'syn');
    title(epochstr{i});
end

saveas(figAUCearly,'H:\2P\summary\M2 AUCsensoryVSchoice_rawFigure.pdf','pdf');
%}
%% assist function
function [Tout]=fCmpSensoryChoiceAUC(Tin,epochstr,pSig,celltypestr)
Tsensory=Tin(logical((Tin.celltype==celltypestr).*(Tin.AUCtype=='sensory')),:);
Tchoice=Tin(logical((Tin.celltype==celltypestr).*(Tin.AUCtype=='choice')),:);
switch epochstr
    case 'early'   
        pAUCsensory=Tsensory.pAUCearly;
        AUCsensory=Tsensory.AUCearly;
        pAUCchoice=Tchoice.pAUCearly;
        AUCchoice=Tchoice.AUCearly;
    case 'late'
        pAUCsensory=Tsensory.pAUClate;
        AUCsensory=Tsensory.AUClate;
        pAUCchoice=Tchoice.pAUClate;
        AUCchoice=Tchoice.AUClate;
    case 'ITI'
        pAUCsensory=Tsensory.pITI;
        AUCsensory=Tsensory.ITI;
        pAUCchoice=Tchoice.pITI;
        AUCchoice=Tchoice.ITI;
    case 'delay'
        pAUCsensory=Tsensory.pdelay;
        AUCsensory=Tsensory.delay;
        pAUCchoice=Tchoice.pdelay;
        AUCchoice=Tchoice.delay;    
    case 'sound'
        pAUCsensory=Tsensory.psound;
        AUCsensory=Tsensory.sound;
        pAUCchoice=Tchoice.psound;
        AUCchoice=Tchoice.sound;   
    case 'response'
        pAUCsensory=Tsensory.presponse;
        AUCsensory=Tsensory.response;
        pAUCchoice=Tchoice.presponse;
        AUCchoice=Tchoice.response;      
    case 'lick'
        pAUCsensory=Tsensory.plick;
        AUCsensory=Tsensory.lick;
        pAUCchoice=Tchoice.plick;
        AUCchoice=Tchoice.lick; 
end
indSigSensory=logical((pAUCsensory<pSig/2)+(pAUCsensory>1-pSig/2));
indSigChoice=logical((pAUCchoice<pSig/2)+(pAUCchoice>1-pSig/2));
% indSigSensory=logical(pAUCsensory<pSig);
% indSigChoice=logical(pAUCchoice<pSig);
indNS=logical(1-logical(indSigSensory+indSigChoice));
indSensory=logical((abs(AUCsensory-0.5)>abs(AUCchoice-0.5)).*(indSigSensory));
indChoice=logical((abs(AUCsensory-0.5)<abs(AUCchoice-0.5)).*(indSigChoice));
colorDot={'FF34E6','2EAF4A','000000'};
colorDot=fHex2RGB(colorDot);
curve1=scatter(AUCsensory(indNS),AUCchoice(indNS),15,colorDot{3});hold on;
curve2=scatter(AUCsensory(indSensory),AUCchoice(indSensory),15,colorDot{1});
curve3=scatter(AUCsensory(indChoice),AUCchoice(indChoice),15,colorDot{2});
set(gca,'Xlim',[0,1],'Ylim',[0,1]);
plot([0,1],[0,1],'k.-.');
plot([0,1],[1,0],'k.-.');
ylabel('choice AUC');
xlabel('sensory AUC');
nSensoryAUC=sum(indSensory);
nChoiceAUC=sum(indChoice);
nNS=sum(indNS);
if length(AUCsensory)>(nNS+nSensoryAUC+nChoiceAUC)
    warning('check larger AUC but not significant');
end
set(gca,'FontSize',12);
text(0.9,0.45,{['sensory',num2str(nSensoryAUC)];['choice',num2str(nChoiceAUC)];['n.s.',num2str(length(AUCsensory)-nChoiceAUC-nSensoryAUC)]},'Unit','Normalized');
%calculate p-value of distribution of seneory v.s. choice selectivity
% %v1 compare all data points
% unsignedAUCsensory=abs(AUCsensory-0.5)+0.5;
% unsignedAUCchoice=abs(AUCchoice-0.5)+0.5;
% if vartest2(unsignedAUCsensory,unsignedAUCchoice)
%     p=ranksum(unsignedAUCsensory,unsignedAUCchoice);
%     pstr={'unequal variance,' ,['ranksum test p=',num2str(p)]};
% else
%     [h,p]=ttest(unsignedAUCsensory,unsignedAUCchoice);
%     pstr={'equal variance, ',['p=',num2str(p)]};
% end
% text(0.9,0.85,pstr,'Unit','Normalized');
%v1.1 compare all data points using ranksum test
%{
unsignedAUCsensory=abs(AUCsensory-0.5)+0.5;
unsignedAUCchoice=abs(AUCchoice-0.5)+0.5;
p=signrank(unsignedAUCsensory,unsignedAUCchoice);
pstr=['Wilcoxon signed rank test p=',num2str(p)];
text(0.9,0.85,pstr,'Unit','Normalized');
%}

% %v2 compare only data points with significant choice/senory selectivity
% unsignedAUCsensory=abs(AUCsensory(~indNS)-0.5)+0.5;
% unsignedAUCchoice=abs(AUCchoice(~indNS)-0.5)+0.5;
% if vartest2(unsignedAUCsensory,unsignedAUCchoice)
%     p=ranksum(unsignedAUCsensory,unsignedAUCchoice);
%     pstr=['unequal variance, p=',num2str(p)];
% else
%     [h,p]=ttest(unsignedAUCsensory,unsignedAUCchoice);
%     pstr=['equal variance, p=',num2str(p)];
% end
% text(0.1,0.4,pstr,'Unit','Normalized');
% h=legend([curve2,curve3,curve1],'stronger sensory AUC','stronger choice AUC','not significant AUC');
% set(h,'box','off');
end
%extend table with variables that replicate a value as categorical
function [Tout]=fExtend(Tin,varcell,namecell)
%varcell and namecell should have same size
for i=1:length(varcell)%for each var as a col
    [vartemp{1:size(Tin,1)}]=deal(varcell{i});
    vartemp=categorical(vartemp);
    vartemp=reshape(vartemp,[],1);
    Tin=addvars(Tin,vartemp,'NewVariableNames',namecell{i});
    clear vartemp;
end
Tout=Tin;
end