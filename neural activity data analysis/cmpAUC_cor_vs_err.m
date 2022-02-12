close all;
[num,txt,raw] =xlsread('C:\Users\PYX\Documents\DataSummary\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
savepath='E:\2P\summary\AUC_cor_vs_err';
clearvars TAUC_combine Tmean_combine;
trialTypeStrPool={'cor','err'};%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
celltypePool={'M2','syn','vglut2','vgat'};
AUCtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
AUCCorrectedMethod='None';%'SensoryChoiceOrthogonalSubtraction';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
manipulation='control';
%% calculate AUC
for i_celltype=1:length(celltypePool)
    celltype=celltypePool{i_celltype};
    %{
    for i_trialtype=1:length(trialTypeStrPool)
        trialTypeStr=trialTypeStrPool{i_trialtype};
        T=cell2table(raw(2:end,1:15));
        T.Properties.VariableNames=strrep(raw(1,1:15),' ','_');%table variable name can't have ' ',so replace them
        ind_session1=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*(~contains(T.behavior_performance,'probe'));%not probe session
        ind_session2=logical(strcmp(T.cell_type,celltype).*contains(T.ROI_type,'soma'));
        % ind_session2=logical(~strcmp(T.cell_type,'M2'));
        ind_session= ind_session1 & ind_session2;
        ind_session=find(ind_session);
        n_session=length(ind_session);
        animal_unique=unique(T.animal(ind_session));
        n_animal=length(animal_unique);
        
        for i_session=1:n_session
            indrow=ind_session(i_session);
            [TAUC_currentSession, Tmean_currentSession] = fGetEpochAUCtableASession(T.file_path{indrow},T.animal{indrow},T.date{indrow},T.field{indrow},T.cell_type{indrow},T.ROI_type{indrow},trialTypeStr,AUCtype,AUCCorrectedMethod);
            %     [TAUC_currentSession, Tmean_currentSession] = fGetEpochAUCtableASession_NPseg(T.file_path{indrow},T.animal{indrow},T.date{indrow},T.field{indrow},T.cell_type{indrow},trialTypeStr,AUCtype,AUCCorrectedMethod);
            if isempty(TAUC_currentSession)
                continue;
            end
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
        save([savepath,filesep,celltype,'-trialType',trialTypeStr,'-',AUCtype,'TepochAUC.mat'],'TAUC_combine');
        clearvars TAUC_combine Tmean_combine;
    end
    %}
    %plot sactter plot comparing cor and err AUC
    clearvars TAUCsession filepath;
    for i_trialtype=1:length(trialTypeStrPool) %usually is 2, cor and err
        trialTypeStr=trialTypeStrPool{i_trialtype};
        filepath{i_trialtype}=[savepath,filesep,celltype,'-trialType',trialTypeStr,'-',AUCtype,'TepochAUC.mat'];
    end
    trialtypestr={'correct','error'};
    for i=1:length(filepath)
        load(filepath{i});
        TAUCsessionTemp=TAUC_combine;
        TAUCsessionTemp=fExtend(TAUCsessionTemp,{trialtypestr{i}},{'trial_type'});
        %combine AUC table
        if exist('TAUCsession','var')
            TAUCsession=vertcat(TAUCsession,TAUCsessionTemp);
        else
            TAUCsession=TAUCsessionTemp;
        end
    end
    TAUCsession.celltype=categorical(TAUCsession.celltype);
    save([savepath,filesep,celltype,'-AUCcorVSerrTable.mat'],'TAUCsession');
    load([savepath,filesep,celltype,'-AUCcorVSerrTable.mat']);
    celltypestr=celltype;
    %plot early AUC comparison
    figAUCearly=figure;
    set(gcf,'position',[100,100,1200,200]);
    pSig=0.05;
    epochstr={'ITI','sound','delay','response','lick'};
    set(gcf,'PaperPosition',[0,0,2*length(epochstr),2]);
    for i=1:length(epochstr)
        subplot(1,length(epochstr),i);
        fCmpCorErrAUC(TAUCsession,epochstr{i},pSig,celltypestr);
        title(epochstr{i});
    end
    set(gcf,'paperPosition',[0,0,11,2]);
    celltypestr_save=strrep(celltypestr,' ','_');
    saveas(figAUCearly,[savepath,filesep,celltypestr,'_AUC_CorVSErr_rawFigure.pdf'],'pdf');

end


%% assistant function
function [Tout]=fCmpCorErrAUC(Tin,epochstr,pSig,celltypestr)
%correct vs. error
Tcor=Tin(logical((Tin.celltype==celltypestr).*(Tin.trial_type=='correct')),:);
Terr=Tin(logical((Tin.celltype==celltypestr).*(Tin.trial_type=='error')),:);
switch epochstr
    case 'early'
        pAUCcor=Tcor.pAUCearly;
        AUCcor=Tcor.AUCearly;
        pAUCerr=Terr.pAUCearly;
        AUCerr=Terr.AUCearly;
    case 'late'
        pAUCcor=Tcor.pAUClate;
        AUCcor=Tcor.AUClate;
        pAUCerr=Terr.pAUClate;
        AUCerr=Terr.AUClate;
    case 'ITI'
        pAUCcor=Tcor.pITI;
        AUCcor=Tcor.ITI;
        pAUCerr=Terr.pITI;
        AUCerr=Terr.ITI;
    case 'delay'
        pAUCcor=Tcor.pdelay;
        AUCcor=Tcor.delay;
        pAUCerr=Terr.pdelay;
        AUCerr=Terr.delay;
    case 'sound'
        pAUCcor=Tcor.psound;
        AUCcor=Tcor.sound;
        pAUCerr=Terr.psound;
        AUCerr=Terr.sound;
    case 'response'
        pAUCcor=Tcor.presponse;
        AUCcor=Tcor.response;
        pAUCerr=Terr.presponse;
        AUCerr=Terr.response;
    case 'lick'
        pAUCcor=Tcor.plick;
        AUCcor=Tcor.lick;
        pAUCerr=Terr.plick;
        AUCerr=Terr.lick;
end
indSigCor=logical((pAUCcor<pSig/2)+(pAUCcor>1-pSig/2));
indSigErr=logical((pAUCerr<pSig/2)+(pAUCerr>1-pSig/2));
% indSigCor=logical(pAUCcor<pSig);
% indSigErr=logical(pAUCerr<pSig);
indNS=logical(1-logical(indSigCor+indSigErr));
indCor=logical((abs(AUCcor-0.5)>abs(AUCerr-0.5)).*(indSigCor));
indErr=logical((abs(AUCcor-0.5)<abs(AUCerr-0.5)).*(indSigErr));
colorDot={'FF34E6','2EAF4A','000000'};
colorDot=fHex2RGB(colorDot);
colorShade={'E7E7E7'};
colorShade=fHex2RGB(colorShade);
xpatch=[0,0.5,1];hold on;
ypatch1=[0,0.5,0];
ypatch2=[1,0.5,1];
patch(xpatch,ypatch1,colorShade{1},'EdgeColor','none','FaceAlpha',0.5);
patch(xpatch,ypatch2,colorShade{1},'EdgeColor','none','FaceAlpha',0.5);
plot([0,1],[0.5,0.5],'k-');
plot([0.5,0.5],[1,0],'k-');
curve1=scatter(AUCcor(indNS),AUCerr(indNS),15,colorDot{3});hold on;
curve2=scatter(AUCcor(indCor),AUCerr(indCor),15,colorDot{1});
curve3=scatter(AUCcor(indErr),AUCerr(indErr),15,colorDot{2});
set(gca,'Xlim',[0,1],'Ylim',[0,1]);

ylabel('error AUC');
xlabel('correct AUC');
nCorAUC=sum(indCor);
nErrAUC=sum(indErr);
nNS=sum(indNS);
if length(AUCcor)>(nNS+nCorAUC+nErrAUC)
    warning('check larger AUC but not significant');
end
set(gca,'FontSize',12);
text(0.9,0.45,{['correct',num2str(nCorAUC)];['error',num2str(nErrAUC)];['n.s.',num2str(length(AUCcor)-nErrAUC-nCorAUC)]},'Unit','Normalized');
%calculate p-value of distribution of seneory v.s. choice selectivity
% %v1 compare all data points
% unsignedAUCcor=abs(AUCcor-0.5)+0.5;
% unsignedAUCerr=abs(AUCerr-0.5)+0.5;
% if vartest2(unsignedAUCcor,unsignedAUCerr)
%     p=ranksum(unsignedAUCcor,unsignedAUCerr);
%     pstr={'unequal variance,' ,['ranksum test p=',num2str(p)]};
% else
%     [h,p]=ttest(unsignedAUCcor,unsignedAUCerr);
%     pstr={'equal variance, ',['p=',num2str(p)]};
% end
% text(0.9,0.85,pstr,'Unit','Normalized');
%v1.1 compare all data points using ranksum test
%
unsignedAUCcor=abs(AUCcor-0.5)+0.5;
unsignedAUCerr=abs(AUCerr-0.5)+0.5;
p=signrank(unsignedAUCcor,unsignedAUCerr);
pstr=['Wilcoxon signed rank test p=',num2str(p)];
text(0.9,0.85,pstr,'Unit','Normalized');
%}

% %v2 compare only data points with significant choice/senory selectivity
% unsignedAUCcor=abs(AUCcor(~indNS)-0.5)+0.5;
% unsignedAUCerr=abs(AUCerr(~indNS)-0.5)+0.5;
% if vartest2(unsignedAUCcor,unsignedAUCerr)
%     p=ranksum(unsignedAUCcor,unsignedAUCerr);
%     pstr=['unequal variance, p=',num2str(p)];
% else
%     [h,p]=ttest(unsignedAUCcor,unsignedAUCerr);
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