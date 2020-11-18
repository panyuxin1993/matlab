%plot each cells' data from all treatment, each row a cell, source of data
%come from xlsx file.
rootpath='D:\xulab\project\slice recording\data3\';
xlsfile=strcat(rootpath,'raw data summary.xlsx');
[num,txt,raw] = xlsread(xlsfile,1);
celltype = {'vglut2+','vglut2-','vgat+','vgat-'};
ncol=size(num,2);
ncell=zeros(1,length(celltype));
ind_celltype=cell(1,length(celltype));%each cell contain index of individual celltype neurons
for n_celltype=1:length(celltype)
    temp=cellfun(@(x) strcmp(x,celltype{n_celltype}),txt);
    temp2=sum(temp,2);
    ind_celltype{1,n_celltype}=find(temp2)-1;
    ncell(1,n_celltype)=sum(temp2);
end
nrow=6;
numfig=ceil(ncell/nrow);
cd(rootpath);
%{
%each column means one experiment
%parameters for EPSP, col-2
Amp_mean=cell(ncol,length(celltype));
amplitude=cell(ncol,length(celltype));%all amp for each trial
Amp_average=cell(ncol,length(celltype));
pResponsive=cell(ncol,length(celltype));
latencyPeak=cell(ncol,length(celltype));
latencyOnset=cell(ncol,length(celltype));
%col-1, parameters for 10P1, many parameters not need
%col-3~5, parameters for PPR
%col-6, parameters for 10P2,
for n_celltype=1:length(celltype) 
    %parameters for EPSP
    amplitude{2,n_celltype}=cell(2,ncell(1,n_celltype));%1d amp for all cell or for significant cells
    Amp_mean{2,n_celltype}=zeros(ncell(1,n_celltype),2);%1d cells, 2d amp_mean for all cell or for significant cells
    Amp_average{2,n_celltype}=zeros(ncell(1,n_celltype),1);
    pResponsive{2,n_celltype}=zeros(ncell(1,n_celltype),2);%1d cells, 2d pRes for all cell or for significant cells
    latencyPeak{2,n_celltype}=zeros(ncell(1,n_celltype),1);%latency of response onset
    latencyOnset{2,n_celltype}=zeros(ncell(1,n_celltype),1);%latency of response peak
    %parameters for 10P1
    Amp_mean{1,n_celltype}=zeros(ncell(1,n_celltype),10);%
    %parameters for 10P2
    Amp_mean{6,n_celltype}=zeros(ncell(1,n_celltype),10);%
    latencyPeak{6,n_celltype}=zeros(ncell(1,n_celltype),10);%latency of response onset
    latencyOnset{6,n_celltype}=zeros(ncell(1,n_celltype),10);%latency of response peak
    amplitude{6,n_celltype}=cell(1,ncell(1,n_celltype));    
    %parameters for PPR
    for i=3:5
        Amp_mean{i,n_celltype}=zeros(ncell(1,n_celltype),2);%
        latencyPeak{i,n_celltype}=zeros(ncell(1,n_celltype),2);%latency of response onset
        latencyOnset{i,n_celltype}=zeros(ncell(1,n_celltype),2);%latency of response peak
    end
    
    for nfig=1:numfig(n_celltype)
        fig=figure;
        set(gcf,'Position',[0,0,400*ncol,300*nrow]);
        nloop=min(nrow,ncell(1,n_celltype)-(nfig-1)*nrow);
        for n_abf=1:nloop
            ylim=zeros(1,2);
            i_cell=n_abf+(nfig-1)*nrow;
            %10p1
            if ~isnan(num(ind_celltype{1,n_celltype}(i_cell),1))%always have file
                abf10p1 = strcat(txt{ind_celltype{1,n_celltype}(i_cell)+1,10},'\AT\',num2str(num(ind_celltype{1,n_celltype}(i_cell),1)),'.abf');
                subplot(nrow,ncol,6*n_abf-6+1);
                [ Amp_mean{1,n_celltype}(i_cell,:),~,~,~,~,~,~ ] = fplotPlasticityCase( abf10p1 );
                ylim_temp=get(gca,'ylim');
                ylim(1)=ylim_temp(1);
                ylim(2)=ylim_temp(end);
            end
            %EPSP
            if ~isnan(num(ind_celltype{1,n_celltype}(i_cell),2))
                abfEPSP = strcat(txt{ind_celltype{1,n_celltype}(i_cell)+1,10},'\AT\',num2str(num(ind_celltype{1,n_celltype}(i_cell),2)),'.abf');
                subplot(nrow,ncol,6*n_abf-6+2);
                [ latencyPeak{2,n_celltype}(i_cell), latencyOnset{2,n_celltype}(i_cell), amplitude{2,n_celltype}{1,i_cell},Amp_mean{2,n_celltype}(i_cell,1),amplitude{2,n_celltype}{2,i_cell},Amp_mean{2,n_celltype}(i_cell,2),Amp_average{2,n_celltype}(i_cell),pResponsive{2,n_celltype}(i_cell,1),pResponsive{2,n_celltype}(i_cell,2)  ] = fplotEPSPCase( abfEPSP );
                ylim_temp=get(gca,'ylim');
                ylim(1)=min(ylim(1),ylim_temp(1));
                ylim(2)=max(ylim(2),ylim_temp(end));
            else
                latencyPeak{2,n_celltype}(i_cell)=nan;
                latencyOnset{2,n_celltype}(i_cell)=nan;
                Amp_mean{2,n_celltype}(i_cell,1:2)=nan;
                Amp_average{2,n_celltype}(i_cell)=nan;
                pResponsive{2,n_celltype}(i_cell,1:2)=nan;
            end
            %PPR
            for i=3:5 %all PPR experiment
                if ~isnan(num(ind_celltype{1,n_celltype}(i_cell),i))
                    abfPPR = strcat(txt{ind_celltype{1,n_celltype}(i_cell)+1,10},'\AT\',num2str(num(ind_celltype{1,n_celltype}(i_cell),i)),'.abf');             
                    subplot(nrow,ncol,6*n_abf-6+i);
                    [ Amp_mean{i,n_celltype}(i_cell,:), currentPPR,latencyOnset{i,n_celltype}(i_cell,:),latencyPeak{i,n_celltype}(i_cell,:),Amp] = fplotPPRCase( abfPPR )
%                     ylim_temp=get(gca,'ylim');
%                     ylim(1)=min(ylim(1),ylim_temp(1));
%                     ylim(2)=max(ylim(2),ylim_temp(end));
                else
                    latencyPeak{i,n_celltype}(i_cell)=nan;
                    latencyOnset{i,n_celltype}(i_cell)=nan;
                    Amp_mean{i,n_celltype}(i_cell)=nan;
                end
            end
            %10p2
            if ~isnan(num(ind_celltype{1,n_celltype}(i_cell),6))
                abf10p2 = strcat(txt{ind_celltype{1,n_celltype}(i_cell)+1,10},'\AT\',num2str(num(ind_celltype{1,n_celltype}(i_cell),6)),'.abf');
                subplot(nrow,ncol,6*n_abf);
                [  Amp_mean{6,n_celltype}(i_cell,:),slope,PPR,plasticity,latencyOnset{6,n_celltype}(i_cell,:),latencyPeak{6,n_celltype}(i_cell,:),amplitude{6,n_celltype}{i_cell} ] = fplotPlasticityCase( abf10p2 );
%                 ylim_temp=get(gca,'ylim');
%                 ylim(1)=min(ylim(1),ylim_temp(1));
%                 ylim(2)=max(ylim(2),ylim_temp(end));
            else
                latencyPeak{6,n_celltype}(i_cell)=nan;
                latencyOnset{6,n_celltype}(i_cell)=nan;
                Amp_mean{6,n_celltype}(i_cell)=nan;
            end
            for i=1:2
                if ~isnan(num(ind_celltype{1,n_celltype}(i_cell),i))
                    subplot(nrow,ncol,6*n_abf-6+i);
                    set(gca,'ylim',[ylim(1), ylim(2)]);
                end
            end
        end
        saveas(fig,strcat(celltype{1,n_celltype},'-fig-',num2str(nfig),'.fig'),'fig');
        saveas(fig,strcat(celltype{1,n_celltype},'-fig-',num2str(nfig),'.pdf'),'pdf');
    end
end
save([rootpath,'summaryAllTreat.mat'],'amplitude','Amp_mean','Amp_average','pResponsive','latencyPeak','latencyOnset');
%}
%%
close all;
cd(rootpath);
clear;
load('summaryAllTreat.mat');
celltype = {'vglut2+','vglut2-','vgat+','vgat-'};
temp_amp=cell(1,length(celltype));
for i=1:length(celltype)
temp_amp{1,i}=Amp_mean{1,i}(:,1);
end
ind_cellRes=cellfun(@(x) ~isnan(x),temp_amp,'UniformOutput',false);
cellNumRes=cellfun(@(x) sum((~isnan(x)),1),temp_amp);
cellNum=cellfun(@(x) size(x,1),temp_amp);

%plot responsive cell proportion
figPResCell=figure;
bar(cellNumRes./cellNum,'w');
set(gca,'XTickLabel',celltype);
ylabel('p(responsive cell)');
set(gca, 'Ylim',[0,1]);
for i=1:length(cellNum)
    text(i-0.2,0.1,[num2str(cellNumRes(i)),'/',num2str(cellNum(i))],'FontSize',14);
end
set(gca,'FontName','Arial','FontSize',14);
%statistics%maybe should use chi-square
responsiveCell=cell(1,4);
for i=1:length(responsiveCell)
    responsiveCell{1,i}=[ones(1,cellNumRes(i)), zeros(1,cellNum(i)-cellNumRes(i))];
end
p12=ranksum(responsiveCell{1,1},responsiveCell{1,2});
p13=ranksum(responsiveCell{1,1},responsiveCell{1,3});
p34=ranksum(responsiveCell{1,3},responsiveCell{1,4});
ylim=get(gca,'Ylim');
text(1.5,0.9*ylim(end),plabelsymbol(p12),'FontSize',14);
text(3.5,0.9*ylim(end),plabelsymbol(p34),'FontSize',14);
text(2,ylim(end),plabelsymbol(p13),'FontSize',14);
hold on;
plot([1,2],[ylim(2)*0.85,ylim(2)*0.85],'k-','LineWidth',2);
plot([3,4],[ylim(2)*0.85,ylim(2)*0.85],'k-','LineWidth',2);
plot([1,3],[ylim(2)*0.95,ylim(2)*0.95],'k-','LineWidth',2);
box off;
saveas(gcf,'p(responsive cell).fig','fig');
%
figure;
temp=cellNumRes./cellNum;
temp(2)=[];
temp(3)=[];
bar(temp,'w');
ylabel('p(responsive cell)');
set(gca,'XTick',[1,2]);
set(gca,'XTickLabel',{celltype{1,1},celltype{1,3}});
set(gca,'FontName','Arial','FontSize',14);
box off;
%
%plot EPSP amplitude,latency,p(response)
amp_aver_epsp=cell(1,length(celltype));
amp_mean_epsp=cell(1,length(celltype));
amp_mean_epsp_sig=cell(1,length(celltype));
latency_onset_epsp=cell(1,length(celltype));
latency_peak_epsp=cell(1,length(celltype));
presponse_epsp=cell(1,length(celltype));
presponse_epsp_sig=cell(1,length(celltype));
for i=1:length(celltype)
    amp_aver_epsp{1,i}=Amp_average{2,i};
    amp_mean_epsp{1,i}=Amp_mean{2,i}(:,1);
    amp_mean_epsp_sig{1,i}=Amp_mean{2,i}(:,2);
    latency_onset_epsp{1,i}=latencyOnset{2,i};
    latency_peak_epsp{1,i}=latencyPeak{2,i};
    presponse_epsp{1,i}=pResponsive{2,i}(:,1);
    presponse_epsp_sig{1,i}=pResponsive{2,i}(:,2);
end
figEPSPAmpMean=bar_error(amp_mean_epsp,[1,3],celltype,'mean EPSP amplitude(mV)'); 
figEPSPAmpMeanSig=bar_error(amp_mean_epsp_sig,[1,3],celltype,'mean EPSP amplitude (sig)'); 
figEPSPAmpAver=bar_error(amp_aver_epsp,[1,3],celltype,'average EPSP amplitude'); 
figEPSPlatencyOnset=bar_error(latency_onset_epsp,[1,3],celltype,'latency of EPSP onset(ms)');
figEPSPlatencyPeak=bar_error(latency_peak_epsp,[1,3],celltype,'latency of EPSP peak');
figEPSPpresponse=bar_error(presponse_epsp,[1,3],celltype,'p(responsive EPSP)');
figEPSPpresponseSig=bar_error(presponse_epsp_sig,[1,3],celltype,'p(responsive EPSP) (sig)');
% amp_mean_epsp_refined_by_10p1=amp_mean_epsp;
% amp_mean_epsp_sig_refined_by_10p1=amp_mean_epsp_sig;
% latency_onset_epsp_refined_by_10p1=latency_onset_epsp;
% latency_peak_epsp_refined_by_10p1=latency_peak_epsp;
% presponse_epsp_refined_by_10p1=presponse_epsp;
% presponse_epsp_sig_refined_by_10p1=presponse_epsp_sig;
% for i=1:length(celltype)
%     amp_mean_epsp_refined_by_10p1{1,i}(~ind_cellRes{1,i})=nan;
%     amp_mean_epsp_sig_refined_by_10p1{1,i}(~ind_cellRes{1,i})=nan;
%     latency_onset_epsp_refined_by_10p1{1,i}(~ind_cellRes{1,i})=nan;
%     latency_peak_epsp_refined_by_10p1{1,i}(~ind_cellRes{1,i})=nan;
%     presponse_epsp_refined_by_10p1{1,i}(~ind_cellRes{1,i})=nan;
%     presponse_epsp_sig_refined_by_10p1{1,i}(~ind_cellRes{1,i})=nan;
% end
% figEPSPAmpMeanRefine=bar_error(amp_mean_epsp_refined_by_10p1,[1,3],celltype,'refined mean EPSP amplitude');
% figEPSPAmpMeanRefineSig=bar_error(amp_mean_epsp_sig_refined_by_10p1,[1,3],celltype,'refined mean EPSP amplitude (sig)');
% figfigEPSPlatencyOnsetRefine=bar_error(latency_onset_epsp_refined_by_10p1,[1,3],celltype,'refined latency of EPSP onset');
% figEPSPlatencyPeakRefine=bar_error(latency_peak_epsp_refined_by_10p1,[1,3],celltype,'refined latency of EPSP peak');
% figEPSPpresponseRefine=bar_error(presponse_epsp_refined_by_10p1,[1,3],celltype,'refined p(responsive EPSP)');
% figEPSPpresponseRefineSig=bar_error(presponse_epsp_sig_refined_by_10p1,[1,3],celltype,'refined p(responsive EPSP) (sig)');
%save vectors
responsive_cell=cell(1,2);
mean_amplitude=cell(1,2);
mean_amplitude_sig=cell(1,2);
latency=cell(1,2);
p_responsive_trial=cell(1,2);
p_responsive_trial_sig=cell(1,2);
for i=1:2
responsive_cell{1,i}=responsiveCell{1,2*i-1};
mean_amplitude{1,i}=amp_mean_epsp_refined_by_10p1{1,2*i-1};
mean_amplitude_sig{1,i}=amp_mean_epsp_sig_refined_by_10p1{1,2*i-1};
latency{1,i}=latency_onset_epsp_refined_by_10p1{1,2*i-1};
p_responsive_trial{1,i}=presponse_epsp_refined_by_10p1{1,2*i-1};
p_responsive_trial_sig{1,i}=presponse_epsp_sig_refined_by_10p1{1,2*i-1};
end
save(['summaryVectors.mat'],'responsive_cell','mean_amplitude','latency','p_responsive_trial','mean_amplitude_sig','p_responsive_trial_sig');
figmeanAmp=histogram_compare(mean_amplitude,[],'Mean Amp distribution');
figmeanAmpSig=histogram_compare(mean_amplitude_sig,[],'Mean Amp (sig) distribution');
figpRes=histogram_compare(p_responsive_trial,[],'P(responsive) distribution');
figpResSig=histogram_compare(p_responsive_trial_sig,[],'P(responsive) (sig) distribution');
figlatency=histogram_compare(latency,[],'Latency distribution');
%% function part
function figname=histogram_compare(dataraw,cellind,titlestr)
if isempty(cellind) % [] means use all raw data
    data=dataraw;
else 
    data=cell(1,length(cellind)); %cellind decide which data to include
    for i=1:length(cellind)
        data{1,i}=dataraw{cellind(i)};           
    end
end
% figure part
color={[64,165,215]/255,[46,175,74]/255};
legendstr={'vglut2+','vgat+'};
figname=figure;
if size(data{1,1},2)==1
    range1=cellfun(@max,data);
    range2=cellfun(@min,data);
    for i=1:length(data)
        width=(max(range1)-min(range2))/20;
        histogram(data{1,i},'DisplayStyle','stairs','EdgeColor',color{i},'LineWidth',2,'BinMethod' ,'fd','BinWidth' ,width);
        hold on;
    end
    hl=legend(legendstr,'Location','Best');
    set(hl,'Box','off');
end
title(titlestr);
box off;
set(gca,'FontName','Arial','FontSize',14);
saveas(gcf,[titlestr,'.fig'],'fig');
saveas(gcf,[titlestr,'.png'],'png');
end

function figname=bar_error(dataraw,cellind,xlabelraw,ylabelstr)
if isempty(cellind) % [] means use all raw data
    data=dataraw;
    xlabel=xlabelraw;
else 
    data=cell(1,length(cellind)); %cellind decide which data to include
    xlabel=cell(1,length(cellind));
    for i=1:length(cellind)
        data{1,i}=dataraw{cellind(i)};      
        xlabel{1,i}=xlabelraw{cellind(i)};      
    end
end
% figure part
figname=figure;
if size(data{1,1},2)==1
    y_pres=cellfun(@(x) nanmean(x),data);
    sem_pres=cellfun(@(x) nanstd(x)/sqrt(length(x)),data);
else
    y_pres=cellfun(@(x) nanmean(x),data,'UniformOutput',false);
    sem_pres=cellfun(@(x) nanstd(x)/sqrt(length(x)),data,'UniformOutput',false);
    nrow=length(y_pres);
    y_pres=cell2mat(y_pres);
    sem_pres=cell2mat(sem_pres);
    y_pres=reshape(y_pres,[],nrow);
    sem_pres=reshape(sem_pres,nrow,[])';
end

if isinf(sum(y_pres))
    warning(strcat('Inf found in',ylabelstr));
end
% c=categorical(celltype);
% bar(c,y_pres,'w');
bfig=bar(y_pres,'FaceColor','flat');
hold on;
if size(y_pres,1)==1
    bfig.FaceColor=[1,1,1];
    set(gca,'XTickLabel',xlabel);
else
    bfig(1).FaceColor=[1,1,1];
    bfig(2).FaceColor=[0,0,0];%second category use another color
    legend(bfig,xlabel,'AutoUpdate','off');
    legend('boxoff');
    set(gca,'XTickLabel',{'50ms','100ms','200ms'});
end

if size(data{1,1},2)==1
    for i=1:1:length(data)
        scatter(i*ones(1,length(data{1,i})),data{1,i},20,[0.5,0.5,0.5],'filled');
    end
errorbar(1:length(data),y_pres,sem_pres,'ok','LineWidth',2);    
else
    for i=1:size(data{1,1},2)
        for j=1:length(data)
            scatter(i*ones(1,size(data{1,j},1))+(j-size(data{1,1},2)/2)*0.3,data{1,j}(:,i),20,[0.5,0.5,0.5],'filled');
            errorbar(i+(j-size(data{1,1},2)/2)*0.3,y_pres(i,j),sem_pres(i,j),'ok','LineWidth',2);
        end
    end
end

ylabel(ylabelstr);
set(gca,'FontName','Arial','FontSize',14);
box off;
ylim=get(gca,'Ylim');
% statistics part
if size(data{1,1},2)==1
    if length(data)==4
        p12=ranksum(data{1,1},data{1,2});
        p34=ranksum(data{1,3},data{1,4});
        p13=ranksum(data{1,1},data{1,3});
        % p12=ttest2(data{1,1},data{1,2});
        % p34=ttest2(data{1,3},data{1,4});
        % p13=ttest2(data{1,1},data{1,3});
        text(1.5,0.9*ylim(end),plabelsymbol(p12),'FontSize',14);
        text(3.5,0.9*ylim(end),plabelsymbol(p34),'FontSize',14);
        text(2,ylim(end),plabelsymbol(p13),'FontSize',14);
        plot([1,2],[ylim(2)*0.85,ylim(2)*0.85],'k-','LineWidth',2);
        plot([3,4],[ylim(2)*0.85,ylim(2)*0.85],'k-','LineWidth',2);
        plot([1,3],[ylim(2)*0.95,ylim(2)*0.95],'k-','LineWidth',2);
    elseif length(data)==2
        p12=ranksum(data{1,1},data{1,2});
        % p12=ttest2(data{1,1},data{1,2});
        text(1.5,0.9*ylim(end),plabelsymbol(p12),'FontSize',14);
        plot([1,2],[ylim(2)*0.85,ylim(2)*0.85],'k-','LineWidth',2);
    end
else
    if length(data)==2
        for i=1:size(data{1,1},2)
            p(i)=ranksum(data{1,1}(i,:),data{1,2}(i,:));
            text(i,0.9*ylim(end),plabelsymbol(p(i)),'FontSize',14);
            plot([-0.2,0.2]+i,[ylim(2)*0.85,ylim(2)*0.85],'k-','LineWidth',2);
        end
    end
end
saveas(gcf,[ylabelstr,'.fig'],'fig');
saveas(gcf,[ylabelstr,'.png'],'png');
end