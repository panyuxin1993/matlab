function [ fig ] = fPlotEphysSummary( fig,position )
%FPLOTEPHYSSUMMARY plot ephys results, including mean amplitude, latency
rootpath='D:\xulab\project\slice recording\data3\';
load(strcat(rootpath,'summaryAllTreat.mat'));
celltype = {'vglut2+','vglut2-','vgat+','vgat-'};
temp_amp=cell(1,length(celltype));
for i=1:length(celltype)
temp_amp{1,i}=Amp_mean{1,i}(:,1);
end
ind_cellRes=cellfun(@(x) ~isnan(x),temp_amp,'UniformOutput',false);
cellNumRes=cellfun(@(x) sum((~isnan(x)),1),temp_amp);
cellNum=cellfun(@(x) size(x,1),temp_amp);
figure(fig);
%plot responsive cell proportion
subplot(2,2,1);
temp=cellNumRes./cellNum;
temp(2)=[];
temp(3)=[];
bar(temp,'w');
ylabel('p(responsive cell)');
set(gca,'XTick',[1,2]);
set(gca,'XTickLabel',{celltype{1,1},celltype{1,3}});
set(gca,'FontName','Arial','FontSize',14);
box off;
set(gca,'Xlim',[0,3]);
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
% figEPSPAmpMean=bar_error(amp_mean_epsp,[1,3],celltype,'mean EPSP amplitude(mV)'); 
% figEPSPAmpMeanSig=bar_error(amp_mean_epsp_sig,[1,3],celltype,'mean EPSP amplitude (sig)'); 
% figEPSPAmpAver=bar_error(amp_aver_epsp,[1,3],celltype,'average EPSP amplitude'); 
% figEPSPlatencyOnset=bar_error(latency_onset_epsp,[1,3],celltype,'latency of EPSP onset(ms)');
% figEPSPlatencyPeak=bar_error(latency_peak_epsp,[1,3],celltype,'latency of EPSP peak');
% figEPSPpresponse=bar_error(presponse_epsp,[1,3],celltype,'p(responsive EPSP)');
% figEPSPpresponseSig=bar_error(presponse_epsp_sig,[1,3],celltype,'p(responsive EPSP) (sig)');
subplot(2,2,2);
bar_error(presponse_epsp,[1,3],celltype,'P(responsive)');
subplot(2,2,3);
bar_error(latency_onset_epsp,[1,3],celltype,'Onset Latency(ms)'); 
subplot(2,2,4);
bar_error(amp_mean_epsp,[1,3],celltype,'Mean Amplitude(mV)'); 
end

%% function part

function []=bar_error(dataraw,cellind,xlabelraw,ylabelstr)
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
% figname=figure;
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
        text(1.5,1.2*ylim(end),plabelsymbol(p12),'FontSize',14);
        plot([1,2],[ylim(2)*1.1,ylim(2)*1.1],'k-','LineWidth',2);
    end
else
    if length(data)==2
        for i=1:size(data{1,1},2)
            p(i)=ranksum(data{1,1}(i,:),data{1,2}(i,:));
            text(i,1.2*ylim(end),plabelsymbol(p(i)),'FontSize',14);
            plot([-0.2,0.2]+i,[ylim(2)*1.05,ylim(2)*1.05],'k-','LineWidth',2);
        end
    end
end
set(gca,'Xlim',[0,3]);
set(gca,'Ylim',[ylim(1),ylim(end)*1.2]);
% saveas(gcf,[ylabelstr,'.fig'],'fig');
% saveas(gcf,[ylabelstr,'.png'],'png');
end


