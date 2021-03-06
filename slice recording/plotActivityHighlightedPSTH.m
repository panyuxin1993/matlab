rootpath='D:\xulab\project\slice recording\data3\EPSP summary\';
celltype = {'vglut2+','vglut2-','vgat+','vgat-'};
Amp_mean=cell(1,length(celltype));
Amp_average=cell(1,length(celltype));
Amp_ttest=cell(1,length(celltype));
pResponsive=cell(1,length(celltype));
SNR=cell(1,length(celltype));
SNR_neuron_mean=cell(1,length(celltype));
SNR_neuron_average=cell(1,length(celltype));
latency=cell(2,length(celltype));%summary of response latency from stim_onset; the 1st row is latency of reponse onset, the  2nd row is latency of response peak
latency1=cell(1,length(celltype));%1st pulse latency
colorCases={[0.7,0.7,0.7],[1,0.9,0.7]};%colorCase{2} highlighted
colorMean={[0,0,0],[1,0.5,0],[1,0,0]};
times=3;
savepath=['D:\xulab\project\slice recording\data3\EPSP summary\SNR',num2str(times),'STD\'];
for n_celltype=1:length(celltype)
    AbfFiles = dir(strcat(rootpath,celltype{n_celltype},'\*.abf'));
    cd(rootpath);    
    amplitude=cell(1,length(AbfFiles));
    Amp_mean{1,n_celltype}=zeros(1,length(AbfFiles));
    Amp_average{1,n_celltype}=zeros(1,length(AbfFiles));
    Amp_ttest{1,n_celltype}=zeros(1,length(AbfFiles));%determine responsive by ttest of amp and basline across trial
    pResponsive{1,n_celltype}=zeros(1,length(AbfFiles));
    SNR{1,n_celltype}=[]; %row vector storing each trial signal-to-noise ratio(response amplitude/ std_baseline)
    SNR_neuron_mean{1,n_celltype}=[];
    SNR_neuron_average{1,n_celltype}=[];
    latency{1,n_celltype}=zeros(length(AbfFiles),1);%latency of response onset
    latency{2,n_celltype}=zeros(length(AbfFiles),1);%latency of response peak
    for n_abf=1:length(AbfFiles)
%             strcat(cellFolder,AbfFiles(n_abf).name);
        [d,si,h]=abfload(strcat(rootpath,celltype{n_celltype},'\',AbfFiles(n_abf).name),'channels','a');
        ts=(1:h.sweepLengthInPts)*h.si/10^6;
        sr=10^6/si;
        xlim=[1,3];
        nchannel=size(d,2);
        amplitude{1,n_abf}=zeros(1,size(d,3)+1);% vector of each individual amplitudes,last one is mean trace peak
        
        figraw=figure;
        for nsubplot=1:nchannel
            figure(figraw);
            subplot(nchannel,1,nsubplot);
            sumStimTrial=0;
            tdata=zeros(size(d,3),2);
            for i=1:size(d,3)
                laser_duration=find(d(:,2,i)>1);%second channel to find when laser on
                if isempty(laser_duration) %if some trials no stimuli, then pass this trial
                    continue;
                else
                    sumStimTrial=sumStimTrial+1;
                end
                laser_start=laser_duration(1);
                STD_baseline = std(md(laser_start(1)-0.5*sr:laser_start(1),1));%0.5s before laser start as baseline
                baseline=mean(md(laser_start(1)-0.5*sr:laser_start(1),1));
                peak=max(d(laser_start:laser_start+0.05*sr,1,i)); %50ms within after stimulus
                amp=abs(peak-baseline);
                if amp<=STD_baseline*times
                    plot(ts,d(:,nsubplot,i),'-','Color',colorCases{1},'LineWidth',0.5);
                    hold on;
                elseif amp>STD_baseline*times
                    amplitude{1,n_abf}(1,i)=amp;
                    plot(ts,d(:,nsubplot,i),'-','Color',colorCases{2},'LineWidth',0.5);
                    hold on;
                end
                if nsubplot==1
                    SNR{1,n_celltype}=[SNR{1,n_celltype} amp/STD_baseline];
                end
                tdata(i,1)=baseline;
                tdata(i,2)=peak;
            end
            if nsubplot==1
                SNR_neuron_average{1,n_celltype}=[SNR_neuron_average{1,n_celltype} mean(SNR{1,n_celltype}(length(SNR{1,n_celltype})-size(d,3)+1:end))];
            end
            figure(figraw);
            md=mean(d,3);
            laser_duration=find(md(:,2)>3);%second channel to find when laser on
            laser_start=laser_duration(1);
            STD_baseline = std(md(:,1));%std(md(laser_start-2*sr:laser_start,1));
            baseline=mean(md(:,1));%mean(md(laser_start-2*sr:laser_start,1));
            peak=max(md(laser_start:laser_start+0.05*sr,1)); %0.05s within after stimulus
            peakT=find(md(laser_start:laser_start+0.05*sr,1)==peak);
            onsetT=find(abs(md(laser_start:laser_start+0.05*sr,1)-baseline)>STD_baseline*times);
            amp=abs(peak-baseline);
            latency{2,n_celltype}(n_abf,1)=peakT(1)/sr*1000;% ms from stim_onset
            if isempty(onsetT)
                latency{1,n_celltype}(n_abf,1)=nan;
            else
                latency{1,n_celltype}(n_abf,1)=onsetT(1)/sr*1000;
            end
            if amp<=STD_baseline*times
                amplitude{1,n_abf}(1,end)=nan;
                plot(ts,md(:,nsubplot),'-','Color',colorMean{1},'LineWidth',2);
                hold on;
            elseif amp>STD_baseline*times
                amplitude{1,n_abf}(1,end)=amp;
                plot(ts,md(:,nsubplot),'-','Color',colorMean{2},'LineWidth',2);
                hold on;
            end
            if nsubplot==1
                peakT=(latency{2,n_celltype}(n_abf,:)*sr/1000+laser_start)/sr;
                peakT(isnan(peakT)==1)=[];
                peaks=md(round(peakT*sr),nsubplot);
                scatter(peakT,peaks,10,'k');
                onsetT=(latency{1,n_celltype}(n_abf,:)*sr/1000+laser_start)/sr;
                onsetT(isnan(onsetT)==1)=[];
                onsetAmp=md(round(onsetT*sr),nsubplot);
                scatter(onsetT,onsetAmp,10,'k','filled');%magenta filled dots as onsets
            end
            if nsubplot==1
                SNR_neuron_mean{1,n_celltype}=[SNR_neuron_mean{1,n_celltype} amp/STD_baseline];
            end
            sum_responsive=sum(double(amplitude{1,n_abf}(1:end-1)>0)); %last one is mean trace peak
            pro_responsive=sum_responsive/sumStimTrial;%rather than assuming that all trials have activity
            amplitude_meantrace=amplitude{1,n_abf}(end);
            if sum(amplitude{1,n_abf}(1:end-1)>0)==0
                amplitude_average=nan;
            else
                amplitude_average=sum(amplitude{1,n_abf}(1,amplitude{1,n_abf}(1:end-1)>0))/sum(amplitude{1,n_abf}(1:end-1)>0); %exclude no response trials and only average among those responsive
            end
            Amp_mean{1,n_celltype}(1,n_abf)=amplitude_meantrace;
            Amp_average{1,n_celltype}(1,n_abf)=amplitude_average;
            pResponsive{1,n_celltype}(1,n_abf)=pro_responsive;
            [~,p]=ttest(tdata(:,1),tdata(:,2));
%            %following plot use ttest significance
%             if p<0.01
%                 Amp_ttest{1,n_celltype}(1,n_abf)=amp;
%                 plot(ts,md(:,nsubplot),'-','Color',colorMean{3},'LineWidth',2);
%                 hold on;
%             else
%                 Amp_ttest{1,n_celltype}(1,n_abf)=nan;
%                 plot(ts,md(:,nsubplot),'-','Color',colorMean{1},'LineWidth',2);
%                 hold on;
%             end      
            strsummary=['pResp=',num2str(pro_responsive),';ampMean=',num2str(amplitude_meantrace),';ampAve',num2str(amplitude_average),'SNRCell',num2str(SNR_neuron_mean{1,n_celltype}(n_abf))];
            set(gca,'xlim',xlim);
            xlabel('time(s)');
            set(gca,'FontName','Arial','FontSize',14);
            xlim=get(gca,'xlim');
            ylim=get(gca,'ylim');
            text(xlim(1),ylim(end),strsummary);
            box off;
            saveas(figraw,[savepath,celltype{n_celltype},'-',AbfFiles(n_abf).name,'.fig'],'fig');
        end
        close all;
    end
end
%
save([savepath,'summary.mat'],'Amp_mean','Amp_average','pResponsive','latency','SNR');

%%

close all;
cd(savepath);
load('summary.mat');
celltype = {'vglut2+','vglut2-','vgat+','vgat-'};
cellNum=cellfun(@(x) length(x),pResponsive);
% cellNumRes=cellfun(@(x) length(x)-sum(x==0),pResponsive);
cellNumRes=cellfun(@(x) length(x)-sum(isnan(x)),Amp_mean);
% cellNumRes=cellfun(@(x) length(x)-sum(isnan(x)),Amp_ttest);

figSNR=figure;%plot signal-to-noise ratio
for i=1:length(SNR)
    subplot(2,2,i)
    histogram(SNR{1,i});
    xlabel('signal to noise ratio of neuronal response to light');
    ylabel('number of trials');
    title(celltype{i});
%     set(gca,'xlim',[0,10]);
    legend(['cell number ',num2str(cellNum(i))]);
%     set(gca,'FontName','Arial','FontSize',14);
end
box off;
saveas(figSNR,'SNR of neuronal response.fig','fig');

%summary of latency 
figlatency=figure;
for i=1:length(latency)
    subplot(2,4,i);
    histogram(latency{1,i});
    xlabel('latency of response onset(ms)');
    ylabel('number of mean responses');
    title(celltype{i});
%     set(gca,'xlim',[0,10]);
%     set(gca,'FontName','Arial','FontSize',14);
    subplot(2,4,i+4);
    histogram(latency{2,i});
    xlabel('latency of response peak(ms)');
    ylabel('number of mean responses');
    title(celltype{i});
end
box off;
saveas(figlatency,'latency of response.fig','fig');

figPResCell=figure;%plot responsive cell proportion
bar(cellNumRes./cellNum,'w');
set(gca,'XTickLabel',celltype);
ylabel('p(responsive cell)');
set(gca, 'Ylim',[0,1]);
for i=1:length(cellNum)
    text(i-0.2,0.1,[num2str(cellNumRes(i)),'/',num2str(cellNum(i))],'FontSize',14);
end
set(gca,'FontName','Arial','FontSize',14);
%statistics
demo=cell(1,4);
for i=1:length(demo)
    demo{1,i}=[ones(1,cellNumRes(i)), zeros(1,cellNum(i)-cellNumRes(i))];
end
p12=ranksum(demo{1,1},demo{1,2});
p13=ranksum(demo{1,1},demo{1,3});
p34=ranksum(demo{1,3},demo{1,4});
ylim=get(gca,'Ylim');
text(1.5,0.9*ylim(end),plabelsymbol(p12),'FontSize',14);
text(3.5,0.9*ylim(end),plabelsymbol(p34),'FontSize',14);
text(2,ylim(end),plabelsymbol(p13),'FontSize',14);
box off;
saveas(gcf,'p(responsive cell).fig','fig');
%plot responive trial proportion 
pResCorrected=pResponsive;%plot only 'responsive cell'
for i=1:length(pResponsive)
    pResCorrected{1,i}(pResCorrected{1,i}==0)=nan;
end
figPRes_nonzero=bar_error(pResCorrected,[1,3],celltype,'p(respond to light stimuli) if cell ever respond');
pResMean=pResponsive;
for i=1:length(pResponsive)
    pResMean{1,i}(isnan(Amp_mean{1,i}))=nan;
end
figResMean=bar_error(pResMean,[1,3],celltype,'p(respond to light stimuli) if mean trace respond');

%plot two kind of amplitude
Amp_averageCor=Amp_average;
for i=1:length(Amp_averageCor)
    %Amp_averageCor{1,i}(isnan(Amp_mean{1,i}))=nan;%Average amplitude respond to light stimuli if mean trace respond
    Amp_averageCor{1,i}(pResponsive{1,i}==0)=nan;%Average amplitude respond to light stimuli if ever respond
end
figAmp_averge=bar_error(Amp_averageCor,[1,3],celltype,'Average amplitude(mV) respond to light stimuli if ever respond');
%figAmp_averge=bar_error(Amp_averageCor,celltype,'Average amplitude respond to light stimuli if mean trace respond');
Amp_meanCor=Amp_mean;
for i=1:length(Amp_meanCor)
    %Amp_meanCor{1,i}(isnan(Amp_mean{1,i}))=nan;%Mean amplitude respond to light stimuli if mean trace respond
    Amp_meanCor{1,i}(pResponsive{1,i}==0)=nan;%Mean amplitude respond to light stimuli if ever respond
end
% figAmp_mean=bar_error(Amp_meanCor,celltype,'Mean amplitude respond to light stimuli if ever respond');
figAmp_mean=bar_error(Amp_mean,[1,3],celltype,'Mean amplitude(mV) respond to light stimuli if mean trace respond');
% figAmp_ttest=bar_error(Amp_ttest,celltype,'Amplitude respond to light stimuli if ttest significant');

%SNR
% figSNR_neuron_mean=bar_error(SNR_neuron_mean,celltype,'Mean SNR of neuron');%plot SNR of neurons
% SNR_neuron_mean_cor=SNR_neuron_mean;
% for i=1:length(SNR_neuron_mean_cor)
%     SNR_neuron_mean_cor{1,i}(SNR_neuron_mean_cor{1,i}<3)=nan;
% end
% figSNR_neuron_mean_cor=bar_error(SNR_neuron_mean_cor,celltype,'Mean SNR of neuron(show significant)');%plot SNR of neurons
% figSNR_neuron_average=bar_error(SNR_neuron_average,celltype,'Average SNR of neuron');%plot SNR of neurons
% SNR_neuron_average_cor=SNR_neuron_average;
% for i=1:length(SNR_neuron_average_cor)
%     SNR_neuron_average_cor{1,i}(SNR_neuron_average_cor{1,i}<3)=nan;
% end
% figSNR_neuron_average_cor=bar_error(SNR_neuron_average_cor,celltype,'Average SNR of neuron(show significant)');%plot SNR of neurons

%plot corrlation of different parameters
figCorr=figure;
set(gcf,'Position',[0,0,1200,400]);
subplot(1,3,1);
for i=1:length(pResponsive)
    scatter(pResponsive{1,i},Amp_average{1,i});
    hold on;
end
legend(celltype);
xlabel('p(respond to light stimuli)')
ylabel('Average amplitude(mV) respond to light stimuli');
set(gca,'FontName','Arial','FontSize',14);
subplot(1,3,2);
for i=1:length(pResponsive)
    scatter(pResponsive{1,i},Amp_mean{1,i});
    hold on;
end
legend(celltype);
xlabel('p(respond to light stimuli)')
ylabel('Mean amplitude(mV) respond to light stimuli');
set(gca,'FontName','Arial','FontSize',14);
subplot(1,3,3);
for i=1:length(Amp_mean)
    scatter(Amp_mean{1,i},Amp_average{1,i});
    hold on;
end
legend(celltype);
xlabel('Mean amplitude(mV) respond to light stimuli')
ylabel('Average amplitude(mV) respond to light stimuli');
set(gca,'FontName','Arial','FontSize',14);
box off;
saveas(gcf,'Correlation of parameters.fig','fig');

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
y_pres=cellfun(@(x) nanmean(x),data);
sem_pres=cellfun(@(x) nanstd(x)/sqrt(length(x)),data);
if isinf(sum(y_pres))
    warning(strcat('Inf found in',ylabelstr));
end
% c=categorical(celltype);
% bar(c,y_pres,'w');
bar(y_pres,'w');
hold on;
for i=1:length(data)
    scatter(i*ones(1,length(data{1,i})),data{1,i},20,[0.5,0.5,0.5],'filled');
end
errorbar(y_pres,sem_pres,'ok','LineWidth',2);
set(gca,'XTickLabel',xlabel);
ylabel(ylabelstr);
set(gca,'FontName','Arial','FontSize',14);
box off;
ylim=get(gca,'Ylim');
% statistics partj
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
saveas(gcf,[ylabelstr,'.fig'],'fig');
saveas(gcf,[ylabelstr,'.png'],'png');
end

% function figname=bar_error(data,xlabel,ylabelstr)
% % statistics part
% if sum(~isnan(data{1,2}))~=0 && sum(~isnan(data{1,4}))~=0 %this means No data remaining after removal of NaNs. and can not test
%     p12=ranksum(data{1,1},data{1,2});
%     p34=ranksum(data{1,3},data{1,4});
% else
%     p12=1;
%     p34=1;
% end
% if sum(isnan(data{1,1}))~=0 || sum(isnan(data{1,3}))~=0 %this means No data remaining after removal of NaNs. and can not test
%     p13=ranksum(data{1,1},data{1,3});
% else
%     p13=1;
% end
% % p12=ttest2(data{1,1},data{1,2});
% % p34=ttest2(data{1,3},data{1,4});
% % p13=ttest2(data{1,1},data{1,3});
% % figure part
% figname=figure;
% y_pres=cellfun(@(x) nanmean(x),data);
% sem_pres=cellfun(@(x) nanstd(x)/sqrt(length(x)),data);
% % c=categorical(celltype);
% % bar(c,y_pres,'w');
% bar(y_pres,'w');
% hold on;
% for i=1:length(data)
%     scatter(i*ones(1,length(data{1,i})),data{1,i},20,[0.5,0.5,0.5],'filled');
% end
% errorbar(y_pres,sem_pres,'ok','LineWidth',2);
% set(gca,'XTickLabel',xlabel);
% ylabel(ylabelstr);
% set(gca,'FontName','Arial','FontSize',14);
% box off;
% ylim=get(gca,'Ylim');
% text(1.5,0.9*ylim(end),plabelsymbol(p12),'FontSize',14);
% text(3.5,0.9*ylim(end),plabelsymbol(p34),'FontSize',14);
% text(2,ylim(end),plabelsymbol(p13),'FontSize',14);
% saveas(gcf,[ylabelstr,'.fig'],'fig');
% end

function [str]=plabelsymbol(pvalue)
if pvalue<0.05 && pvalue>=0.01
    str=' *';
elseif pvalue<0.01 && pvalue>=0.001
    str=' **';
elseif pvalue<0.001
    str=' ***';
elseif pvalue>=0.05 && pvalue<0.01
    pvalue=round(pvalue,2);%保留两位即可
    str=strcat('p=',num2str(pvalue));
else
    pvalue=round(pvalue,1);%保留1位即可
    str=strcat('p=',num2str(pvalue));
end
end