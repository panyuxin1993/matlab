%first decide which cell respond to stimuli using amplitude>3 std baseline,
%then no response cells' amplitude is nan and responsive cells' amplitude
%is |peak-baseline| no matter whether its SNR>3 
rootpath='D:\xulab\project\slice recording\data3\PPR summary\';
celltype = {'vglut2+','vglut2-','vgat+','vgat-'};
Amp=cell(1,length(celltype)); %response amp for each pulse of each trials
Amp_mean=cell(1,length(celltype));
Amp_mean_sig=cell(1,length(celltype));
SNR=cell(1,length(celltype));
colorCases={[1,0.7,0.7],[0.7,1,0.7],[0.7,0.7,1],[0.7,0.7,0.7]};%colorCase{1}-{4} increase, decrease, maintain/contradict, no response
colorMean={[1,0,0],[0,1,0],[0,0,1],[0,0,0]};
PPR=cell(1,length(celltype));%only compare the first 2 response, ratio >1 or <1
slope=cell(1,length(celltype));%linearly fit the whole trend of response change, slope positive or negative
plasticity=cell(1,length(celltype));%summary results from PPR and slope
latency=cell(2,length(celltype));%summary of response latency from stim_onset; the 1st row is latency of reponse onset, the  2nd row is latency of response peak
latency1=cell(1,length(celltype));%1st pulse latency
times=3;
savepathPPR=['D:\xulab\project\slice recording\XUNL\summary of short term plasticity\PPR\'];
savepathRefine=['D:\xulab\project\slice recording\data3\PPR summary\PPR\'];
for n_celltype=1:length(celltype)
    AbfFiles = dir(strcat(rootpath,celltype{n_celltype},'\*.abf'));
    cd(rootpath);    
    Amp{1,n_celltype}=cell(1,length(AbfFiles));%each cell a matrix containing individual amp
    Amp_mean{1,n_celltype}=zeros(length(AbfFiles),10);%2d matrix, 1d-cell,2d 10 pulse elicit response
    Amp_mean_sig{1,n_celltype}=zeros(length(AbfFiles),10);
    PPR{1,n_celltype}=zeros(length(AbfFiles),1);
    slope{1,n_celltype}=zeros(length(AbfFiles),1);
    plasticity{1,n_celltype}=zeros(length(AbfFiles),1);%1-increase,2-decrease,3-contradict,4-no response
    SNR{1,n_celltype}=[]; %row vector storing each trial signal-to-noise ratio(response amplitude/ std_baseline)
    latency{1,n_celltype}=zeros(length(AbfFiles),10);%latency of response onset
    latency{2,n_celltype}=zeros(length(AbfFiles),10);%latency of response peak
    latency1{1,n_celltype}=[];
    for n_abf=1:length(AbfFiles)
%             strcat(cellFolder,AbfFiles(n_abf).name);
        [d,si,h]=abfload(strcat(rootpath,celltype{n_celltype},'\',AbfFiles(n_abf).name),'channels','a');
        ts=(1:h.sweepLengthInPts)*h.si/10^6;
        sr=10^6/si;
        xlim=[1.5,3.5];
        nchannel=size(d,2);
        ntrial=size(d,3);
        md=mean(d,3);
        laser_duration=find(md(:,2)>3);%second channel to find when laser on
        laser_start=laser_duration(1):0.1*sr:(laser_duration(1)+sr*0.9);%first 10 pulses
        STD_baseline = std(md(laser_start(1)-1*sr:laser_start(1),1));%1s before laser start as baseline
        baseline=mean(md(laser_start(1)-1*sr:laser_start(1),1));
        for npulse=1:length(laser_start)
            peak=max(md(laser_start(npulse):laser_start(npulse)+0.05*sr,1)); %0.05s within after stimulus
            amp=abs(peak-baseline);
            if amp<=STD_baseline*times
                Amp_mean_sig{1,n_celltype}(n_abf,npulse)=0;
            elseif amp>STD_baseline*times
                Amp_mean_sig{1,n_celltype}(n_abf,npulse)=amp;
            end
        end
        if sum(Amp_mean_sig{1,n_celltype}(n_abf,:)==0)>8 % not respond if more than 8 (in 10) stimuli do not respond
            slope{1,n_celltype}(n_abf,1)=nan;
            PPR{1,n_celltype}(n_abf,1)=nan;
            Amp_mean{1,n_celltype}(n_abf,:)=nan;
            latency{1,n_celltype}(n_abf,:)=nan;
            latency{2,n_celltype}(n_abf,:)=nan;
            Amp{1,n_celltype}{1,n_abf}=zeros(1,10);
        else %for responsive cell, re-calculate responses
            Amp{1,n_celltype}{1,n_abf}=zeros(ntrial,10);%each row a trial, each col a pulse
            for itrial=1:ntrial
                for npulse=1:length(laser_start)
                    ibaseline=mean(d(laser_start(1)-1*sr:laser_start(1),1,itrial));%1s before laser start as baseline
                    ipeak=max(d(laser_start(npulse):laser_start(npulse)+0.05*sr,1,itrial)); %0.05s within after stimulus
                    Amp{1,n_celltype}{1,n_abf}(itrial,npulse)=abs(ipeak-ibaseline);
                end
            end
            for npulse=1:length(laser_start)
                peak=max(md(laser_start(npulse):laser_start(npulse)+0.05*sr,1)); %0.05s within after stimulus
                peakT=find(md(laser_start(npulse):laser_start(npulse)+0.05*sr,1)==peak);
                onsetT=find(abs(md(laser_start(npulse):laser_start(npulse)+0.05*sr,1)-baseline)>STD_baseline*times);
                Amp_mean{1,n_celltype}(n_abf,npulse)=abs(peak-baseline);
                latency{2,n_celltype}(n_abf,npulse)=peakT(1)/sr*1000;% ms from stim_onset
                if isempty(onsetT)
                    latency{1,n_celltype}(n_abf,npulse)=nan;
                else
                    latency{1,n_celltype}(n_abf,npulse)=onsetT(1)/sr*1000;
                end
            end
            p=polyfit(1:length(laser_start),Amp_mean{1,n_celltype}(n_abf,:),1);%linear fitting,y=p(1)x+p(2)
            slope{1,n_celltype}(n_abf,1)=p(1);
            PPR{1,n_celltype}(n_abf,1)=Amp_mean{1,n_celltype}(n_abf,2)/Amp_mean{1,n_celltype}(n_abf,1);
            if slope{1,n_celltype}(n_abf,1)*PPR{1,n_celltype}(n_abf,1)==0 ||isnan(slope{1,n_celltype}(n_abf,1)*PPR{1,n_celltype}(n_abf,1))
                warning(strcat('Amp=',num2str(Amp_mean{1,n_celltype}),',PPR=',num2str(PPR{1,n_celltype}),',slope=',num2str(slope{1,n_celltype})));
            end
        end       
        temp=latency{1,n_celltype}(n_abf,:);
        temp(isnan(temp)==1)=[];      
        if ~isempty(temp)
            if temp(1)>10%refined by latency
                slope{1,n_celltype}(n_abf,1)=nan;
                PPR{1,n_celltype}(n_abf,1)=nan;
                Amp_mean{1,n_celltype}(n_abf,:)=nan;
                latency{1,n_celltype}(n_abf,:)=nan;
                latency{2,n_celltype}(n_abf,:)=nan;
            else
                latency1{1,n_celltype}=[latency1{1,n_celltype} temp(1)];
            end     
        end
        
        figraw=figure;
        for nsubplot=1:nchannel
            figure(figraw);
            subplot(2,1,nsubplot);
            if slope{1,n_celltype}(n_abf,1)>0 && (PPR{1,n_celltype}(n_abf,1)>1 || isnan(PPR{1,n_celltype}(n_abf,1)))%2 methods same results; or when PPR is nan, only decided by slope
                plot(ts,md(:,nsubplot),'-','Color',colorMean{1},'LineWidth',2); 
                plasticity{1,n_celltype}(n_abf,1)=1;
            elseif slope{1,n_celltype}(n_abf,1)<0 && (PPR{1,n_celltype}(n_abf,1)<1 || isnan(PPR{1,n_celltype}(n_abf,1)))
                plot(ts,md(:,nsubplot),'-','Color',colorMean{2},'LineWidth',2);
                plasticity{1,n_celltype}(n_abf,1)=2;
            elseif slope{1,n_celltype}(n_abf,1)*(PPR{1,n_celltype}(n_abf,1)-1)<0 %2 methods contradict results
                plot(ts,md(:,nsubplot),'-','Color',colorMean{3},'LineWidth',2);
                plasticity{1,n_celltype}(n_abf,1)=3;
            else% no response
                plot(ts,md(:,nsubplot),'-','Color',colorMean{4},'LineWidth',2);
                plasticity{1,n_celltype}(n_abf,1)=4;
            end
            hold on;
            if nsubplot==1
                plot([ts(1),ts(end)],[baseline,baseline],'-k','LineWidth',0.1);
                plot([ts(1),ts(end)],[baseline+STD_baseline*times,baseline+STD_baseline*times],'--k','LineWidth',0.1);
                plot([ts(1),ts(end)],[baseline+STD_baseline,baseline+STD_baseline],'-.k','LineWidth',0.1);
                plot([ts(1),ts(end)],[baseline-STD_baseline*times,baseline-STD_baseline*times],'--k','LineWidth',0.1);
                plot([ts(1),ts(end)],[baseline-STD_baseline,baseline-STD_baseline],'-.k','LineWidth',0.1);
                peakT=(latency{2,n_celltype}(n_abf,:)*sr/1000+laser_start)/sr;
                peakT(isnan(peakT)==1)=[];
                peaks=md(round(peakT*sr),nsubplot);
                scatter(peakT,peaks,10,'k');
                onsetT=(latency{1,n_celltype}(n_abf,:)*sr/1000+laser_start)/sr;
                onsetT(isnan(onsetT)==1)=[];
                onsetAmp=md(round(onsetT*sr),nsubplot);
                scatter(onsetT,onsetAmp,10,'k','filled');%magenta filled dots as onsets
            end
            strsummary=['PPR=',num2str(PPR{1,n_celltype}(n_abf,1)),';slope=',num2str(slope{1,n_celltype}(n_abf,1))];
            set(gca,'xlim',xlim);
            xlabel('time(s)');
            set(gca,'FontName','Arial','FontSize',14);
            xlim=get(gca,'xlim');
            ylim=get(gca,'ylim');
            x_resp=laser_start(Amp_mean_sig{1,n_celltype}(n_abf,:)>0)/sr;
            ind_resp=find(Amp_mean_sig{1,n_celltype}(n_abf,:)>0);
            if nsubplot==1
                if ~isempty(x_resp)
                    plot(x_resp,(ylim(end))*ones(length(x_resp),1),'*k');
                    for iPulse=1:length(x_resp)
                        text(x_resp(iPulse),ylim(end)+0.5,[num2str(round(Amp_mean_sig{1,n_celltype}(n_abf,ind_resp(iPulse)),2,'significant')),',',num2str(round(latency{2,n_celltype}(n_abf,ind_resp(iPulse)),2,'significant'))]);
                    end
                end
            else
                text(xlim(1),ylim(end),strsummary);
            end
            box off;
        end
        saveas(figraw,[savepathRefine,celltype{n_celltype},'-',AbfFiles(n_abf).name,'.fig'],'fig');
        close all;
    end
end
save([savepathRefine,'summary.mat'],'Amp','Amp_mean','PPR','slope','plasticity','Amp_mean_sig','latency','latency1');
%% for calculation demo
% a=Amp_mean{1,1}(6,:);
% scatter(1:10,a,30,'K');
% hold on;
% set(gca,'FontName','Arial','FontSize',14);
% p=polyfit(1:10,a,1);
% x=1:10;
% y=p(1)*x+p(2);
% plot(x,y,'k-');
% xlabel('Pulse number');
% ylabel('Response peak(mV)');
% saveas(gcf,'calculation demo.fig','fig');
%%
close all;
cd(savepathRefine);
clear;
load('summary.mat');
celltype = {'vglut2+','vglut2-','vgat+','vgat-'};
%summary figure
cellNum=cellfun(@(x) length(x),Amp_mean_sig);
cellNumRes=cellfun(@(x) size(x,2)-sum(x==0,2),Amp_mean_sig,'UniformOutput',false);
pRes=cellfun(@(x) x/10,cellNumRes,'UniformOutput',false);
pResCor=pRes;
for i=1:length(pRes)
    pResCor{i}(pRes{i}<0.2,1)=nan;
end
figPRes=bar_error(pResCor,[1,3],celltype,'pRes');

% cellNumAmpIncrease=cellfun(@(x) sum(x>0),slope);
%cellNumAmpIncrease=cellfun(@(x) sum(x>1),PPR);

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
figlatency1=figure;
group=[repmat({celltype{1,1}},1,length(latency1{1,1})),repmat({celltype{1,3}},1,length(latency1{1,3}))];
boxplot([latency1{1,1},latency1{1,3}],group);%plot first pulse latency
title('latency for first response');
saveas(figlatency1,'latency of first response.fig','fig');
figlatency2=bar_error(latency1,[1,3],celltype,'latency for first response(ms)');

%summary of amplitude ratio
AmpR=cell(1,length(celltype));
AmpRMean=cell(1,length(celltype));
for i=1:length(celltype)
    AmpR{1,i}=cellfun(@(x) sum(x(:,6:10),2)./sum(x(:,1:5),2),Amp{1,i},'UniformOutput',false);
    AmpRMean{1,i}=cellfun(@mean,AmpR{1,i});
end
figAmpR=bar_error(AmpRMean,[1,3],celltype,'Amp ratio');
MeanAmpR=cellfun(@(x) sum(x(:,6:10),2)./sum(x(:,1:5),2),Amp_mean,'UniformOutput',false);
figMeanAmpR=bar_error(MeanAmpR,[1,3],celltype,'Mean Amp ratio');
%summary of PPR
figPPR=bar_error(PPR,[1,3],celltype,'PPR');
figure(figPPR);
xlim=get(gca,'Xlim');
ylim=get(gca,'Ylim');
plot([xlim(1),xlim(2)],[1,1],'k--');
hold on;
% plot([1,2],[ylim(2)*0.85,ylim(2)*0.85],'k-','LineWidth',2);
saveas(figPPR,'PPR.png','png');
%summary of slope
figslope=bar_error(slope,[1,3],celltype,'slope');
%summary of first response amplitude
Amp_mean_first=cellfun(@(x) x(:,1),Amp_mean,'UniformOutput',false);
figFirstAmp=bar_error(Amp_mean_first,[1,3],celltype,'first Amplitude(mV)');
%summary of average Amp
Amp_mean_average=cellfun(@(x) mean(x,2),Amp_mean,'UniformOutput',false);
figAmpAverage=bar_error(Amp_mean_average,[1,3],celltype,'Average Amplitude(mV)');
%%%% compare consistency of two method(PPR v.s. slope)
% plasticitySummary=zeros(length(celltype),4);
% x=zeros(1,4);
% for i=1:4
%     plasticitySummary(:,i)=cellfun(@(x) sum(x==i), plasticity);
%     x(i)=sum(plasticitySummary(:,i));
% end
% labels = {'Increase','Decrease','Contradict','No response'};
% cm=colorMean{1};
% for i=2:length(colorMean)
%     cm=[cm;colorMean{i}];
% end
% figsummaryAll=summarypie(x,labels,'All cells',cm);
% for i=1:4
%     figsummarytype=summarypie(plasticitySummary(i,:),labels,celltype{i},cm);
% end

%function part
function figsummary=summarypie(x,labels,titlestr,cm)
for i=1:4
    labels{i}=strcat(labels{i},'(',num2str(x(i)),'/',num2str(sum(x)),')');
end
figsummary=figure;
p=pie(x,labels);
for i=1:length(p)/2
    t=p(i*2);
    t.FontSize=14;
end
colormap(cm);
title(['Summary of ',titlestr,' plasticity calculated by 2 methods']);
saveas(figsummary,['plasticity summary of ',titlestr,'.fig'],'fig');
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
