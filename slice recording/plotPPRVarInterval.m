%first decide which cell respond to stimuli using amplitude>3 std baseline,
%then no response cells' amplitude is nan and responsive cells' amplitude
%is |peak-baseline| no matter whether its SNR>3 
rootpath='D:\xulab\project\slice recording\data3\PPR summary\';
celltype = {'vglut2+','vglut2-','vgat+','vgat-'};
Amp=cell(1,length(celltype)); %response amp for each pulse of each trials
Amp_mean=cell(1,length(celltype));
Amp_mean_sig=cell(1,length(celltype));
SNR=cell(1,length(celltype));
colorCases={[1,0.7,0.7],[0.7,1,0.7],[0.7,0.7,0.7]};%colorCase{1}-{3} increase, decrease, no response
colorMean={[1,0,0],[0,1,0],[0,0,0]};
PPR=cell(1,length(celltype));%4 cell types
plasticity=cell(1,length(celltype));%summary results from PPR and slope
latency=cell(2,length(celltype));%summary of response latency from stim_onset; the 1st row is latency of reponse onset, the  2nd row is latency of response peak
latency1=cell(1,length(celltype));%1st pulse latency
times=3;
savepathPPR=['D:\xulab\project\slice recording\data3\PPR summary\PPR\'];

%{
for n_celltype=1:length(celltype)
    AbfFiles = dir(strcat(rootpath,celltype{n_celltype},'\*.abf'));
    cd(rootpath);    
    Amp{1,n_celltype}=cell(1,length(AbfFiles));%each cell a matrix containing individual amp
    Amp_mean{1,n_celltype}=zeros(length(AbfFiles),3);%2d matrix, 1d-cell,2d 2 pulse elicit response and interval
    Amp_mean_sig{1,n_celltype}=zeros(length(AbfFiles),3);
    PPR{1,n_celltype}=zeros(length(AbfFiles)/3,3);%1d different cells,2d three intervals-50ms,100ms,200ms
    plasticity{1,n_celltype}=zeros(length(AbfFiles),1);%1-increase,2-decrease,0-no response
    SNR{1,n_celltype}=[]; %row vector storing each trial signal-to-noise ratio(response amplitude/ std_baseline)
    latency{1,n_celltype}=zeros(length(AbfFiles),2);%latency of response onset
    latency{2,n_celltype}=zeros(length(AbfFiles),2);%latency of response peak
    latency1{1,n_celltype}=[];
    for n_abf=1:length(AbfFiles)
%             strcat(cellFolder,AbfFiles(n_abf).name);
        [d,si,h]=abfload(strcat(rootpath,celltype{n_celltype},'\',AbfFiles(n_abf).name),'channels','a');
        ts=(1:h.sweepLengthInPts)*h.si/10^6;
        sr=10^6/si;
        xlim=[1.5,2.5];
        nchannel=size(d,2);
        ntrial=size(d,3);
        dRefine=d;
        dRefine(dRefine(:,1,:)>0)=nan;%some traces fire so set nan and not mean them
        md=nanmean(d,3);
        laser_duration=find(md(:,2)>3);%second channel to find when laser on
        laser_start=[laser_duration(1) laser_duration(ceil(length(laser_duration)/2)+1)];%2 pulses, so 2 start; for odd number,may use the 2nd time point at 2nd stimulus
%         stim_interval=laser_start(2)-laser_start(1);
        tempSI=laser_duration(2:end)-laser_duration(1:end-1);
        stim_interval=max(tempSI);%time from end of 1st stimulus to start of 2nd stimulus
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
        Amp_mean_sig{1,n_celltype}(n_abf,3) =stim_interval*h.si/10^6*1000;%should be either 50,100 or 200(ms)
        if sum(Amp_mean_sig{1,n_celltype}(n_abf,1:2))==0 % not respond if both responses are zero
            if Amp_mean_sig{1,n_celltype}(n_abf,3)<60
                PPR{1,n_celltype}(ceil(n_abf/3),1)=nan;
            elseif Amp_mean_sig{1,n_celltype}(n_abf,3)>60 && Amp_mean_sig{1,n_celltype}(n_abf,3)<110
                PPR{1,n_celltype}(ceil(n_abf/3),2)=nan;
            else
                PPR{1,n_celltype}(ceil(n_abf/3),3)=nan;
            end
            currentPPR=nan;
            Amp_mean{1,n_celltype}(ceil(n_abf/3),1:2)=nan;
            Amp_mean{1,n_celltype}(ceil(n_abf/3),3)=Amp_mean_sig{1,n_celltype}(n_abf,3);
            latency{1,n_celltype}(n_abf,:)=nan;
            latency{2,n_celltype}(n_abf,:)=nan;
            Amp{1,n_celltype}{1,n_abf}=zeros(1,2);
        else %for responsive cell, re-calculate responses
            Amp{1,n_celltype}{1,n_abf}=zeros(ntrial,2);%each row a trial, each col a pulse
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
            currentPPR=Amp_mean{1,n_celltype}(n_abf,2)/Amp_mean{1,n_celltype}(n_abf,1);
            if Amp_mean_sig{1,n_celltype}(n_abf,3)<60
                PPR{1,n_celltype}(ceil(n_abf/3),1)=currentPPR;
            elseif Amp_mean_sig{1,n_celltype}(n_abf,3)>60 && Amp_mean_sig{1,n_celltype}(n_abf,3)<110
                PPR{1,n_celltype}(ceil(n_abf/3),2)=currentPPR;
            else
                PPR{1,n_celltype}(ceil(n_abf/3),3)=currentPPR;
            end  
        end       
%         %refined by latency
%         temp=latency{1,n_celltype}(n_abf,:);
%         temp(isnan(temp)==1)=[];      
%         if ~isempty(temp)
%             if temp(1)>10%refined by latency
%                 PPR{1,n_celltype}(n_abf,1)=nan;
%                 Amp_mean{1,n_celltype}(n_abf,:)=nan;
%                 latency{1,n_celltype}(n_abf,:)=nan;
%                 latency{2,n_celltype}(n_abf,:)=nan;
%             else
%                 latency1{1,n_celltype}=[latency1{1,n_celltype} temp(1)];
%             end     
%         end
        
        figraw=figure;
        for nsubplot=1:nchannel
            figure(figraw);
            subplot(2,1,nsubplot);
            if  currentPPR>1 %increase
                plot(ts,md(:,nsubplot),'-','Color',colorMean{1},'LineWidth',2); 
                plasticity{1,n_celltype}(n_abf,1)=1;
            elseif currentPPR<1 %decrease
                plot(ts,md(:,nsubplot),'-','Color',colorMean{2},'LineWidth',2);
                plasticity{1,n_celltype}(n_abf,1)=2;
            else% no response, here currentPPR=nan
                plot(ts,md(:,nsubplot),'-','Color',colorMean{3},'LineWidth',2);
                plasticity{1,n_celltype}(n_abf,1)=0;
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
            strsummary=['PPR=',num2str(currentPPR)];
            set(gca,'xlim',xlim);
            xlabel('time(s)');
            set(gca,'FontName','Arial','FontSize',14);
            xlim=get(gca,'xlim');
            ylim=get(gca,'ylim');
            x_resp=laser_start(Amp_mean_sig{1,n_celltype}(n_abf,1:end-1)>0)/sr;
            ind_resp=find(Amp_mean_sig{1,n_celltype}(n_abf,1:end-1)>0);
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
        saveas(figraw,[savepathPPR,celltype{n_celltype},'-',AbfFiles(n_abf).name,'.fig'],'fig');
        close all;
    end
end
save([savepathPPR,'summary.mat'],'Amp','Amp_mean','PPR','plasticity','Amp_mean_sig','latency','latency1');
%}
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
cd(savepathPPR);
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

% figlatency1=figure;
% group=[repmat({celltype{1,1}},1,length(latency1{1,1})),repmat({celltype{1,3}},1,length(latency1{1,3}))];
% boxplot([latency1{1,1},latency1{1,3}],group);%plot first pulse latency
% title('latency for first response');
% saveas(figlatency1,'latency of first response.fig','fig');
% figlatency2=bar_error(latency1,[1,3],celltype,'latency for first response(ms)');

%summary of PPR
% PPR50=cellfun(@(x) x(:,1),PPR,'UniformOutput',false);
% PPR100=cellfun(@(x) x(:,2),PPR,'UniformOutput',false);
% PPR200=cellfun(@(x) x(:,3),PPR,'UniformOutput',false);
figPPR=bar_error(PPR,[1,3],celltype,'PPR');
figure(figPPR);
xlim=get(gca,'Xlim');
ylim=get(gca,'Ylim');
plot([xlim(1),xlim(2)],[1,1],'k--');
hold on;
% plot([1,2],[ylim(2)*0.85,ylim(2)*0.85],'k-','LineWidth',2);
saveas(figPPR,'PPR.png','png');

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

function [str]=plabelsymbol(pvalue)
if pvalue<0.05 && pvalue>=0.01
    str=' *';
elseif pvalue<0.01 && pvalue>=0.001
    str=' **';
elseif pvalue<0.001
    str=' ***';
else
    pvalue=round(pvalue,2);%������λ����
    str=strcat('p=',num2str(pvalue));
end
end