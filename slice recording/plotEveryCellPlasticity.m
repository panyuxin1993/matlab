%this script is for plotting example into one figure, and save as eps 

%first decide which cell respond to stimuli using amplitude>3 std baseline,
%then no response cells' amplitude is nan and responsive cells' amplitude
%is |peak-baseline| no matter whether its SNR>3 
%rootpath='D:\xulab\project\slice recording\XUNL\summary of short term plasticity\';%for plasticity
rootpath='D:\xulab\project\slice recording\data3\response summary\';%for EPSP
celltype = {'vglut2+','vglut2-','vgat+','vgat-'};
Amp_mean=cell(1,length(celltype));
Amp_mean_sig=cell(1,length(celltype));
SNR=cell(1,length(celltype));
colorCases={[1,0.7,0.7],[0.7,1,0.7],[0.7,0.7,1],[0.7,0.7,0.7]};%colorCase{1}-{4} increase, decrease, maintain/contradict, no response
colorMean={[1,0,0],[0,1,0],[0,0,1],[0,0,0]};
PPR=cell(1,length(celltype));%only compare the first 2 response, ratio >1 or <1
slope=cell(1,length(celltype));%linearly fit the whole trend of response change, slope positive or negative
plasticity=cell(1,length(celltype));%summary results from PPR and slope
latency=cell(2,length(celltype));%summary of response latency from stim_onset; the 1st row is latency of reponse onset, the  2nd row is latency of response peak
times=3;
savepath='D:\xulab\project\slice recording\data3\response summary\response cases\';%for EPSP
%rootpath='D:\xulab\project\slice recording\XUNL\summary of short term plasticity\example\';%for plasticity
for n_celltype=1:length(celltype)
    AbfFiles = dir(strcat(rootpath,celltype{n_celltype},'\*.abf'));
    cd(rootpath);    
    cellNum=length(AbfFiles);%currently, <=22
    ncol=5;
    nrow=ceil(cellNum/ncol);
    Amp_mean{1,n_celltype}=zeros(length(AbfFiles),10);%2d matrix, 1d-cell,2d 10 pulse elicit response
    Amp_mean_sig{1,n_celltype}=zeros(length(AbfFiles),10);
    PPR{1,n_celltype}=zeros(length(AbfFiles),1);
    slope{1,n_celltype}=zeros(length(AbfFiles),1);
    plasticity{1,n_celltype}=zeros(length(AbfFiles),1);%1-increase,2-decrease,3-contradict,4-no response
    SNR{1,n_celltype}=[]; %row vector storing each trial signal-to-noise ratio(response amplitude/ std_baseline)
    latency{1,n_celltype}=zeros(length(AbfFiles),10);%latency of response onset
    latency{2,n_celltype}=zeros(length(AbfFiles),10);%latency of response peak
    figexample=figure;
    set(gcf,'Position',[0,0,400*ncol,300*nrow]);
    for n_abf=1:length(AbfFiles)
%             strcat(cellFolder,AbfFiles(n_abf).name);
        [d,si,h]=abfload(strcat(rootpath,celltype{n_celltype},'\',AbfFiles(n_abf).name),'channels','a');
        ts=(1:h.sweepLengthInPts)*h.si/10^6;
        sr=10^6/si;
        xlim=[1.5,3.5];
        nchannel=size(d,2);
        md=mean(d,3);
        laser_duration=find(md(:,2)>3);%second channel to find when laser on
        laser_start=laser_duration(1):0.1*sr:(laser_duration(1)+sr*0.9);%first 10 pulses
%         STD_baseline = std(md(laser_start(1)-1*sr:laser_start(1),1));%1s before laser start as baseline
%         baseline=mean(md(laser_start(1)-1*sr:laser_start(1),1));
        STD_baseline= std(md(:,1));
        baseline=mean(md(:,1));
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
        else %for responsive cell, re-calculate responses
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
        figure(figexample);
        subplot(nrow,ncol,n_abf);
        if slope{1,n_celltype}(n_abf,1)>0 && (PPR{1,n_celltype}(n_abf,1)>1 || isnan(PPR{1,n_celltype}(n_abf,1)))%2 methods same results; or when PPR is nan, only decided by slope
            plot(ts,md(:,1),'-','Color',colorMean{1},'LineWidth',2);
            plasticity{1,n_celltype}(n_abf,1)=1;
        elseif slope{1,n_celltype}(n_abf,1)<0 && (PPR{1,n_celltype}(n_abf,1)<1 || isnan(PPR{1,n_celltype}(n_abf,1)))
            plot(ts,md(:,1),'-','Color',colorMean{2},'LineWidth',2);
            plasticity{1,n_celltype}(n_abf,1)=2;
        elseif slope{1,n_celltype}(n_abf,1)*(PPR{1,n_celltype}(n_abf,1)-1)<0 %2 methods contradict results
            plot(ts,md(:,1),'-','Color',colorMean{3},'LineWidth',2);
            plasticity{1,n_celltype}(n_abf,1)=3;
        else
            plot(ts,md(:,1),'-','Color',colorMean{4},'LineWidth',2);
            plasticity{1,n_celltype}(n_abf,1)=4;
        end
        hold on;
        strsummary=[AbfFiles(n_abf).name,'-PPR=',num2str(PPR{1,n_celltype}(n_abf,1)),'; slope=',num2str(slope{1,n_celltype}(n_abf,1))];
        set(gca,'xlim',xlim);
        xlabel('time(s)');
        set(gca,'FontName','Arial','FontSize',14);
        xlim=get(gca,'xlim');
        ylim=get(gca,'ylim');
        x_resp=laser_start(Amp_mean_sig{1,n_celltype}(n_abf,:)>0)/sr;
        ind_resp=find(Amp_mean_sig{1,n_celltype}(n_abf,:)>0);
        if ~isempty(x_resp)
            plot(x_resp,(ylim(end))*ones(length(x_resp),1),'*k');
%             for iPulse=1:length(x_resp)
%                 text(x_resp(iPulse),ylim(end)+0.5,num2str(round(Amp_mean_sig{1,n_celltype}(n_abf,ind_resp(iPulse)),2,'significant')));
%             end
        end
        text(xlim(1),ylim(end)+1,strsummary);
        set(gca,'ylim',[ylim(1),ylim(2)+1]);
        clear ylim;
        box off;
    end
    suptitle(celltype{n_celltype});
    saveas(figexample,[savepath,celltype{n_celltype},'.pdf'],'pdf');
    saveas(figexample,[savepath,celltype{n_celltype},'.png'],'png');
end
save([savepath,'summary.mat'],'Amp_mean','PPR','slope','plasticity','Amp_mean_sig','latency');
