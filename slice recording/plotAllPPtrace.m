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
savepathPPR=['D:\xulab\project\slice recording\data3\PPR summary\PPR cases\'];

for n_celltype=1:length(celltype)
    AbfFiles = dir(strcat(rootpath,celltype{n_celltype},'\*.abf'));
    cd(rootpath);    
    ncell=length(AbfFiles);
    Amp{1,n_celltype}=cell(1,length(AbfFiles));%each cell a matrix containing individual amp
    Amp_mean{1,n_celltype}=zeros(length(AbfFiles),3);%2d matrix, 1d-cell,2d 2 pulse elicit response and interval
    Amp_mean_sig{1,n_celltype}=zeros(length(AbfFiles),3);
    PPR{1,n_celltype}=zeros(length(AbfFiles)/3,3);%1d different cells,2d three intervals-50ms,100ms,200ms
    plasticity{1,n_celltype}=zeros(length(AbfFiles),1);%1-increase,2-decrease,0-no response
    SNR{1,n_celltype}=[]; %row vector storing each trial signal-to-noise ratio(response amplitude/ std_baseline)
    latency{1,n_celltype}=zeros(length(AbfFiles),2);%latency of response onset
    latency{2,n_celltype}=zeros(length(AbfFiles),2);%latency of response peak
    latency1{1,n_celltype}=[];
    figraw=figure;
    ncol=6;
    nrow=ceil(ncell/ncol);
    set(gcf,'position',[0,0,300*ncol,300*nrow]);
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
        
            figure(figraw);
            subplot(nrow,ncol,n_abf);
            if  currentPPR>1 %increase
                plot(ts,md(:,1),'-','Color',colorMean{1},'LineWidth',2); 
                plasticity{1,n_celltype}(n_abf,1)=1;
            elseif currentPPR<1 %decrease
                plot(ts,md(:,1),'-','Color',colorMean{2},'LineWidth',2);
                plasticity{1,n_celltype}(n_abf,1)=2;
            else% no response, here currentPPR=nan
                plot(ts,md(:,1),'-','Color',colorMean{3},'LineWidth',2);
                plasticity{1,n_celltype}(n_abf,1)=0;
            end
            hold on;
            
%             plot([ts(1),ts(end)],[baseline,baseline],'-k','LineWidth',0.1);
%             plot([ts(1),ts(end)],[baseline+STD_baseline*times,baseline+STD_baseline*times],'--k','LineWidth',0.1);
%             plot([ts(1),ts(end)],[baseline+STD_baseline,baseline+STD_baseline],'-.k','LineWidth',0.1);
%             plot([ts(1),ts(end)],[baseline-STD_baseline*times,baseline-STD_baseline*times],'--k','LineWidth',0.1);
%             plot([ts(1),ts(end)],[baseline-STD_baseline,baseline-STD_baseline],'-.k','LineWidth',0.1);
            peakT=(latency{2,n_celltype}(n_abf,:)*sr/1000+laser_start)/sr;
            peakT(isnan(peakT)==1)=[];
            peaks=md(round(peakT*sr),1);
%             scatter(peakT,peaks,10,'k');
            onsetT=(latency{1,n_celltype}(n_abf,:)*sr/1000+laser_start)/sr;
            onsetT(isnan(onsetT)==1)=[];
            onsetAmp=md(round(onsetT*sr),1);
%             scatter(onsetT,onsetAmp,10,'k','filled');%magenta filled dots as onsets
            
            strsummary=['PPR=',num2str(currentPPR)];
            set(gca,'xlim',xlim);
            if ceil(n_abf/ncol)==nrow
                xlabel('time(s)');
            end
%             set(gca,'FontName','Arial','FontSize',14);
            xlim=get(gca,'xlim');
            ylim=get(gca,'ylim');
            x_resp=laser_start(Amp_mean_sig{1,n_celltype}(n_abf,1:end-1)>0)/sr;
            ind_resp=find(Amp_mean_sig{1,n_celltype}(n_abf,1:end-1)>0);
            if ~isempty(x_resp)
                plot(x_resp,(ylim(end))*ones(length(x_resp),1),'*k');
%                 for iPulse=1:length(x_resp)
%                     text(x_resp(iPulse),ylim(end)+0.5,[num2str(round(Amp_mean_sig{1,n_celltype}(n_abf,ind_resp(iPulse)),2,'significant')),',',num2str(round(latency{2,n_celltype}(n_abf,ind_resp(iPulse)),2,'significant'))]);
%                 end
            end
%             text(xlim(1),ylim(end),strsummary);
            box off;  
    end
    suptitle(celltype{n_celltype});
    saveas(figraw,[savepathPPR,celltype{n_celltype},'-all cases-no label.fig'],'fig');
    saveas(figraw,[savepathPPR,celltype{n_celltype},'-all cases-no label.png'],'png');
%     saveas(figraw,[savepathPPR,celltype{n_celltype},'-all cases-no label.eps'],'epsc');
end
save([savepathPPR,'summary.mat'],'Amp','Amp_mean','PPR','plasticity','Amp_mean_sig','latency','latency1');