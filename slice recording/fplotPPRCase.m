function [ Amp_mean, currentPPR,latencyOnset,latencyPeak,Amp] = fplotPPRCase( abffile )
%FPLOTPPRCASE use single abf file to plot mean trace of paired pulse
%stimuli and some values
%   input is abf file path
%   output include
colorMean={[1,0,0],[0,1,0],[0,0,0]};
times=3;
[d,si,h]=abfload(abffile,'channels','a');
ts=(1:h.sweepLengthInPts)*h.si/10^6;
sr=10^6/si;
xlim=[1.5,2.5];
nchannel=size(d,2);
ntrial=size(d,3);
md=nanmean(d,3);
Amp_mean_sig=zeros(1,2);
Amp_mean=zeros(1,2);%2d matrix, 1d-cell,2d 2 pulse elicit response
% PPR=zeros(1,1);
latencyOnset=zeros(1,2);%latency of response onset
latencyPeak=zeros(1,2);%latency of response peak
% Amp=zeros(ntrial,2);%each row a trial, each col a pulse, restore all trial all pulse responses
laser_duration=find(md(:,2)>3);%second channel to find when laser on
laser_start=[laser_duration(1) laser_duration(ceil(length(laser_duration)/2)+1)];%2 pulses, so 2 start; for odd number,may use the 2nd time point at 2nd stimulus
%         stim_interval=laser_start(2)-laser_start(1);
% tempSI=laser_duration(2:end)-laser_duration(1:end-1);
% stim_interval=max(tempSI);%time from end of 1st stimulus to start of 2nd stimulus
STD_baseline = std(md(laser_start(1)-1*sr:laser_start(1),1));%1s before laser start as baseline
baseline=mean(md(laser_start(1)-1*sr:laser_start(1),1));
for npulse=1:length(laser_start)
    peak=max(md(laser_start(npulse):laser_start(npulse)+0.05*sr,1)); %0.05s within after stimulus
    amp=abs(peak-baseline);
    if amp<=STD_baseline*times
        Amp_mean_sig(1,npulse)=0;
    elseif amp>STD_baseline*times
        Amp_mean_sig(1,npulse)=amp;
    end
end

if sum(Amp_mean_sig(1,1:2))==0 % not respond if both responses are zero
    currentPPR=nan;
    Amp_mean(1,:)=nan;
    latencyOnset(1,:)=nan;
    latencyPeak(1,:)=nan;
    Amp=zeros(ntrial,2);%each row a trial, each col a pulse
else %for responsive cell, re-calculate responses
    Amp=zeros(ntrial,2);%each row a trial, each col a pulse
    for itrial=1:ntrial
        for npulse=1:length(laser_start)
            ibaseline=mean(d(laser_start(1)-1*sr:laser_start(1),1,itrial));%1s before laser start as baseline
            ipeak=max(d(laser_start(npulse):laser_start(npulse)+0.05*sr,1,itrial)); %0.05s within after stimulus
            Amp(itrial,npulse)=abs(ipeak-ibaseline);
        end
    end
    for npulse=1:length(laser_start)
        peak=max(md(laser_start(npulse):laser_start(npulse)+0.05*sr,1)); %0.05s within after stimulus
        peakT=find(md(laser_start(npulse):laser_start(npulse)+0.05*sr,1)==peak);
        onsetT=find(abs(md(laser_start(npulse):laser_start(npulse)+0.05*sr,1)-baseline)>STD_baseline*times);
        Amp_mean(1,npulse)=abs(peak-baseline);
        latencyPeak(1,npulse)=peakT(1)/sr*1000;% ms from stim_onset
        if isempty(onsetT)
            latencyOnset(1,npulse)=nan;
        else
            latencyOnset(1,npulse)=onsetT(1)/sr*1000;
        end
    end
    currentPPR=Amp_mean(1,2)/Amp_mean(1,1);
end
%         %refined by latency
%         temp=latencyOnset(1,:);
%         temp(isnan(temp)==1)=[];
%         if ~isempty(temp)
%             if temp(1)>10%refined by latency
%                 PPR{1,n_celltype}(n_abf,1)=nan;
%                 Amp_mean(1,:)=nan;
%                 latencyOnset(1,:)=nan;
%                 latencyPeak(1,:)=nan;
%             else
%                 latency1{1,n_celltype}=[latency1{1,n_celltype} temp(1)];
%             end
%         end

%plot
if  currentPPR>1 %increase
    plot(ts,md(:,1),'-','Color',colorMean{1},'LineWidth',2);
elseif currentPPR<1 %decrease
    plot(ts,md(:,1),'-','Color',colorMean{2},'LineWidth',2);
else% no response, here currentPPR=nan
    plot(ts,md(:,1),'-','Color',colorMean{3},'LineWidth',2);
end
hold on;
% plot([ts(1),ts(end)],[baseline,baseline],'-k','LineWidth',0.1);
% plot([ts(1),ts(end)],[baseline+STD_baseline*times,baseline+STD_baseline*times],'--k','LineWidth',0.1);
% plot([ts(1),ts(end)],[baseline+STD_baseline,baseline+STD_baseline],'-.k','LineWidth',0.1);
% plot([ts(1),ts(end)],[baseline-STD_baseline*times,baseline-STD_baseline*times],'--k','LineWidth',0.1);
% plot([ts(1),ts(end)],[baseline-STD_baseline,baseline-STD_baseline],'-.k','LineWidth',0.1);
peakT=(latencyPeak(1,:)*sr/1000+laser_start)/sr;
peakT(isnan(peakT)==1)=[];
peaks=md(round(peakT*sr),1);
scatter(peakT,peaks,10,'k');
onsetT=(latencyOnset(1,:)*sr/1000+laser_start)/sr;
onsetT(isnan(onsetT)==1)=[];
onsetAmp=md(round(onsetT*sr),1);
scatter(onsetT,onsetAmp,10,'k','filled');%magenta filled dots as onsets
strsummary=['PPR=',num2str(currentPPR)];
set(gca,'xlim',xlim);
xlabel('time(s)');
set(gca,'FontName','Arial','FontSize',14);
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
x_resp=laser_start(Amp_mean_sig(1,:)>0)/sr;
ind_resp=find(Amp_mean_sig(1,:)>0);
if ~isempty(x_resp)
    plot(x_resp,(ylim(end))*ones(length(x_resp),1),'*k');
%     for iPulse=1:length(x_resp)
%         text(x_resp(iPulse),ylim(end)+0.5,[num2str(round(Amp_mean_sig(1,ind_resp(iPulse)),2,'significant')),',',num2str(round(latencyPeak(1,ind_resp(iPulse)),2,'significant'))]);
%     end
end
text(xlim(1),ylim(end),strsummary);
box off;

end

