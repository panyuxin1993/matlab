function [ Amp_mean,slope,PPR,plasticity,latencyOnset,latencyPeak,Amp ] = fplotPlasticityCase( abffile )
%FPLOTPLASTICITYCASE Summary of this function using single abf file to plot
%response to a train of stimuli, and return some values.
%   Input-abffile
%   Output-PPR,slop,1st latency, amplitude, etc.
times=3;
colorMean={[1,0,0],[0,1,0],[0,0,1],[0,0,0]};%{1}-{4} increase, decrease, maintain/contradict, no response
[d,si,h]=abfload(abffile,'channels','a');
ts=(1:h.sweepLengthInPts)*h.si/10^6;
sr=10^6/si;
xlim=[1.5,3.5];
nchannel=size(d,2);
ntrial=size(d,3);
md=mean(d,3);
Amp_mean_sig=zeros(1,10);
Amp_mean=zeros(1,10);%2d matrix, 1d-cell,2d 10 pulse elicit response
PPR=zeros(1,1);
slope=zeros(1,1);
plasticity=zeros(1,1);%1-increase,2-decrease,3-contradict,4-no response
% SNR=[]; %row vector storing each trial signal-to-noise ratio(response amplitude/ std_baseline)
latencyOnset=zeros(1,10);%latency of response onset
latencyPeak=zeros(1,10);%latency of response peak
Amp=zeros(ntrial,10);%each row a trial, each col a pulse, restore all trial all pulse responses
laser_duration=find(md(:,2)>3);%second channel to find when laser on
laser_start=laser_duration(1):0.1*sr:(laser_duration(1)+sr*0.9);%first 10 pulses
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
if sum(Amp_mean_sig(1,:)==0)>=7 % not respond if more than 8 (in 10) stimuli do not respond
    slope(1,1)=nan;
    PPR(1,1)=nan;
    Amp_mean(1,:)=nan;
    latencyOnset(1,:)=nan;
    latencyPeak(1,:)=nan;
else %for responsive cell, re-calculate responses
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
    p=polyfit(1:length(laser_start),Amp_mean(1,:),1);%linear fitting,y=p(1)x+p(2)
    slope(1,1)=p(1);
    PPR(1,1)=Amp_mean(1,2)/Amp_mean(1,1);
    if slope(1,1)*PPR(1,1)==0 ||isnan(slope(1,1)*PPR(1,1))
        warning(strcat('Amp=',num2str(Amp_mean),',PPR=',num2str(PPR),',slope=',num2str(slope)));
    end
end
% temp=latencyOnset(1,:);
% temp(isnan(temp)==1)=[];
% if ~isempty(temp)
%     if temp(1)>10%refined by latency
%         slope(1,1)=nan;
%         PPR(1,1)=nan;
%         Amp_mean(1,:)=nan;
%         latencyOnset(1,:)=nan;
%         latencyPeak(1,:)=nan;
%     end
% end

if slope(1,1)>0 && (PPR(1,1)>1 || isnan(PPR(1,1)))%2 methods same results; or when PPR is nan, only decided by slope
    plot(ts,md(:,1),'-','Color',colorMean{1},'LineWidth',2);
    plasticity(1,1)=1;
elseif slope(1,1)<0 && (PPR(1,1)<1 || isnan(PPR(1,1)))
    plot(ts,md(:,1),'-','Color',colorMean{2},'LineWidth',2);
    plasticity(1,1)=2;
elseif slope(1,1)*(PPR(1,1)-1)<0 %2 methods contradict results
    plot(ts,md(:,1),'-','Color',colorMean{3},'LineWidth',2);
    plasticity(1,1)=3;
else% no response
    plot(ts,md(:,1),'-','Color',colorMean{4},'LineWidth',2);
    plasticity(1,1)=4;
end
hold on;
if 1==1
%     plot([ts(1),ts(end)],[baseline,baseline],'-k','LineWidth',0.1);
%     plot([ts(1),ts(end)],[baseline+STD_baseline*times,baseline+STD_baseline*times],'--k','LineWidth',0.1);
%     plot([ts(1),ts(end)],[baseline-STD_baseline*times,baseline-STD_baseline*times],'--k','LineWidth',0.1);
    peakT=(latencyPeak(1,:)*sr/1000+laser_start)/sr;
    peakT(isnan(peakT)==1)=[];
    peaks=md(round(peakT*sr),1);
    scatter(peakT,peaks,10,'k');
    onsetT=(latencyOnset(1,:)*sr/1000+laser_start)/sr;
    onsetT(isnan(onsetT)==1)=[];
    onsetAmp=md(round(onsetT*sr),1);
    scatter(onsetT,onsetAmp,10,'k','filled');%magenta filled dots as onsets
end
strsummary=['PPR=',num2str(PPR(1,1)),';slope=',num2str(slope(1,1))];
temp=strsplit(abffile,'\');
strfile=temp{end};
set(gca,'xlim',xlim);
xlabel('time(s)');
set(gca,'FontName','Arial','FontSize',14);
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
x_resp=laser_start(Amp_mean_sig(1,:)>0)/sr;
ind_resp=find(Amp_mean_sig(1,:)>0);
if ~isempty(x_resp)
    plot(x_resp,(ylim(end))*ones(length(x_resp),1),'*k');
%         for iPulse=1:length(x_resp)
%             text(x_resp(iPulse),ylim(end)+0.5,[num2str(round(Amp_mean_sig(1,ind_resp(iPulse)),2,'significant')),',',num2str(round(latencyPeak(1,ind_resp(iPulse)),2,'significant'))]);
%         end
end
% text(xlim(1),ylim(end)+1,strsummary);
% set(gca,'ylim',[ylim(1),ylim(2)+1]);
% clear ylim;
text(xlim(1),ylim(end),strfile);
box off;
end

