function [ latencyPeak, latencyOnset, amplitude,amplitude_meantrace, amplitude_sig,amplitude_sig_meantrace,amplitude_sig_average,pro_responsive,pro_responsive_sig ] = fplotEPSPCase( abffile )
%FPLOTEPSPCASE Summary of this function using single abf file to plot a
%mean trace of EPSP(only mean trace) and return some values
%   Input is abf file path
%   Output is latency, amplitude_sig,..
colorCases={[0.7,0.7,0.7],[1,0.9,0.7]};%colorCase{2} highlighted
colorMean={[0,0,0],[1,0.5,0],[1,0,0]};
times=3;
[d,si,h]=abfload(abffile,'channels','a');
ts=(1:h.sweepLengthInPts)*h.si/10^6;
sr=10^6/si;
xlim=[1.5,2.5];
md=mean(d,3);
amplitude_sig=zeros(size(d,3)+1,1);% vector of each individual amplitude_sigs,last one is mean trace peak
amplitude=zeros(size(d,3)+1,1);%max as peak even if not exceed 3std
sumStimTrial=0;
for i=1:size(d,3)
    laser_duration=find(d(:,2,i)>1);%second channel to find when laser on
    if isempty(laser_duration) %if some trials no stimuli, then pass this trial
        continue;
    else
        sumStimTrial=sumStimTrial+1;
    end
    laser_start=laser_duration(1);
    STD_baseline = std(md(laser_start(1)-0.5*sr:laser_start(1),1));%0.5s before laser start as baseline
    baseline=mean(md(laser_start(1)-1*sr:laser_start(1),1));%1s before stimuli as baseline
    peak=max(d(laser_start:laser_start+0.05*sr,1,i)); %50ms within after stimulus
    amp=abs(peak-baseline);
    amplitude(i)=amp;
    if amp>STD_baseline*times %else remain 0
        amplitude_sig(i)=amp;
        hold on;
    end
end
%for mean trace
laser_duration=find(md(:,2)>3);%second channel to find when laser on
laser_start=laser_duration(1);
STD_baseline = std(md(laser_start(1)-1*sr:laser_start(1),1));%1s before laser start as baseline
baseline=mean(md(laser_start(1)-1*sr:laser_start(1),1));
peak=max(md(laser_start:laser_start+0.05*sr,1)); %0.1s within after stimulus
peakT=find(md(laser_start:laser_start+0.05*sr,1)==peak);
onsetT=find(abs(md(laser_start:laser_start+0.05*sr,1)-baseline)>STD_baseline*times);
latencyPeak=peakT(1)/sr*1000;% ms from stim_onset
if isempty(onsetT)
    latencyOnset=nan;
else
    latencyOnset=onsetT(1)/sr*1000;
end
if peak<=-50
    amp=abs(peak-baseline);
    isspike=0;
else
    peak=-50;
    amp=abs(peak-baseline);
    isspike=1;
end
amplitude(end)=amp;
sum_responsive=sum(double(amplitude_sig(1:end-1)>0)); %last one is mean trace peak
pro_responsive=sum_responsive/sumStimTrial;%rather than assuming that all trials have activity
pro_responsive_sig=pro_responsive;
%plot mean trace
if isspike==1%this is a spike
    amplitude_sig(end)=amp;
    plot(ts,md(:,1),'-','Color',colorMean{3},'LineWidth',2);
    hold on;
elseif amp<=STD_baseline*times %|| latencyOnset>10% latency>10 data be got rid off
    amplitude_sig(end)=nan;
    latencyOnset=nan;
    latencyPeak=nan;
    pro_responsive_sig=nan;
    plot(ts,md(:,1),'-','Color',colorMean{1},'LineWidth',2);
    hold on;
elseif amp>STD_baseline*times
    amplitude_sig(end)=amp;
    plot(ts,md(:,1),'-','Color',colorMean{2},'LineWidth',2);
    hold on;
end
peakT=(latencyPeak*sr/1000+laser_start)/sr;
peakT(isnan(peakT)==1)=[];
peaks=md(round(peakT*sr),1);
scatter(peakT,peaks,10,'k');
onsetT=(latencyOnset*sr/1000+laser_start)/sr;
onsetT(isnan(onsetT)==1)=[];
onsetAmp=md(round(onsetT*sr),1);
scatter(onsetT,onsetAmp,10,'k','filled');%magenta filled dots as onsets
amplitude_sig_meantrace=amplitude_sig(end);
amplitude_meantrace=amplitude(end);
if sum(amplitude_sig(1:end-1)>0)==0
    amplitude_sig_average=nan;
else
    amplitude_sig_average=sum(amplitude_sig(amplitude_sig(1:end-1)>0))/sum(amplitude_sig(1:end-1)>0); %exclude no response trials and only average among those responsive
end

strsummary=['ampMean=',num2str(amplitude_sig_meantrace),'mV,resOnsetLatency',num2str(latencyOnset),'ms'];
set(gca,'xlim',xlim);
xlabel('time(s)');
ylabel('V(mV)');
set(gca,'FontName','Arial','FontSize',14);
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');

text(xlim(1),ylim(end)+1,strsummary);
set(gca,'ylim',[ylim(1),ylim(2)+1]);
clear ylim;
box off;
end


