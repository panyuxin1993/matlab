function [spike_vec,spike_unit,mid_data,figPCA,figTrace] = fGetSpikes(data,ts,threSTD,crit)
%FGETSPIKES Using votage data recorded using open-ephys to deduce spikes
%Input- 
%   data- n-by-1 vector
%   threSTD- threshod for extracting data for analysis
%   crit- criterial for determine spikes
%Output-
%   spike_vec- n-by-1 vector made of 1/0
fr=30000;%Hz
frameNumTime=[0.0003,0.0007];
frameNum=fr*frameNumTime;
figTrace=figure;
plot(ts,data);
hold on;
data_std=zeros(size(data));
for i=1:3000:length(data)%non-overlapping bin
    ind_end=min(i+3000,length(data));
    data_std(i:ind_end)=std(data(i:ind_end));
end
plot([ts(1),ts(end)],(mean(data)-std(data)*threSTD)*ones(1,2),'k--');
plot([ts(1),ts(end)],(mean(data)+std(data)*threSTD)*ones(1,2),'k--');
plot(ts,mean(data)-data_std*threSTD,'k:');
plot(ts,mean(data)+data_std*threSTD,'k:');
% ind_exceed=data<(mean(data)-std(data)*threSTD);
ind_exceed=data<(mean(data)-data_std*threSTD);
exceed_start=diff(ind_exceed);
exceed_start=[0;exceed_start];
ind_exceed_start=find(exceed_start>0);
ind_exceed_end=find(exceed_start<0);
ind_exceed_end(ind_exceed_end<ind_exceed_start(1))=[];%in case data start from a case where first frame data exceed threshold 
if length(ind_exceed_end) ~= length(ind_exceed_start)
    warning(['start - end = ',num2str(length(ind_exceed_start)-length(ind_exceed_end))]);
end
n_spk=min(length(ind_exceed_end),length(ind_exceed_start));
ind_exceed_start=ind_exceed_start(1:n_spk);
ind_exceed_end=ind_exceed_end(1:n_spk);
ind_peak_relative=zeros(size(ind_exceed_start));%relative index of peak to ind_exceed_start
flag_spk=zeros(size(ind_exceed_start));
spkdata=zeros(n_spk,sum(frameNum)+1);
spkdata_aligned=zeros(n_spk,sum(frameNum)+1);
for i=1:n_spk
    ind1=max(1,ind_exceed_start(i)-frameNum(1));
    ind2=min(ind_exceed_start(i)+frameNum(2),length(data));
    if ind_exceed_start(i)-frameNum(1)<1
        spkdata(i,frameNum(1)-ind_exceed_start(i)+2:end)=data(ind1:ind2)';
    elseif ind_exceed_start(i)+frameNum(2)>length(data)
        spkdata(i,1:ind2-ind1+1)=data(ind1:ind2)';
    else
        spkdata(i,:)=data(ind1:ind2)';
    end
    exceed_data=data(ind_exceed_start(i):ind_exceed_end(i));
    ind_peak_temp=find(min(exceed_data)==exceed_data);
    ind_peak_relative(i)=ind_peak_temp(1)-1;%in case peak last a long time
    ind1=max(1,ind_exceed_start(i)+ind_peak_relative(i)-frameNum(1));
    ind2=min(ind_exceed_start(i)+ind_peak_relative(i)+frameNum(2),length(data));
    %align spkdata to their peaks
    if ind_exceed_start(i)+ind_peak_relative(i)-frameNum(1)<1
        spkdata_aligned(i,frameNum(1)-ind_exceed_start(i)-ind_peak_relative(i)+2:end)=data(ind1:ind2)';
    elseif ind_exceed_start(i)+ind_peak_relative(i)+frameNum(2)>length(data)
        spkdata_aligned(i,1:ind2-ind1+1)=data(ind1:ind2)';
    else
        spkdata_aligned(i,:)=data(ind1:ind2)';
    end
end
ts_aligned=(-frameNum(1):frameNum(2))/fr;
flag_spk(min(spkdata(:,frameNum(1):frameNum(1)+15),[],2)<mean(data)-std(data)*crit)=1;%within 0.5ms, should reach peak if it had
flag_spk=logical(flag_spk);

[coeff,score,latent,tsquared,explained,mu] = pca(spkdata_aligned);
fig_noise=figure;%see noise level across time, where large noise level indicates animal movements. 
histogram(data_std);
hold on;
y_lim=get(gca,'Ylim');
plot(std(data)*ones(1,2),y_lim,'r-');


figPCA=figure;
ax1 = subplot(2,3,1);
yyaxis(ax1,'left');
% plot(explained,'r-');hold on;
plot(cumsum(explained(1:10)),'k-');
ylabel('Explained varience');
title('nPCs/nClusters');
%clustering to separate data, group noise to one group, ref https://blog.csdn.net/qq_44646352/article/details/124266723
idx_n=1:10;
eva = evalclusters(score(:,1:3),'kmeans','CalinskiHarabasz','KList',idx_n);
yyaxis(ax1,'right');
plot(eva);
n_cluster=eva.OptimalK;
[idx] = kmeans(score(:,1:3),n_cluster);

tb_group=tabulate(idx);
idx_noise=tb_group(tb_group(:,1)==max(tb_group(:,1)),1);%n-by-3 matrix; assumption may not be the case, all spikes may be real spikes
n_unit=size(tb_group,1);
color_3d=cell(n_unit,1);
for i=1:n_unit
    color_3d{i}=[1-i/n_unit,0,i/n_unit];
end
color_3d{idx_noise}=[0.7,0.7,0.7];
color_mat=cell2mat(color_3d);

subplot(2,3,2);
hold on;
for i=1:size(tb_group,1)    
    idx_unit=tb_group(i,1);
    ind_unit=ind_exceed_start(idx==idx_unit);
    i_spks=(idx==idx_unit);
    %plot3(score(i_spks,1),score(i_spks,2),score(i_spks,3),'.','Color',color_3d{i});
    %plot waveform together for each unit
    plot(repmat(ts_aligned',1,sum(i_spks)),spkdata_aligned(i_spks,:)','Color',color_3d{i});
end

subplot(2,3,3);%plot ISI of each unit, see whether all >2ms
ind_peak=ind_peak_relative+ind_exceed_start;
set(gca,'xscale','log')
hold on;
for i=1:size(tb_group,1)
    idx_unit=tb_group(i,1);
    i_spks=(idx==idx_unit);
    ISI=diff(ind_peak(i_spks));
    histogram(ISI/fr,'EdgeColor',color_3d{i},'DisplayStyle','stairs');
end
y_lim=get(gca,'Ylim');
plot([0.002,0.002],y_lim,'k-');

subplot(2,3,4);
gscatter(score(:,1),score(:,2),idx,color_mat);
xlabel('PC1');
ylabel('PC2');
subplot(2,3,5);
gscatter(score(:,1),score(:,3),idx,color_mat);
xlabel('PC1');
ylabel('PC3');
subplot(2,3,6);
gscatter(score(:,2),score(:,3),idx,color_mat);
xlabel('PC2');
ylabel('PC3');

figure(figTrace);
hold on;
y_lim=get(gca,'Ylim');
y_label_pos=(0.1*y_lim(1)+0.9*y_lim(2));
y_label_pos2=(0.2*y_lim(1)+0.8*y_lim(2));
spike_unit=[];
for i=1:size(tb_group,1)    
    idx_unit=tb_group(i,1);
    ind_unit=ind_exceed_start(logical((idx==idx_unit).*flag_spk));
    plot(ts(ind_unit),y_label_pos*ones(size(ind_unit)),'o','Color',color_3d{i});
%     if idx_unit~=idx_noise
        spike_i=zeros(size(data));
        spike_i(ind_unit)=1;
        spike_unit=[spike_unit,spike_i];
%     end
end
spike_unit=logical(spike_unit);
spike_vec=logical(sum(spike_unit,2));
plot(ts(ind_exceed_start(logical(flag_spk))),y_label_pos2*ones(size(ind_exceed_start(logical(flag_spk)))),'or');
plot(ts(ind_exceed_start(logical(~flag_spk))),y_label_pos2*ones(size(ind_exceed_start(logical(~flag_spk)))),'ok');
mid_data.flag_spk=flag_spk;
mid_data.data_std=data_std;


end

