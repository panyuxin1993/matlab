%plot raster plot and PSTH for opto stimuli
close all;
rootpath='E:\ephys\20220705';
dirs=dir(rootpath);
dircell=struct2cell(dirs);
filenames=dircell(1,:);   
ind_tg=cellfun(@(x) contains(x,'2022'), filenames);
ind_tg=find(ind_tg);
n_recordings=length(ind_tg);
%know how many conditions in total
metadata_cell = cellfun(@fGetMetadata,filenames(ind_tg),'UniformOutput',true);
metadata_table= struct2table(metadata_cell);
Tmetadata = removevars(metadata_table,{'cell_info','date','time'});
Tuni_condition = unique(Tmetadata);
n_condition = size(Tuni_condition,1);
[mean_spk_PSTH,mean_LFP_PSTH,laser_dur_cell,ts_cell]=deal(cell(n_condition,1));%each cell a group of multiunit activities for similar stim condition 
%some global parameters
fr=30000;%Hz
savepath='E:\ephys\mid_data';
%
for i_file=1:n_recordings
    name_recording=filenames{ind_tg(i_file)};
    metadata_i=fGetMetadata(name_recording);
    ind_condition=logical((metadata_i.stim_dur == Tuni_condition.stim_dur).*...
        (metadata_i.power == Tuni_condition.power).*...
        strcmp(metadata_i.stim_type,Tuni_condition.stim_type));
    datapath=strcat(rootpath,filesep,name_recording);
    expID=strsplit(datapath,'\');
    threSTD=3;
    crit=3;
    mid_data_name=[savepath,filesep,expID{end},'para-movingthreSTD',num2str(threSTD),'crit',num2str(crit),'.mat'];
    mid_fig_name=[savepath,filesep,expID{end}];
    if exist(mid_data_name,'file')
        [data, timestamps, info]=load_open_ephys_data([datapath,'\Record Node 108\100_2.continuous']);
        load(mid_data_name);
    else
        [data, timestamps, info]=load_open_ephys_data([datapath,'\Record Node 108\100_2.continuous']);
        data_trial_start=load_open_ephys_data([datapath,'\Record Node 108\100_33.continuous']);
        data_laser=load_open_ephys_data([datapath,'\Record Node 108\100_34.continuous']);
        [LFPdata, timestampsLFP, infoLFP]=load_open_ephys_data([datapath,'\Record Node 111\100_2.continuous']);
        % plot(timestamps, data_laser);
        %from votage data to spike raster
        [spike_vec,spike_unit,mid_data,figPCA,figTrace] = fGetSpikes(data,timestamps,threSTD,crit);
        save(mid_data_name,'spike_vec','spike_unit','data','mid_data','timestamps','info','data_trial_start','data_laser','LFPdata');
        saveas(figPCA,[mid_fig_name,'-PCA.png'],'png');
        saveas(figTrace,[mid_fig_name,'-trace.png'],'png');
    end
    %aligned to laser event
    smooth_laser=smooth(data_laser,3);
    %plot(timestamps,data_laser);
    laser_on=smooth_laser>mode(smooth_laser);
    laser_on_start=[0;diff(laser_on)];
    laser_start=find(laser_on_start>0);
    laser_start=laser_start+1;%due to smooth, find the real time of laser on
    laser_end=find(laser_on_start<0);
    laser_end=laser_end-1;
    n_laser=min(length(laser_start),length(laser_end));
    laser_end=laser_end(1:n_laser);
    laser_start=laser_start(1:n_laser);
    laser_dur=laser_end-laser_start;
    ind_artifact=laser_dur<300;%<10ms
    laser_end=laser_end(~ind_artifact);
    laser_start=laser_start(~ind_artifact);
    laser_dur_mean=mean(laser_dur(~ind_artifact));
    if isempty( laser_dur_cell{ind_condition})%make same length of each cases
        laser_dur_cell{ind_condition}=laser_dur_mean;
        ts_cell{ind_condition}= -2:1/fr:2+laser_dur_mean/fr;
    else
        laser_dur_cell{ind_condition}=min(laser_dur_cell{ind_condition},laser_dur_mean);
        ts_cell{ind_condition}= -2:1/fr:2+laser_dur_cell{ind_condition}/fr;
    end
    %recheck if laser number is 20
    if length(laser_end) ~= 20
        warning('laser number ~= 20, check again');
    end
    timewindow=[timestamps(laser_start)-2,timestamps(laser_start)+2+laser_dur_cell{ind_condition}/fr];
    tw_aligned=timewindow-timewindow(:,1)-2;%x label, 0 means laser onset
    framewindow=[laser_start-fr*2,laser_start+laser_dur_cell{ind_condition}+fr*2];%ref settings for the exp
    
    %----------
    %plot multi-unit firing rate
    %----------
    %get data matrix for spike rater and LFP
    spk_unit_mat=[];
    data_mat=[];
    LFPdata_mat=[];
    for i=1:length(laser_end)
        spk_unit_mat=[spk_unit_mat,spike_vec(framewindow(i,1):framewindow(i,2))];
        data_mat=[data_mat,data(framewindow(i,1):framewindow(i,2))];
        LFPdata_mat=[LFPdata_mat,LFPdata(framewindow(i,1):framewindow(i,2))];
    end
    if size(mean_LFP_PSTH{ind_condition},1)==0
        n_min_frame=size(LFPdata_mat,1);
        mean_LFP_PSTH{ind_condition}=mean(LFPdata_mat,2);
        mean_spk_PSTH{ind_condition}=mean(spk_unit_mat,2);    
    else
        n_min_frame=min(size(mean_LFP_PSTH{ind_condition},1),size(LFPdata_mat,1));
        mean_LFP_PSTH{ind_condition}=cat(2,mean_LFP_PSTH{ind_condition}(1:n_min_frame,:),mean(LFPdata_mat(1:n_min_frame,:),2));
        mean_spk_PSTH{ind_condition}=cat(2,mean_spk_PSTH{ind_condition}(1:n_min_frame,:),mean(spk_unit_mat(1:n_min_frame,:),2));    
    end

    %plot multiunit spiking raster with raw data, LFP, and PSTH 
    %
    figRasterLFP=figure;
    set(gcf,'Position',[50,50,500,450]);
    suptitle([metadata_i.cell_info,'-power',num2str(metadata_i.power),'uw-',num2str(metadata_i.stim_dur),'-',metadata_i.stim_type]);
    subplot(3,2,[1,3]);%spike raster with raw data
    hold on;
    for i=1:length(laser_end)
        trial_vec=spike_vec(framewindow(i,1):framewindow(i,2));
        spk_train=find(trial_vec);
        spk_unit_mat=[spk_unit_mat,trial_vec];
        ts=tw_aligned(i,1):1/fr:tw_aligned(i,2);
        range=[min(data_mat,[],'all'),max(data,[],'all')];
        plot(ts,(i-1)*(range(2)-range(1))+data_mat(:,i)-range(1),'Color',[0.5,0.5,0.5]);
        plot(ts(spk_train),((i-1)*(range(2)-range(1))+mean(range))*ones(size(spk_train)),'b.');
    end
    %plot opto epochs using bar above
    set(gca,'Ylim',[0,length(laser_end)*(range(2)-range(1))]);
    y_lim=get(gca,'Ylim');
    plot([0,0],[y_lim(1),y_lim(2),],'k-');
    plot(laser_dur_mean/fr*ones(1,2),[y_lim(1),y_lim(2),],'k-');
    plot([0,laser_dur_mean/fr],y_lim(2)*1*ones(1,2),'b-','LineWidth',3);
    set(gca,'YTick',(range(2)-range(1))*10:(range(2)-range(1))*10:y_lim(2),'YTickLabel',{'10','20'});
    ylabel('Trial number');
    set(gca,'FontSize',12);

    subplot(3,2,[2,4]);%LFP
    hold on;
    for i=1:length(laser_end)
        ts=tw_aligned(i,1):1/fr:tw_aligned(i,2);
        range=[min(LFPdata_mat,[],'all'),max(LFPdata,[],'all')];
        plot(ts,(i-1)*(range(2)-range(1))+LFPdata_mat(:,i)-range(1),'k');
    end
    %plot opto epochs using bar above
    set(gca,'Ylim',[0,length(laser_end)*(range(2)-range(1))]);
    y_lim=get(gca,'Ylim');
    plot([0,0],[y_lim(1),y_lim(2),],'k-');
    plot(laser_dur_mean/fr*ones(1,2),[y_lim(1),y_lim(2),],'k-');
    plot([0,laser_dur_mean/fr],y_lim(2)*1*ones(1,2),'b-','LineWidth',3);
    set(gca,'YTick',(range(2)-range(1))*10:(range(2)-range(1))*10:y_lim(2),'YTickLabel',{'10','20'});
    ylabel('Trial number');
    set(gca,'FontSize',12);
    
    subplot(3,2,5);%PSTH of spikes
    hold on;
    spk_rate=mean(spk_unit_mat,2);
    binsize=0.1;%s
    spk_rate_bin=zeros(ceil(size(spk_rate,1)/(binsize*fr)),1);
    ts_bin=ts(0.5*binsize*fr:binsize*fr:size(spk_rate,1));
    j=1;
    for i=1:binsize*fr:size(spk_rate,1)
        ind2=min(size(spk_rate,1),i+binsize*fr-1);
        spk_rate_bin(j,1)=sum(spk_rate(i:ind2,1))/binsize;
        j=j+1;
    end
    n_bin=min(length(spk_rate_bin),length(ts_bin));
    plot(ts_bin(1:n_bin),spk_rate_bin(1:n_bin,1),'k-');
    %plot opto epochs using patch
    y_lim=get(gca,'Ylim');
    plot([0,0],[y_lim(1),y_lim(2),],'k-');
    plot(laser_dur_mean/fr*ones(1,2),[y_lim(1),y_lim(2),],'k-');
    plot([0,laser_dur_mean/fr],y_lim(2)*0.95*ones(1,2),'b-','LineWidth',3);
    xlabel('Laser onset (s)');
    ylabel('Firing rate (spikes/s)');
    set(gca,'FontSize',12);

    subplot(3,2,6);%PSTH of LFP
    hold on;
    LFPdata_mat_mean=mean(LFPdata_mat,2);
    binsize=0.01;%s
    LFPdata_bin=zeros(ceil(size(LFPdata_mat_mean,1)/(binsize*fr)),1);
    ts_bin=ts(0.5*binsize*fr:binsize*fr:size(LFPdata_mat_mean,1));
    j=1;
    for i=1:binsize*fr:size(LFPdata_mat_mean,1)
        ind2=min(size(LFPdata_mat_mean,1),i+binsize*fr-1);
        LFPdata_bin(j,1)=sum(LFPdata_mat_mean(i:ind2,1))/binsize;
        j=j+1;
    end
    n_bin=min(length(ts_bin),length(LFPdata_bin));
    plot(ts_bin(1:n_bin),LFPdata_bin(1:n_bin,1),'k-');
    %plot opto epochs using patch
    y_lim=get(gca,'Ylim');
    plot([0,0],[y_lim(1),y_lim(2),],'k-');
    plot(laser_dur_mean/fr*ones(1,2),[y_lim(1),y_lim(2),],'k-');
    plot([0,laser_dur_mean/fr],y_lim(2)*0.95*ones(1,2),'b-','LineWidth',3);
    xlabel('Laser onset (s)');
    ylabel('LFP');
    set(gca,'FontSize',12);
    %

    %--------------
    %plot multiunit spiking raster and PSTH, concise version
    %--------------
    %{
    figRaster=figure;
    set(gcf,'Position',[100,100,250,400]);
    temp=strsplit(expID{end},'_');
    suptitle([temp{1},'-',temp{end}]);
    subplot(2,1,1);%spike raster
    hold on;
    for i=1:length(laser_end)
        trial_vec=spike_vec(framewindow(i,1):framewindow(i,2));
        spk_train=find(trial_vec);
        spk_unit_mat=[spk_unit_mat,trial_vec];
        ts=tw_aligned(i,1):1/fr:tw_aligned(i,2);
        plot(ts(spk_train),i*ones(size(spk_train)),'k.');
    end

    %plot opto epochs using bar above
    y_lim=get(gca,'Ylim');
    plot([0,0],[y_lim(1),y_lim(2),],'k-');
    plot(laser_dur_mean/fr*ones(1,2),[y_lim(1),y_lim(2),],'k-');
    plot([0,laser_dur_mean/fr],y_lim(2)*1.1*ones(1,2),'b-','LineWidth',3);
    ylabel('Trial number');
    set(gca,'FontSize',12);
    
    subplot(2,1,2);%PSTH
    hold on;
    spk_rate=mean(spk_unit_mat,2);
    binsize=0.2;%s
    spk_rate_bin=zeros(ceil(size(spk_rate,1)/(binsize*fr)),1);
    ts_bin=ts(0.5*binsize*fr:binsize*fr:size(spk_rate,1));
    
    j=1;
    for i=1:binsize*fr:size(spk_rate,1)
        ind2=min(size(spk_rate,1),i+binsize*fr-1);
        spk_rate_bin(j,1)=sum(spk_rate(i:ind2,1))/binsize;
        j=j+1;
    end
    plot(ts_bin,spk_rate_bin(:,1),'k-');

    %plot opto epochs using patch
    y_lim=get(gca,'Ylim');
    plot([0,0],[y_lim(1),y_lim(2),],'k-');
    plot(laser_dur_mean/fr*ones(1,2),[y_lim(1),y_lim(2),],'k-');
    plot([0,laser_dur_mean/fr],y_lim(2)*1.1*ones(1,2),'b-','LineWidth',3);
    xlabel('Laser onset (s)');
    ylabel('Firing rate (spikes/s)');
    set(gca,'FontSize',12);
    %}

end
save([savepath,filesep,'grouped-data-20220705.mat'],'mean_spk_PSTH','mean_LFP_PSTH','laser_dur_cell','ts_cell');
%}
load([savepath,filesep,'grouped-data-20220705.mat']);
%% plot population data
binsize=0.1;%s
datatype_str='firing rate';
figMeanSpk= fPlotNormMean(mean_spk_PSTH,ts,binsize,fr,laser_dur_cell,Tuni_condition,datatype_str);

binsize=0.01;%s
datatype_str='LFP';
figMeanLFP= fPlotNormMean(mean_LFP_PSTH,ts,binsize,fr,laser_dur_cell,Tuni_condition,datatype_str);



function metadata = fGetMetadata(name_recording)
temp=strsplit(name_recording,'_');
metadata.cell_info=temp{1};
metadata.date=temp{2};
metadata.time=temp{3};
stim_info=temp{4};
temp2=strsplit(stim_info,'-');
powerstr=strrep(temp2{1},'power','');
powernum=str2num(strrep(powerstr,'uw',''));
metadata.power=powernum;
metadata.stim_dur=0.5;%default
metadata.stim_type='const_ramp';
if length(temp2)>1
    stim_dur_str=temp2{2};
    metadata.stim_dur=str2num(strrep(stim_dur_str,'s',''));
    if length(temp2)>2
        metadata.stim_type=temp2{3};
    end
end
end

function figout= fPlotNormMean(data,ts,binsize,fr,laser_dur_cell,Tmetadata,datatype_str)
%plot normalized firing rate across cells 
%binsize in s, fr frames/s
n_cells=cellfun(@(x) size(x,2), data,'UniformOutput',true);
n_col=6;
n_row=ceil(length(n_cells)/n_col);
ind_baseline=ts<0;
figout=figure;
set(gcf,'Position',[50,50,n_col*150, n_row*150]);
if strcmp(datatype_str,'LFP')
    y_str='Norm. LFP';
else
    y_str='Norm. firing rate';
end

for i_condition=1:length(n_cells)
    data_mat_mean=data{i_condition};
    data_mean=mean(data_mat_mean(ind_baseline,:),1);
    data_mat_norm=data_mat_mean./repmat(data_mean,size(data_mat_mean,1),1);
    data_bin=zeros(ceil(size(data_mat_norm,1)/(binsize*fr)),size(data_mat_norm,2));
    ts=-2:1/fr:2+laser_dur_cell{i_condition}/fr;
    ts_bin=ts(0.5*binsize*fr:binsize*fr:size(data_mat_norm,1));
    j=1;
    for i=1:binsize*fr:size(data_mat_norm,1)
        ind2=min(size(data_mat_norm,1),i+binsize*fr-1);
        data_bin(j,:)=sum(data_mat_norm(i:ind2,:),1)/binsize;
        j=j+1;
    end
    n_bin=min(length(ts_bin),size(data_bin,1));
    [data_bin_mean, data_bin_SE]=fMean_SE(data_bin(1:n_bin,:)');
    subplot(n_row,n_col,i_condition);
    xpatch=[ts_bin(1:n_bin),fliplr(ts_bin(1:n_bin))];
    ypatch=[data_bin_SE(1,:),fliplr(data_bin_SE(2,:))];
    p=patch(xpatch,ypatch,[0.5,0.5,0.5]);
    p.EdgeColor='none';
    hold on;
    plot(ts_bin(1:n_bin),data_bin_mean,'k-');
    %plot opto epochs using patch
    laser_dur_mean=laser_dur_cell{i_condition};
    y_lim=get(gca,'Ylim');
    plot([0,0],[y_lim(1),y_lim(2),],'k-');
    plot(laser_dur_mean/fr*ones(1,2),[y_lim(1),y_lim(2),],'k-');
    plot([0,laser_dur_mean/fr],y_lim(2)*0.95*ones(1,2),'b-','LineWidth',3);
    xlabel('Laser onset (s)');
    ylabel(y_str);
    power =Tmetadata.power;
    stim_type= Tmetadata.stim_type;
    titlestr=['power',num2str(power(i_condition)),'uw-',stim_type{i_condition}];
    title(strrep(titlestr,'_','\_'));
    text(0.9,0.2,['n=',num2str(n_cells(i_condition))],'Unit','Normalized');
    set(gca,'FontSize',12);
end

end

