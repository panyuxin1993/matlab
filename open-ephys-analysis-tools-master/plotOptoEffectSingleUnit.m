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
metadata_cell = cellfun(@fGetMetadata,name_recording(ind_tg),'UniformOutput',1);
metadata_table= struct2table(metadata_cell);
Tmetadata = removevars(metadata_table,'cell_info');
Tuni_condition = unique(Tmetadata);
n_condition = size(Tuni_condition,1);
[mean_spk_PSTH,mean_LFP_PSTH]=deal(cell(n_condition,1));%each cell a group of multiunit activities for similar stim condition 
for i_file=1:n_recordings
    name_recording=filenames{ind_tg(i_file)};
    metadata_i=fGetMetadata(name_recording);
    ind_condition=logical(strcmp(metadata_i.date,Tuni_condition.date).*...
        (metadata_i.stim_dur == Tuni_condition.stim_dur).*
        (metadata_i.power == Tuni_condition.power).*
        strcmp(metadata_i.stim_type,Tuni_condition.stim_type));
    datapath=strcat(rootpath,filesep,name_recording);
    savepath='E:\ephys\mid_data';
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
    %recheck if laser number is 20
    if length(laser_end) ~= 20
        warning('laser number ~= 20, check again');
    end
    fr=30000;%Hz
    timewindow=[timestamps(laser_start)-2,timestamps(laser_start)+2+laser_dur_mean/fr];
    tw_aligned=timewindow-timewindow(:,1)-2;%x label, 0 means laser onset
    framewindow=[laser_start-fr*2,laser_start+laser_dur_mean+fr*2];%ref settings for the exp
    

    %----------
    %plot each unit firing rate, but note the PCA results may not good
    %----------
    %
    spk_mat=[];
    figRaster=figure;
    set(gcf,'Position',[100,100,250,400]);
    subplot(2,1,1);%spike raster
    hold on;
    n_unit=size(spike_unit,2);
    color_unit=cell(n_unit,1);
    for i=1:n_unit
        color_unit{i}=[1-i/n_unit,0,i/n_unit];
    end
    for i_unit=1:size(spike_unit,2)
        spk_unit_mat=[];
        for i=1:length(laser_end)
            trial_vec=spike_unit(framewindow(i,1):framewindow(i,2),i_unit);
            spk_train=find(trial_vec);
            spk_unit_mat=[spk_unit_mat,trial_vec];
            ts=tw_aligned(i,1):1/fr:tw_aligned(i,2);
            plot(ts(spk_train),i*ones(size(spk_train))+(i_unit-1)*length(laser_end),'.','Color',color_unit{i_unit});
        end
        spk_mat=cat(3,spk_mat,spk_unit_mat);
    end
    %plot opto epochs using bar above
    y_lim=get(gca,'Ylim');
    plot([0,0],[y_lim(1),y_lim(2),],'k-');
    plot([1,1],[y_lim(1),y_lim(2),],'k-');
    plot([0,1],y_lim(2)*1.1*ones(1,2),'b-','LineWidth',3);
    ylabel('Trial number');
    
    subplot(2,1,2);%PSTH
    hold on;
    spk_rate=squeeze(mean(spk_mat,2));
    binsize=0.2;%s
    spk_rate_bin=zeros(ceil(size(spk_rate,1)/(binsize*fr)),size(spk_rate,2));
    ts_bin=ts(1+0.5*binsize:binsize*fr:size(spk_rate,1));
    for i_unit=1:size(spike_unit,2)
        j=1;
        for i=1:binsize*fr:size(spk_rate,1)
            ind2=min(size(spk_rate,1),i+binsize*fr-1);
            spk_rate_bin(j,i_unit)=sum(spk_rate(i:ind2,i_unit))/binsize;
            j=j+1;
        end
        plot(ts_bin,spk_rate_bin(:,i_unit),'-','Color',color_unit{i_unit});
    end
    %plot opto epochs using patch
    y_lim=get(gca,'Ylim');
    plot([0,0],[y_lim(1),y_lim(2),],'k-');
    plot([1,1],[y_lim(1),y_lim(2),],'k-');
    plot([0,1],y_lim(2)*1.1*ones(1,2),'b-','LineWidth',3);
    xlabel('Laser onset (s)');
    ylabel('Firing rate (spikes/s)');
    %}
end




function metadata = fGetMetadata(name_recording)
temp=strsplit(name_recording,'_');
metadata.cell_info=temp{1};
metadata.date=temp{2};
metadata.time=temp{3};
stim_info=temp{4};
temp2=strsplit(stim_info,'-');
powerstr=strrep(temp2{1},'power','');
powernum=str2num(strrep(powerstr,'um',''));
metadata.power=powernum;
metadata.stim_dur=0.5;%default
metadata.stim_type='const_ramp';
if length(temp2)>1
    stim_dur_str=temp2{2};
    metadata.stim_dur=str2num(strrep(stim_dur_str,'s',''));
    if length(temp2)>2
        metadata.stim_type=temp{3};
    end
end
end

function figout= fPlotNormMean(data,ts,binsize,fr,laser_dur_cell,datatype_str)
%plot normalized firing rate across cells 
%binsize in s, fr frames/s
n_cells=cellfun(@(x) size(x,2), data,'UniformOutput',true);
n_col=6;
n_row=ceil(n_cells/n_col);
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
    data_mean=mean(data_mat_mean(ind_baseline),1);
    data_mat_norm=(data_mat_mean-repmat(data_mean,size(data_mat_mean,1),1))./repmat(data_mean,size(data_mat_mean,1),1);
    data_bin=zeros(ceil(size(data_mat_norm,1)/(binsize*fr)),1);
    ts_bin=ts(0.5*binsize*fr:binsize*fr:size(data_mat_norm,1));
    j=1;
    for i=1:binsize*fr:size(data_mat_norm,1)
        ind2=min(size(data_mat_norm,1),i+binsize*fr-1);
        data_bin(j,1)=sum(data_mat_norm(i:ind2,1))/binsize;
        j=j+1;
    end
    n_bin=min(length(ts_bin),length(data_bin));
    subplot(n_row,n_col,i_condition);
    plot(ts_bin(1:n_bin),data_bin(1:n_bin,1),'k-');
    %plot opto epochs using patch
    laser_dur_mean=laser_dur_cell{i_condition};
    y_lim=get(gca,'Ylim');
    plot([0,0],[y_lim(1),y_lim(2),],'k-');
    plot(laser_dur_mean/fr*ones(1,2),[y_lim(1),y_lim(2),],'k-');
    plot([0,laser_dur_mean/fr],y_lim(2)*0.95*ones(1,2),'b-','LineWidth',3);
    xlabel('Laser onset (s)');
    set(gca,'FontSize',12);
end

end

