function [spiking_rate, denoised_trace, zscored_spkr,ROIFitCoefs] = fxy_GetDeconvolvedFR(SessionTraceData,Paras,sessionname,saving_path, flag_plot)
%FXY_GETDECONVOLVEDFR from dff to get deconvolved firing rate, code by Xin,
%Yu based on CalmAn tool box
%   

[spkSNR,lamPr,fr,DecayTime] = deal(Paras{:});
%default parameters
if isempty(spkSNR)
    spkSNR = 0.5;
end
if isempty(lamPr)
    lamPr = 0.99;
end
if isempty(fr)
    fr = 30;
end
if isempty(DecayTime)
    DecayTime = 2;
end

[spiking_rate, denoised_trace, ROIFitCoefs] = MtxTrace2SP(SessionTraceData,{spkSNR,lamPr,fr,DecayTime});

%smooth the spiking rate and calculated a z-scoredfiring rate
nROI=size(SessionTraceData,1);
zscored_spkr=zeros(size(spiking_rate));
for roiNo = 1:nROI
    smoothed_sr=smooth(spiking_rate(roiNo,:),3);%smooth the data
    sr_mean=nanmean(smoothed_sr);%mean for each roi
    sr_std=nanstd(smoothed_sr);%std for each roi
    zscored_spkr(roiNo,:)=(smoothed_sr-sr_mean)/sr_std;
end
save([saving_path, filesep,'deconvolution.mat'],  'denoised_trace','spiking_rate','zscored_spkr');

%plot each trace
if strcmp(flag_plot,'plot_result')
    fig_deconv=figure();
    ncol=8;
    nrow=ceil(nROI/ncol);
    for roiNo = 1:size(SessionTraceData,1)
        subplot(nrow,ncol,roiNo);
        curve_dff=plot(SessionTraceData(roiNo,:),'k-');%dff
        hold on;
        curve_dn=plot(denoised_trace(roiNo,:),'b-');
        curve_sr=plot(spiking_rate(roiNo,:),'r-');
        curve_zssr= plot(zscored_spkr(roiNo,:),'y-');
        text(0.1,0.9,['ROI-',num2str(roiNo)],'Units','Normalized');
        if roiNo>ncol*(nrow-1)
            xlabel('frames')
        end
        if roiNo==1
            legend([curve_dff,curve_dn,curve_sr,curve_zssr],{'f_raw','denoised trace','spiking rate','z-scored spiking rate'});
        end
    end
    titlestr=strrep(sessionname,'_','\_');
    suptitle(titlestr);
    if ~exist([saving_path, filesep,sessionname])
        mkdir([saving_path, filesep,sessionname]);
    end
    saveas(fig_deconv,[saving_path, filesep,sessionname, filesep,sessionname, '_deconvolution.jpg'],'jpg');
end
end

