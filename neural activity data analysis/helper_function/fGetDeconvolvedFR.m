function [denoised_trace,spiking_rate,zscored_spkr] = fGetDeconvolvedFR(SavedCaTrials,sessionname,saving_path, flag_plot,varargin)
%FGETDECONVOLVEDFR deconvolve using fluorescence signals
%Inputs:
%   f_raw: nROI x T  matrix, fluorescence trace, 
%   sessionname: animal_date_field, used as title and saving name of figure
%   varargin: other input to deconvolveCa function
%outputs:
%   denoised_trace: nROI x T matrix, denoised trace
%   spiking_rate: nROI xT matrix, deconvolved signal

%concanetate f signals together
f_cat  = [];

ind_tr_1=1;
ntr = length(SavedCaTrials.f_raw); %cancate all trials together, or  just set the number of trials to use
for i = ind_tr_1:ntr
    f_cat = [f_cat SavedCaTrials.f_raw{i}];
end

options=parseinputs();
nROI=size(f_cat,1);
denoised_trace=zeros(size(f_cat));
spiking_rate=zeros(size(f_cat));
zscored_spkr=zeros(size(f_cat));
for roiNo = 1:nROI
    [c, s] = deconvolveCa(f_cat(roiNo,:)');
    denoised_trace(roiNo,:)=c;
    spiking_rate(roiNo,:)=s;
end
%calculate z-scored firing rate
sr_mean=nanmean(spiking_rate,2);%mean for each roi
sr_std=nanstd(spiking_rate,0,2);%std for each roi
for roiNo = 1:nROI
    zscored_spkr(roiNo,:)=(spiking_rate(roiNo,:)-sr_mean(roiNo))/sr_std(roiNo);
end
%plot each trace
if strcmp(flag_plot,'plot_result')
    fig_deconv=figure();
    ncol=10;
    nrow=nROI/ncol;
    for roiNo = 1:nROI
        subplot(nrow,ncol,roiNo);
        plot(f_cat(roiNo,:),'k-');
        hold on;
        plot(denoised_trace(roiNo,:),'b-');
        plot(spiking_rate(roiNo,:),'r-');
        plot(zscored_spkr(roiNo,:),'y-');
        text(0.1,0.9,['ROI-',num2str(roiNo)],'Units','Normalized');
        if roiNo>ncol*(nrow-1)
            xlabel('frames')
        end
        if roiNo==1
            legend('f_raw','denoised trace','spiking rate','z-scored spiking rate');
        end
    end
    suptitle(sessionname);
    saveas(fig_deconv,[saving_path, filesep,sessionname, filesep,sessionname, '_deconvolution.jpg'],'jpg');
end
save([saving_path, filesep,'deconvolution.mat'],  'denoised_trace','spiking_rate','zscored_spkr');
end
function options=parseinputs(varargin)
%% parse input variables

%% default options
options.type = 'ar2';
options.pars = [];
options.sn = [];
options.b = 0;
options.lambda = 0;
options.optimize_b = false;
options.optimize_pars = false;
options.optimize_smin = false;
options.method = 'constrained';
options.window = 200;
options.shift = 100;
options.smin = 0;
options.maxIter = 10;
options.thresh_factor = 1.0;
options.extra_params = [];

%% parse all input arguments
k = 1;
while k<=nargin
    
    switch lower(varargin{k})
        case {'ar1', 'ar2', 'exp2', 'kernel'}
            % convolution kernel type
            options.type = lower(varargin{k});
            if (k<nargin) && (isnumeric(varargin{k+1}))
                options.pars = varargin{k+1};
                k = k + 1;
            end
            k = k + 1;
        case 'pars'
            % parameters for the kernel
            options.pars = varargin{k+1};
            k = k+2;
        case 'sn'
            % noise
            options.sn = varargin{k+1};
            k = k+2;
        case 'b'
            % baseline
            options.b = varargin{k+1};
            k = k+2;
        case 'optimize_b'
            % optimize the baseline
            options.optimize_b = true;
            k = k+1;
        case 'optimize_pars'
            % optimize the parameters of the convolution kernel
            options.optimize_pars = true;
            k = k+1;
        case 'optimize_smin'
            % optimize the parameters of the convolution kernel
            options.optimize_smin = true;
            k = k+1;
        case 'lambda'
            % penalty
            options.lambda = varargin{k+1};
            k = k+2;
        case {'foopsi', 'constrained', 'thresholded', 'mcmc'}
            % method for running the deconvolution
            options.method = lower(varargin{k});
            k = k+1;
            if strcmpi(options.method, 'mcmc') && (k<=length(varargin)) && (~ischar(varargin{k}))
                options.extra_params = varargin{k};
                k = k+1;
            end
        case 'window'
            % maximum length of the kernel
            options.window = varargin{k+1};
            k = k+2;
        case 'shift'
            % number of frames by which to shift window from on run of NNLS
            % to the next
            options.shift = varargin{k+1};
            k = k+2;
        case 'smin'
            % number of frames by which to shift window from on run of NNLS
            % to the next
            options.smin = varargin{k+1};
            k = k+2;
        case 'maxiter'
            % number of frames by which to shift window from on run of NNLS
            % to the next
            options.maxIter = varargin{k+1};
            k = k+2;
        case 'thresh_factor'
            % number of frames by which to shift window from on run of NNLS
            % to the next
            options.thresh_factor = varargin{k+1};
            k = k+2;
        otherwise
            k = k+1;
    end
end
end