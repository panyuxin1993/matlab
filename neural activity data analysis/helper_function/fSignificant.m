function [significant_matrix,start_sig_mat] = fSignificant(signals,criteria,threshold_dur,varargin)
%fSignificant change signals/time series into logical 0/1 matrix of same
%size, based on the criteria
%Input-
%   signals - m-by-n matrix, m ROIs and n time points
%   criteria- char, eg. 2STD,
%   threshold_dur- number of time (unit) needed to be included as
%   significant
%   varargin- name-value pair
%       'BaselineIndex'-ones(size(signals))(default)| 1-0 logical
%       matrix of same size/vector along 2nd dimension, the time dimension

%set default value
indBaseline=true(size(signals));
criteria_off=[];
if ~isempty(varargin)
    if mod(length(varargin),2)~=0
        warning('Input Name-value pair for fSignificat function not match');
    else
        for i=1:2:length(varargin)
            switch varargin{i}
                case 'BaselineIndex'
                    if strcmp(varargin{i+1},'auto')
                        indBaseline=fIndexBaselineAroundMode(signals);
                    else
                        indBaseline=varargin{i+1};
                    end
                case 'DataForm'
                    activity_form=varargin{i+1};
                case 'EventOffCriteria'
                    criteria_off=varargin{i+1};
                otherwise
                    warning(['Invalid name-value pair named ',varargin{i}]);
            end
        end
    end
end
if size(indBaseline,1) == 1
    indBaseline=repmat(indBaseline,size(signals,1),1);
elseif size(indBaseline,1) ~= size(signals,1)
    warning('inconsistent size of BaselineIndex matrix and signal matrix');
    indBaseline=repmat(indBaseline(1,:),size(signals,1),1);
elseif size(indBaseline,2) ~= size(signals,2)
    warning('inconsistent size of BaselineIndex matrix and signal matrix, include all signals data as baseline');
end

[significant_matrix,significant_matrix_off]=deal(zeros(size(signals)));
if contains(criteria,'STD')%using several folds of standard deviation as criteral
    nfolds_str=strrep(criteria,'STD','');
    nfolds=str2double(nfolds_str);
    nfolds_off=nfolds;%by default, using same threshold for event off with event rising, for dff, may ref (Dombeck et al.,2007, neuron) 2std as rising threshold and 0.5std as event ending.
    if ~isempty(criteria_off) && contains(criteria_off,'STD')
        nfolds_off_str=strrep(criteria_off,'STD','');
        nfolds_off=str2double(nfolds_off_str);
    end        
    
    for iROI=1:size(signals,1)
        data_mean=nanmean(signals(iROI,indBaseline(iROI,:)),2);
        data_std=nanstd(signals(iROI,indBaseline(iROI,:)),0,2);
        significant_matrix(iROI,:)=(signals(iROI,:)>data_mean+nfolds*data_std)+(signals(iROI,:)<data_mean-nfolds*data_std);
        significant_matrix_off(iROI,:)=(signals(iROI,:)>data_mean+nfolds_off*data_std)+(signals(iROI,:)<data_mean-nfolds_off*data_std);
    end
    
elseif contains(criteria,'resampling')
    disp('The resampling method still under development');
end
%threshold for continues duration
for i=1:size(significant_matrix,1)
    j=1;
    while j<=size(significant_matrix,2)
        if significant_matrix(i,j)==1
            current_dur=1;
            while j+current_dur<=size(significant_matrix,2)
                if significant_matrix_off(i,j+current_dur)==1 %once above the threshold, then use offset threshold
                    current_dur=current_dur+1;
                else
                    break;
                end
            end
            if current_dur< threshold_dur
                significant_matrix(i,j:j+current_dur-1)=0;
            else
                significant_matrix(i,j:j+current_dur-1)=1;
            end
            j=j+current_dur;
        else
            j=j+1;
        end
    end
end
      
significant_matrix=logical(significant_matrix);
%find the start of significance series
diff_sig_mat=[zeros(size(significant_matrix,1),1), diff(significant_matrix,1,2)];
start_sig_mat=(diff_sig_mat==1);

end

