function [outputcurve] = fPSTHBinDelayEndPoint(neuralActivity,color,varargin)
%FPSTHBINDELAYENDPOINT plot PSTH(mean+individual) with CI, middle points, 
%especially useful for sessions with variable delay length
% Input- neuralActivity, idealy is a m-by-n matrix,where m is trial
%   number, n is frames of each trial, each row have variable nan
%   values,which indicates the variable delay(end of delay is start of nan)
%   color- the mean trace color of the PSTH to be plotted
%   varargin can be these name-value pairs
%   'Style' = 'meanCIEnd'(default)|'meanCI'|'meanEnd', indicating whether
%   to plot CI or end point lines
%   'PCI'=0.05(default)|0.01|0.001 or other values indicating the p value of CI
%   'EndSmooth'='smooth'(default)|'raw', indicating whether smooth the
%   end point curve
%   'EndSmoothMethod'='bin'(default)|'span', indicating the method to smooth
%   'ParaEndSmoothMethod'= 3(defalt for 'bin' method),
%   indicating the parameter need for corresponding method, see fSmooth in fGetEndDff

% Output- outputcurve is the handle of curve, that need for legend
%using default value
PCI=0.05;
EndSmooth='smooth';
Style='meanCIEnd';
if nargin>2
    nvar=nargin-2;
    if mod(nvar,2)==1
        warning('input variables number error');
    end
    for i=1:nvar/2
        switch varargin{i*2-1}
            case 'PCI'
                PCI=varargin{i*2};
            case 'Style'
                Style=varargin{i*2};
            case 'EndSmooth'
                EndSmooth=varargin{i*2};
            case 'EndSmoothMethod'
                EndSmoothMethod=varargin{i*2};
                if strcmp(varargin{i*2+1},'ParaEndSmoothMethod')
                    ParaEndSmoothMethod=varargin{i*2+2};
                    i=i+2;
                else
                    warning('input endSmoothMethod should be followed by paraEndSmoothMethod');
                end
            otherwise
        end
    end
end
[neuralActivityMean, neuralActivityCI]=fMean_CI(neuralActivity,PCI);
%plot individual traces
color_case=(1+color)/2;
plot(1:size(neuralActivity,2),neuralActivity,'Color',color_case,'linewidth',1);
hold on;
%shadow as ci
if contains(Style,'CI')
    xpatch=[1:size(neuralActivity,2), fliplr(1:size(neuralActivity,2))];
    ypatch=[neuralActivityCI(1,:),fliplr(neuralActivityCI(2,:))];
    p=patch(xpatch,ypatch,color);%plot confidence interval
    p.FaceAlpha=0.1;
    p.EdgeColor=color;%'none';
    hold on;
end
outputcurve=plot(1:size(neuralActivity,2),neuralActivityMean,'Color',color,'linewidth',2);
hold on;

%plot endpoint dff
if contains(Style,'End')
    [ ~,~,neuralActivityGroupedtemp ] = fGetEndDff( neuralActivity ,EndSmooth,EndSmoothMethod,ParaEndSmoothMethod);
    [ timepoint,endDff,neuralActivityGrouped ] = fGetEndDff( neuralActivityGroupedtemp ,'raw');
    for i=1:size(neuralActivityGrouped,1)
        plot(neuralActivityGrouped(i,:),'Color',color,'linewidth',1.5);
    end
    scatter(timepoint,endDff,50,color,'MarkerFaceColor','w');
end
end

