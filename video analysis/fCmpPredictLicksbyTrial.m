function [ str_out ] = fCmpPredictLicksbyTrial( indTrial,stimOnset,delayOff,legendstr,color,varargin )
%FCMPPREDICTLICKS compare licks predicted by different methods
%   input-index of trial number; variable licking data from different methods
figure(1);
clf;
n_method=length(varargin);
if n_method==0
    warning('No available data of predicted lickings');
elseif length(legendstr)~=n_method
    warning('legend string length should be consistent with number of input data');
elseif length(color)~=n_method
    warning('color length should be consistent with number of input data');
elseif length(stimOnset)~=n_method ||length(delayOff)~=n_method
    warning('behavioral event vectors length should be consistent with number of input data');
else
    str_out=['Trial', num2str(indTrial),' of ', num2str(n_method),'methods'];
    %plot delay duration
    x_delay=[stimOnset(indTrial),delayOff(indTrial),delayOff(indTrial),stimOnset(indTrial)];
    y_delay=[0,0,n_method,n_method];
    p=patch(x_delay,y_delay,0.5*ones(1,3));
    p.FaceAlpha=0.1;
    hold on;
    %plot licks predicted by different methods
    for i_method=1:n_method
        lick=varargin{i_method}{indTrial};
        lick(isnan(lick))=[];%when combining left and right lick together, if either is nan, then cause something like [1,2,3,nan]
        if sum(size(lick)>1)>1
            warning('Each element of input cell should be a vector rather than a matrix');
        elseif isempty(lick)
            curve(i_method)=plot([0,0],[i_method-1,i_method],'Color',color{i_method},'LineWidth',1);
            hold on;
        else
            for i_lick=1:length(i_method)
                curve(i_method)=plot([lick(i_lick),lick(i_lick)],[i_method-1,i_method],'Color',color{i_method},'LineWidth',1);
                hold on;
            end
        end
    end

    xlabel('time(s) from trial start');
    h=legend(curve,legendstr);
    set(h,'box','off');   
    set(gca,'FontName','Arial','FontSize',12);
    title(str_out);
end



end

