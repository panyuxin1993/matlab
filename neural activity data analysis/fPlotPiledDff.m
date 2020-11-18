function [figPiled] = fPlotPiledDff(neuralActivity,ts,indTrial,curvecolor)
%FPLOTPILEDDFF plot activities piled together one figure.
%Input-
%   neuralActivity-n-by-m matrix, n trials, each m frames
%   ts- time stamps, 1-by-m;
%   indTrial- index of trial to plot, if is empty, then plot all trials.
%   xlabelstr- string for xlabel
if isempty(indTrial)
    indTrial=1:size(neuralActivity,1);
end
figPiled=figure;
range=[];
for i=1:length(indTrial)
    current_f=neuralActivity(indTrial(i),:);
    current_f=current_f-min(current_f)+sum(range);
    plot(ts,current_f,'color',curvecolor,'LineWidth',1.5);
    hold on;
    xlim=get(gca,'Xlim');
    range=[range,max(current_f)-min(current_f)];
    text(xlim(end),sum(range),strcat('Trial-',num2str(indTrial(i))));
end
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[]);
% set(gca,'xtick',[ceil(ts(1)):0.5:ts(2)],'xticklabel',[ceil(ts(1)):0.5:ts(2)]);
ylim=get(gca,'Ylim');
plot([0,0],[ylim(1),ylim(2)],'k--','LineWidth',1);
%plot scale bar and text
if length(range)>=10
    scalebar_lengthy=fRound(prctile(range,10));
else
    scalebar_lengthy=0.1;
end
scalebar_lengthx=fRound((ts(end)-ts(1))/5);
%n_digit=fix(log10(scalebar_lengthy));
%scalebar_lengthy=round(scalebar_lengthy/10^n_digit)*10^n_digit;
xlim=get(gca,'Xlim');
plot([xlim(end),xlim(end)],[0,scalebar_lengthy],'k-','LineWidth',2);
text(xlim(end),scalebar_lengthy/2,['\it\DeltaF/F ','\rm',num2str(scalebar_lengthy)],'FontName','Arial','FontSize',10);
plot(xlim(end)-[0,scalebar_lengthx],zeros(2,1),'k-','LineWidth',2);
text(xlim(end)-scalebar_lengthx/2,0,[num2str(scalebar_lengthx),'s'],'FontName','Arial','FontSize',10);
end

%fuction for getting a round number, example: 166-->100, 25-->10
function [out]=fRound(in)
    n_digit=floor(log10(in));
    out=round(in/10^n_digit)*10^n_digit;
end
