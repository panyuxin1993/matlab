function [ fig ] = fPlotIHCinFig( fig,position )
%FPLOTIHCINFIG plot IHC result in bar graph into specific position in whole
% figure, 
%   input is figure to plot, position 
%   output is the input fig

n=[160,171];
figure(fig);
set(gcf,'PaperPosition',position);
c=bar(1,n(1)/sum(n),'y');hold on;
c=bar(2,n(2)/sum(n),'g');hold on;
set(gca,'XTick',[1,2],'XTickLabel',{'GABA+','GABA-'});
% text(position(1),position(2)+position(4)*0.9,panel_label,'FontName','Arial','FontSize',14);
box off;
set(gca,'FontName','Arial','FontSize',14);
ylabel('Proportion')


end

