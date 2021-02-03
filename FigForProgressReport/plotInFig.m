rootpath='C:\Users\PYX\Documents\MATLAB\FigForProgressReport\';
% fig11=figure;
% fig11=fPlotBehaviorExample(fig11,[0.1,0.2,8,4]);
% print(fig11,[rootpath,'Behavior example-1A'],'-dpdf');
% fig12=figure;
% fig12=fPlotCNOSummary(fig12,[0.1,0.1,8,6],'ip');
% print(fig12,[rootpath,'M2 inhibition-1B'],'-dpdf');
% fig13=figure;
% fig13=fPlotCNOSummary(fig13,[0.1,0.1,8,6],'infusion');
% print(fig13,[rootpath,'M2-SC inhibition-1D'],'-dpdf');
% fig14=figure;
% fig14=fPlotCNOCurve(fig14,[0.1,0.1,7,2.5],'infusion');
% print(fig14,[rootpath,'M2-SC inhibition curve-1E'],'-dpdf');
% fig15=figure;
% fig15=fPlotCNOExample(fig15,[0.1,0.1,3,3],'ip');
% print(fig15,[rootpath,'M2 inhibition example-1F'],'-dpdf');
% fig16=figure;
% fig16=fPlotCNOExample(fig16,[0.1,0.1,3,3],'infusion');
% print(fig16,[rootpath,'M2-SC inhibition example-1G'],'-dpdf');
%  fig17=figure;
%  fig17=fPlotCNOSummaryError(fig17,[0.1,0.1,6,6],'infusion');
% print(fig17,[rootpath,'M2-SC inhibition summary by animal-1H'],'-dpdf');
% fig18=figure;
% fig18=fPlotCNOSummaryError(fig18,[0.1,0.1,6,6],'ip');
% print(fig18,[rootpath,'M2 inhibition summary by animal-1I'],'-dpdf');
% fig22=figure;
% fig22=fPlotImagingExample(fig22,[0.1,0.1,4,6]);
% print(fig22,[rootpath,'Example imaging cell-2B'],'-dpdf');
fig33=figure;
fig33=fPlotIHCinFig(fig33,[1,1,3,3]);
% print(fig33,[rootpath,'Cell type of SC neurons-3C'],'-dpdf')
fig34=figure;
fig34=fPlotEphysExample(fig34,[4,1,4,3]);
% saveas(fig34,[rootpath,'Cell type of SC neurons-3D.pdf'],'pdf');

% fig35=figure;
% fig35=fPlotEphysSummary(fig35,[1,1,4,3]);
% print(fig35,[rootpath,'Cell type of SC neurons-3E'],'-dpdf')
% fig42=figure;
% fig42=fPlotFPExample(fig42,[0.1,0.1,4,6]);
% print(fig42,[rootpath,'Example FP site-4B'],'-dpdf');

% saveas(fig3,'Cell type of SC neurons.pdf','pdf');
% close all;