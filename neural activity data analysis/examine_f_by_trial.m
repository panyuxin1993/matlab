inROI=[1,1];
indTrial=130:136;%35:39;
savepath='F:\FP\example';
figExample=fPlotF_ROI('F:\FP\pyx241_20191201','dff','FP',inROI,indTrial);%example for vglut2
set(figExample,'PaperPosition',[0,0,2,2]);
saveas(figExample,[savepath,filesep,'vglut2 example.pdf'],'pdf');

% inROI=[2,2];
% indTrial=11:20;
% savepath='F:\FP\example';
% figExample=fPlotF_ROI('F:\FP\pyx237_20191130','dff','FP',inROI,indTrial);%example for vglut2
% set(figExample,'PaperPosition',[0,0,2.5,2]);
% saveas(figExample,[savepath,filesep,'vgat example.pdf'],'pdf');

inROI=[1,1];
indTrial=3:9;
savepath='F:\FP\example';
figExample=fPlotF_ROI('F:\FP\pyx214_20190925','dff','FP',inROI,indTrial);%example for vglut2
set(figExample,'PaperPosition',[0,0,2,2]);
saveas(figExample,[savepath,filesep,'vgat example.pdf'],'pdf');

% %case used in 2020 progress report
% inROI=13:15;
% indTrial=15;%in 11:20
% figExample=fPlotF_ROI('H:\2P\pyx290_20200528\im_data_reg\result_save','dff','2P',inROI,indTrial);%example for 2P contralateral preference during delay
% set(figExample,'PaperPosition',[0,0,4,2.5]);
% savepath='F:\2P\example';
% saveas(figExample,[savepath,filesep,'contralteral delay selectivity example.pdf'],'pdf');