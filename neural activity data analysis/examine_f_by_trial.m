% inROI=[1,1];
% indTrial=130:136;%35:39;
% savepath='F:\FP\example';
% figExample=fPlotF_ROI('F:\FP\pyx241_20191201','dff','FP',inROI,indTrial);%example for vglut2
% set(figExample,'PaperPosition',[0,0,2,2]);
% saveas(figExample,[savepath,filesep,'vglut2 example.pdf'],'pdf');

% inROI=[2,2];
% indTrial=11:20;
% savepath='F:\FP\example';
% figExample=fPlotF_ROI('F:\FP\pyx237_20191130','dff','FP',inROI,indTrial);%example for vglut2
% set(figExample,'PaperPosition',[0,0,2.5,2]);
% saveas(figExample,[savepath,filesep,'vgat example.pdf'],'pdf');

% inROI=[1,1];
% indTrial=3:9;
% savepath='F:\FP\example';
% figExample=fPlotF_ROI('F:\FP\pyx214_20190925','dff','FP',inROI,indTrial);%example for vglut2
% set(figExample,'PaperPosition',[0,0,2,2]);
% saveas(figExample,[savepath,filesep,'vgat example.pdf'],'pdf');

% %case used in 2020 progress report
% inROI=13:15;
% indTrial=15;%in 11:20
% figExample=fPlotF_ROI('H:\2P\pyx290_20200528\im_data_reg\result_save','dff','2P',inROI,indTrial);%example for 2P contralateral preference during delay
% set(figExample,'PaperPosition',[0,0,4,2.5]);
% savepath='F:\2P\example';
% saveas(figExample,[savepath,filesep,'contralteral delay selectivity example.pdf'],'pdf');

% %case used in 2020 annual meeting
% inROI=[3,35,56];%1:103;
% indTrial=41:50;%in 11:20
% figExample=fPlotF_ROI('H:\2P\pyx285_20200520\im_data_reg\result_save','dff','2P',inROI,indTrial);%example for 2P contralateral preference during delay
% set(figExample,'PaperPosition',[0,0,4,2.5]);
% savepath='H:\2P\example';
% saveas(figExample,[savepath,filesep,'contralteral delay selectivity example.pdf'],'pdf');
%     %zoom in
% indTrial=49;%in 11:20
% figExample=fPlotF_ROI('H:\2P\pyx285_20200520\im_data_reg\result_save','dff','2P',inROI,indTrial);%example for 2P contralateral preference during delay
% set(figExample,'PaperPosition',[0,0,4,2.5]);
% savepath='H:\2P\example';
% saveas(figExample,[savepath,filesep,'contralteral delay selectivity example-zoom in.pdf'],'pdf');

%% try cases
inROI=1:3:50;%1:103;
indTrial=[];%in 11:20
% file_path='E:\2P\pyx326_20201223\im_data_reg\result_save';
% file_path='E:\2P\CD058_20180207\im_data_reg\result_save';
file_path='E:\2P\pyx397_20210928\im_data_reg\result_save';
figExample=fPlotF_ROI(file_path,'dff','2P',inROI,indTrial,'sigThreshSTD',2,'sigShowingStyle','patch','BaselineIndex','auto');%example for 2P contralateral preference during delay
figExample=fPlotF_ROI(file_path,'spkr','2P',inROI,indTrial,'sigThreshSTD',3,'sigShowingStyle','patch');%example for 2P contralateral preference during delay

