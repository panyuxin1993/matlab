%% Example Title
% Summary of example objective
rootpath='C:\Data';
raw88=imread([rootpath,filesep,'pyx397_20211102\pyx397_20211102_deep20um_920nm_power30_12x_00002.tif'],88);
raw89=imread([rootpath,filesep,'pyx397_20211102\pyx397_20211102_deep20um_920nm_power30_12x_00002.tif'],89);
pre88=imread([rootpath,filesep,'pyx397_20211102\corrected\pre_warp\pyx397_20211102_deep20um_920nm_power30_12x_00002_corrected.tif'],88);
pre89=imread([rootpath,filesep,'pyx397_20211102\corrected\pre_warp\pyx397_20211102_deep20um_920nm_power30_12x_00002_corrected.tif'],89);
post8_88=imread([rootpath,filesep,'pyx397_20211102\corrected\post_warp_blocksize8\pyx397_20211102_deep20um_920nm_power30_12x_00002_corrected_warped.tif'],88);
post8_89=imread([rootpath,filesep,'pyx397_20211102\corrected\post_warp_blocksize8\pyx397_20211102_deep20um_920nm_power30_12x_00002_corrected_warped.tif'],89);
post15_88=imread([rootpath,filesep,'pyx397_20211102\corrected\post_warp_blocksize15\pyx397_20211102_deep20um_920nm_power30_12x_00002_corrected_warped.tif'],88);
post15_89=imread([rootpath,filesep,'pyx397_20211102\corrected\post_warp_blocksize15\pyx397_20211102_deep20um_920nm_power30_12x_00002_corrected_warped.tif'],89);
reg88=imread([rootpath,filesep,'pyx397_20211102\im_data_reg\pyx397_20211102_deep20um_920nm_power30_12x_00dftReg_002.tif'],88);
reg89=imread([rootpath,filesep,'pyx397_20211102\im_data_reg\pyx397_20211102_deep20um_920nm_power30_12x_00dftReg_002.tif'],89);
%% Section 1 plot example frame mismatch
% Description of first code block
figCmp=figure;
set(gcf,'Position',[100,100,800,800]);
subplot(3,2,1);
obj=imshowpair(raw88,raw89,'Scaling','joint');
title('Raw');
% legend('frame88','frame89');
subplot(3,2,3);
obj=imshowpair(pre88,pre89,'Scaling','joint');
title('Pre-warp');
subplot(3,2,4);
obj=imshowpair(reg88,reg89,'Scaling','joint');
title('Lab code');

subplot(3,2,5);
obj=imshowpair(post8_88,post8_89,'Scaling','joint');
title('Post-warp-blocksize8');
subplot(3,2,6);
obj=imshowpair(post15_88,post15_89,'Scaling','joint');
title('Post-warp-blocksize15');
saveas(figCmp,[rootpath,filesep,'pyx397_20211102\corrected\visualizeExamplePatchWarp.pdf'],'pdf');
save([rootpath,filesep,'pyx397_20211102\corrected\visualizeExamplePatchWarp.mat']);
%% Section 2 Title
% Description of second code block
