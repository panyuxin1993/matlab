%finding the target image for alignment
%cd('//Users/cduan/Desktop/Xu_Lab/Imaging/data/CD058/20180201');
cd('D:\pyx338_20210203');

[im, ~] = load_scim_data('pyx338_20210203_deep50um_920nm_power70_3X_001.tif');  
% [im, ~] = load_scim_data_pyx('pyx298_20200608_920nm_power50_3X_001.tif',[],0,2);  
colormap gray;
selectframe=im(:,:,250:320);   % change here for different target frame
imagesc(mean(selectframe,3),[-50 200]);
im_reg_target = mean(selectframe,3); 

%%
dir_imreg_src = pwd;
% save TargetImage.mat im_reg_target
dir_imreg_dest = [dir_imreg_src filesep 'im_data_reg'];
t_total=tic;
dft_reg_dir_2_zy(dir_imreg_src, dir_imreg_dest, [], im_reg_target)
% dft_reg_dir_2_pyx(dir_imreg_src, dir_imreg_dest, [], im_reg_target,[],1,2);
%[im_dft, shift] = dft_reg(20160125_rf_242.28_p100%_5.9x_001dir_imreg_src, im_reg_target);
t=toc(t_total);
disp(t);
save TargetImage.mat im_reg_target

%% save all the figures of shift 
h = get(0,'children');
for i=1:length(h)
  saveas(h(i), ['shift in trials' num2str(i)], 'png');
end