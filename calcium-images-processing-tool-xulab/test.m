%%
filename = 'E:\data_of_different_brian_area\xinyu\20180910_b55a0706\b55a07_test06_2x_2afc_160um_20180910_dftReg_001.tif';
info = imfinfo(filename);
Key = 'scanFrameRate = ';
Str = info(1).ImageDescription;
Index = strfind(Str, Key);
startIndex = Index + length(Key);
endIndex = strfind(Str(startIndex:end), newline);
endIndex = endIndex(1);
Value = str2double(Str(startIndex: startIndex + endIndex -1))

%%
name = categorical({'mrcnn', 'STNeuroNet', 'HNCcorr + conv2d', 'Sourcery (Suite2p)', 'Unet2DS'});
name = reordercats(name, {'mrcnn', 'STNeuroNet', 'HNCcorr + conv2d', 'Sourcery (Suite2p)', 'Unet2DS'});
f_score = [0.692, 0.688, 0.617, 0.583, 0.569];
recall = [0.641, 0.678, 0.602, 0.629, 0.609];
precision = [0.764, 0.727, 0.702, 0.599, 0.618];
%%
all = vertcat(f_score, recall, precision);
%%
bar(all)
xticks([1, 2, 3])
xticklabels({'f-score', 'recall', 'precision'})
legend('mrcnn', 'STNeuroNet', 'HNCcorr + conv2d', 'Sourcery (Suite2p)', 'Unet2DS')




%% check saved hdf5 file
h5_filename = 'E:\data_of_different_brian_area\xinyu\20180910_b55a0706\ROIInfo\ROIInfo_4.h5';
masks = uint8(h5read(h5_filename, '/masks'));
trial_No = int16(h5read(h5_filename, '/trial_No'));
frame_No = int16(h5read(h5_filename, '/frame_No'));
%%
imshow(double(masks(:, :, 1)))



%%
tif_filename = 'E:\data_of_different_brian_area\xinyu\20180910_b55a0706\b55a07_test06_2x_2afc_160um_20180910_dftReg_001.tif';
% filename = 'E:\data_of_different_brian_area\zhonglin\h36_20161127\ZL_h36_20161127_fied1_d150_3x_dftReg_001.tif';
raw_data = load_raw_data(tif_filename);
%%
bin_size = 4;
% for i = 1:floor(size(raw_data, 3) / bin_size) - 2
i = 62;
temp_1 = norm_img(mean(raw_data(:, :, (i-1)*bin_size+1:i*bin_size), 3));
temp_2 = norm_img(mean(raw_data(:, :, (i)*bin_size+1:(i+1)*bin_size), 3));
temp_3 = norm_img(mean(raw_data(:, :, (i+1)*bin_size+1:(i+2)*bin_size), 3));
temp = cat(3, temp_1, temp_2, temp_3);
imshow(temp)
% end

%%
function img = norm_img(img)
	img = img - min(img(:));
	img = img / max(img(:));
end
function image_data = load_raw_data(tiff_filename)
	info = imfinfo(tiff_filename);
	numImages = length(info);
	tiff = Tiff(tiff_filename, 'r');
	temp = double(read(tiff));
	image_data = zeros(size(temp, 1), size(temp, 2), numImages);
	image_data(:, :, 1) = temp;
	for j = 2:numImages
		tiff.setDirectory(j);
		temp = double(read(tiff));
		image_data(:, :, j) = temp;
	end
end 





