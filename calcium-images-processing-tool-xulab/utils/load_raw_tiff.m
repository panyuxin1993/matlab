function [raw_images, mean_image, max_image] = load_raw_tiff(filename)
	info = imfinfo(filename);
	numImages = length(info);
	img_h = info(1).Height;
	img_w = info(1).Width;
	raw_images = zeros(img_h, img_w, numImages);
	tiff = Tiff(filename, 'r');
	for j = 1:numImages
		tiff.setDirectory(j);
		temp = double(read(tiff));
		raw_images(:, :, j) = temp;
	end
	raw_images = raw_images - min(raw_images(:));
	raw_images = raw_images / max(raw_images(:));
	
	mean_image = mean(raw_images, 3);
	max_image = max(raw_images, [], 3);
end