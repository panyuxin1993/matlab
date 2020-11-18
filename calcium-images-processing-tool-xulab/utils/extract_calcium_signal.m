function extract_calcium_signal(CaSignal)

if CaSignal.ROI_num > 0
	% convert local mask to while image mask
	masks = zeros(CaSignal.image_height, CaSignal.image_width, CaSignal.ROI_num);
	for i = 1:CaSignal.ROI_num 
		y_start = CaSignal.ROIs{i}{1};
		y_end = CaSignal.ROIs{i}{2};
		x_start = CaSignal.ROIs{i}{3};
		x_end = CaSignal.ROIs{i}{4};
		tempRoi = CaSignal.ROIs{i}{5};
		masks(y_start:y_end, x_start:x_end, i) = tempRoi(1:y_end - y_start + 1, 1:x_end - x_start + 1);
	end
	temp_mask = reshape(masks, [size(masks, 1)*size(masks, 2), size(masks, 3)]);
	% load image data
	cal_signal = {};
	parfor i = 1:numel(CaSignal.imageFilenames)
        fprintf('%d/%d\n', [i, numel(CaSignal.imageFilenames)])
		tif_data = [];
		tiff_filename = CaSignal.imageFilenames{i};
		info = imfinfo(tiff_filename);
		numImages = length(info);
		tiff = Tiff(tiff_filename, 'r');
		for k = 1:numImages
			tiff.setDirectory(k);
			tif_data(:, :, k) = double(read(tiff));
		end
		tif_data = reshape(tif_data, [size(tif_data, 1)*size(tif_data, 2), size(tif_data, 3)]);
		cal_signal{i} = temp_mask'*tif_data;
	end
	
	SavedCaTrials.f_raw = cal_signal;
	SavedCaTrials.nROIs = CaSignal.ROI_num;
	SavedCaTrials.TrialNum = numel(CaSignal.imageFilenames);
	% get frame time
	Key = 'scanFrameRate = ';
	info = imfinfo(CaSignal.imageFilenames{1});
	Str = info(1).ImageDescription;
	Index = strfind(Str, Key);
	startIndex = Index + length(Key);
	endIndex = strfind(Str(startIndex:end), newline);
	scanFrameRate = str2double(Str(startIndex: startIndex + endIndex(1) -1));
	SavedCaTrials.FrameTime = 1000 / scanFrameRate;
	
	path = fullfile(CaSignal.imagePathName, 'result_save');
	% make dir if not exist
	if ~exist(path, 'dir')
		[status, msg, msgID] = mkdir(path);
		if status == 0
			error(msgID, msg);
		end
	end
	
	% automated generate saved filename based on tiff filename
	if numel(CaSignal.imageFilenames) > 0
		[~, tif_filename, ~] = fileparts(CaSignal.imageFilenames{1});
		temp = split(tif_filename, '_');
		temp = [{'CaTrialsSIM'}, temp(1:end-1)', {'.mat'}];
		mat_filename = strjoin(temp, '_');
	else
		error('No Tif images founded!')
	end
	
	save(fullfile(path, mat_filename), 'SavedCaTrials')
	
end
end