function [CaSignal, ROIs] = load_roi(filename, CaSignal)
	ROIs = {};
	[ROImasks, trial_No] = load_ROImasks(filename);
	for i = 1:length(ROImasks)
		tempROI = {};
		mask = ROImasks{i};
		[ys, xs] = find(mask);
		x = uint16(median(xs, 'all'));
		y = uint16(median(ys, 'all'));
		[x_start, x_end, y_start, y_end] = generate_loaction_boxes(CaSignal, x, y);
		tempROI{1} = y_start;
		tempROI{2} = y_end;
		tempROI{3} = x_start;
		tempROI{4} = x_end;
		sub_mask = uint8(zeros(2 * CaSignal.ROIDiameter + 1,  2 * CaSignal.ROIDiameter + 1));
		sub_mask(1:y_end - y_start + 1, 1:x_end - x_start + 1) = mask(y_start:y_end, x_start:x_end);
		tempROI{5} = sub_mask;
		B = bwboundaries(sub_mask, 'noholes');
		tempROI{6} = B{1};
		tempROI{7} = trial_No(i);
		ROIs{i} = tempROI;
	end
	CaSignal = generate_summarizedMask(CaSignal);
end

function [ROImasks, trial_No] = load_ROImasks(ROI_file)
	ROImasks = {};
	[~, ~, ext] = fileparts(ROI_file);
	if strcmp(ext, '.mat')
		ROIinfo = load(ROI_file);
		if isfield(ROIinfo, 'ROIinfoBU')
			ROImasks = ROIinfo.ROIinfoBU.ROImask;
			trial_No = ROIinfo.ROIinfoBU.ROI_def_trialNo;
		else
			errordlg(['Not find any ROIinfo field in ', ROI_file], 'File Error');
			return;
		end
	elseif strcmp(ext, '.hdf5')
		masks = h5read(ROI_file, '/masks');
		for i = 1:size(masks, 3)
			ROImasks{i} = masks(:, :, i);
		end
		trial_No = int16(h5read(ROI_file, '/trial_No'));
	end
end


function CaSignal = generate_summarizedMask(CaSignal)
% 	ROIs = CaSignal.ROIs;
% 	assignin('base','ROIs',ROIs)
	CaSignal.SummarizedMask = zeros(CaSignal.image_height, CaSignal.image_width);
	for i = 1:CaSignal.ROI_num
		tempMask =  zeros(CaSignal.image_height, CaSignal.image_width);
		y_start = CaSignal.ROIs{i}{1};
		y_end = CaSignal.ROIs{i}{2};
		x_start = CaSignal.ROIs{i}{3};
		x_end = CaSignal.ROIs{i}{4};
		tempRoi = CaSignal.ROIs{i}{5};
		tempMask(y_start:y_end, x_start:x_end) = tempRoi(1:y_end - y_start + 1, 1:x_end - x_start + 1);
		CaSignal.SummarizedMask(tempMask == 1) = i;
	end
end