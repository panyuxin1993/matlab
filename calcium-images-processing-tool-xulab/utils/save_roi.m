function save_roi(CaSignal)

% filter = {'*.h5; *.hdf5; *.mat', ...
% 	'ROIInfo Files (*.h5, *.hdf5, *.mat)'};
% [file, path] = uiputfile(filter, 'Save ROI information', CaSignal.imagePathName);
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
	temp = [{'ROIinfoBU'}, temp(1:end-1)', {'.mat'}];
	mat_filename = strjoin(temp, '_');
else
	error('No Tif images founded!')
end

% save ROI information
masks = zeros(CaSignal.image_height, CaSignal.image_width, CaSignal.ROI_num);
trial_No = zeros(CaSignal.ROI_num, 1);

ROImasksAll = cell(CaSignal.ROI_num,1);
ROIposAll = cell(CaSignal.ROI_num,1);
ROItypeAll = cell(CaSignal.ROI_num,1);
for i = 1:CaSignal.ROI_num 
	y_start = CaSignal.ROIs{i}{1};
	y_end = CaSignal.ROIs{i}{2};
	x_start = CaSignal.ROIs{i}{3};
	x_end = CaSignal.ROIs{i}{4};
	tempRoi = CaSignal.ROIs{i}{5};
	masks(y_start:y_end, x_start:x_end, i) = tempRoi(1:y_end - y_start + 1, 1:x_end - x_start + 1);
	trial_No(i) = CaSignal.ROIs{i}{7};

	ROImasksAll{i} = masks(:, :, i);
	B=bwboundaries(ROImasksAll{i});
	ROIposAll{i} = B{1}(:,[2,1]);
	ROItypeAll{i} = 'Soma';
end

h5_filename = strrep(mat_filename, 'mat', 'h5');

if exist(fullfile(path, h5_filename), 'file') == 2
	delete(fullfile(path, h5_filename))
end

disp('Save ROIinfo file');
h5create(fullfile(path, h5_filename), '/masks', size(masks))
h5write(fullfile(path, h5_filename), '/masks', uint8(masks));
h5create(fullfile(path, h5_filename), '/trial_No', size(trial_No))
h5write(fullfile(path, h5_filename), '/trial_No', int16(trial_No));

ROIinfoBU.ROImask = ROImasksAll;
ROIinfoBU.ROIpos = ROIposAll;
ROIinfoBU.ROItype = ROItypeAll;
ROIinfoBU.ROI_def_trialNo = trial_No;
ROIinfoBU.BGpos = [];
ROIinfoBU.BGmask = {};
ROIinfoBU.Method = '';
save(fullfile(path, mat_filename), 'ROIinfoBU');

end

