function CaSignal = update_signal_show(handles, CaSignal)
	if CaSignal.CurrentROINo <= numel(CaSignal.ROIs)
		[y_start, y_end, x_start, x_end] = CaSignal.ROIs{CaSignal.CurrentROINo}{1: 4};
		BW = CaSignal.ROIs{CaSignal.CurrentROINo}{5};
		tempMask = zeros(CaSignal.image_height, CaSignal.image_width);
		tempMask(y_start:y_end, x_start:x_end) = BW(1:y_end - y_start + 1, 1:x_end - x_start + 1);
		tempMask = tempMask / sum(tempMask(:));
		signal = reshape(tempMask, [1, CaSignal.image_height * CaSignal.image_width]) * ...
			reshape(CaSignal.raw_images, [CaSignal.image_height * CaSignal.image_width, size(CaSignal.raw_images, 3)]);
		axes(handles.SignalShowAxes)
		plot(signal)
	else
		axes(handles.SignalShowAxes)
		plot(0)
	end
end