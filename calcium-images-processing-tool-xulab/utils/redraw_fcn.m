function CaSignal = redraw_fcn(handles, CaSignal)
	x = CaSignal.TempXY(1);
	y = CaSignal.TempXY(2);
	CaSignal = update_subimage_show(handles, CaSignal, false);
	[x_start, x_end, y_start, y_end] = generate_loaction_boxes(CaSignal, x, y);
	h_draw = drawfreehand;
	if numel(h_draw) == 0
		return;
	end
	% pos = h_draw.getPosition;
	BW = createMask(h_draw);
	B = bwboundaries(BW, 'noholes');
	boundary = B{1};
	if numel(B) > 0
		CaSignal.ROIs{CaSignal.CurrentROINo} = {y_start, y_end, x_start, x_end, BW, boundary, CaSignal.current_trial};
		tempMask = zeros(CaSignal.image_height, CaSignal.image_width);
		tempMask(y_start:y_end, x_start:x_end) = BW(1:y_end - y_start + 1, 1:x_end - x_start + 1);
		CaSignal.SummarizedMask(CaSignal.SummarizedMask == CaSignal.CurrentROINo) = 0;
		CaSignal.SummarizedMask(tempMask > 0) = CaSignal.CurrentROINo;
		CaSignal.ROI_num = numel(CaSignal.ROIs);
	end
	CaSignal = update_signal_show(handles, CaSignal);
	CaSignal = update_subimage_show(handles, CaSignal, true);
	set(handles.CurrentROINoEdit, 'String', num2str(CaSignal.CurrentROINo));
	set(handles.ROINumShowText, 'String', num2str(CaSignal.ROI_num));
end
