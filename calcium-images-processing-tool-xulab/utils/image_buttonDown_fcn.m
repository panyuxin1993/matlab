function CaSignal = Image_buttonDown_fcn(hObject,eventdata, handles, CaSignal)

	x = round(eventdata.IntersectionPoint(1));
	y = round(eventdata.IntersectionPoint(2));
	CaSignal.TempXY = [x, y];
	if ~isequal(CaSignal.SummarizedMask, []) && CaSignal.SummarizedMask(y, x) > 0
		disp('select ROI');
		CaSignal.CurrentROINo = CaSignal.SummarizedMask(y, x);
		set(handles.CurrentROINoEdit, 'String', num2str(CaSignal.CurrentROINo));
		CaSignal = update_subimage_show(handles, CaSignal, true);
		CaSignal = update_signal_show(handles, CaSignal);
	else
		disp('update sub_image');
		CaSignal.CurrentROINo = CaSignal.ROI_num + 1;
		CaSignal = update_subimage_show(handles, CaSignal, false);
	end
end