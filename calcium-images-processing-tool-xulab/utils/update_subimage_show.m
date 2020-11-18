function CaSignal = update_subimage_show(handles, CaSignal, with_ROI)
	
	if CaSignal.CurrentROINo <= CaSignal.ROI_num && with_ROI
		axes(handles.SubimageShowAxes);
		[y_start, y_end, x_start, x_end] = CaSignal.ROIs{CaSignal.CurrentROINo}{1:4};
		img = CaSignal.showing_image(y_start:y_end, x_start:x_end, :);
		p_bottom = prctile(img, CaSignal.bottom_percentile, 'all');
		p_top = prctile(img, CaSignal.top_percentile, 'all');
		img = imadjust(img, [p_bottom, p_top]);
		CaSignal.h_subimage = imshow(img);
		hold on; 
		if get(handles.ShowROINoCheckbox, 'Value') == 0    
			boundary = CaSignal.ROIs{CaSignal.CurrentROINo}{6};
			plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 1)
		else	    
			boundary = CaSignal.ROIs{CaSignal.CurrentROINo}{6};
			plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 1)
			text(min(boundary(:,2)), min(boundary(:,1)), num2str(CaSignal.CurrentROINo), 'Color', 'r', 'FontSize', 10);	
		end
		hold off;
	else
		x = CaSignal.TempXY(1);
		y = CaSignal.TempXY(2);
		axes(handles.SubimageShowAxes);
		[x_start, x_end, y_start, y_end] = generate_loaction_boxes(CaSignal, x, y);
		img = CaSignal.showing_image(y_start:y_end, x_start:x_end, :);
		p_bottom = prctile(img, CaSignal.bottom_percentile, 'all');
		p_top = prctile(img, CaSignal.top_percentile, 'all');
		img = imadjust(img, [p_bottom, p_top]);
		CaSignal.h_subimage = imshow(img);
	end
end