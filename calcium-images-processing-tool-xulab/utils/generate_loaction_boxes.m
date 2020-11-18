 function [x_start, x_end, y_start, y_end] = generate_loaction_boxes(CaSignal, x, y)
	x_start = x - CaSignal.ROIDiameter;
	if x_start < 1 
		x_start = 1;
	end
	x_end = x_start + CaSignal.ROIDiameter * 2;
	if x_end > CaSignal.image_width
		x_end = CaSignal.image_width;
	end
	y_start = y - CaSignal.ROIDiameter;
	if y_start < 1
		y_start = 1;
	end
	y_end = y_start + CaSignal.ROIDiameter * 2;
	if y_end > CaSignal.image_height
		y_end = CaSignal.image_height;
	end
end