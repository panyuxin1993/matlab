function CaSignal = delete_all_roi(CaSignal)
	CaSignal.ROIs = {};
	CaSignal.ROI_num = 0;
	CaSignal.SummarizedMask = zeros(CaSignal.image_height, CaSignal.image_width);
end