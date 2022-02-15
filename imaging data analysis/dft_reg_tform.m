function [im_dft] = dft_reg_tform(im_s, im_tg, pix_range)
% using imregdemons function to do image registration
% im_s, source image data, can be multi-frame, m x n x p matrix.
% im_tg, target image for registration, a single frame image, m x n matrix.
%        If left empty, will generate default target image.
% 
% pix_range, specify the fraction of image used for registraiton.
%            If with 2 elements [startRow endRow], take as line_range,
%            If with 4 elements [startRow endRow startCol endCol], use a rectangle fraction of
%            the image.
% - NX
% Updated Aug, 2012, - NX

if nargin < 2 || isempty(im_tg)
    % By default, use mean of the last 10 frames as the target image.
    im_tg = mean(im_s(:,:,end-9:end),3);
end
if nargin < 3 || isempty(pix_range)
    row_nums = 1 : size(im_tg,1);
    col_nums = 1: size(im_tg, 2);
else
    if numel(pix_range) == 2
        row_nums = pix_range(1) : pix_range(2);
        col_nums = 1:size(im_tg, 2);
    end
    if numel(pix_range) == 4
        row_nums = pix_range(1):pix_range(2);
        col_nums = pix_range(3):pix_range(4);
    end
end

fixed=im_tg;
class_str = class(im_s);
im_dft = ones(size(im_s) + [0, 0, 0], class_str);
im_dft_gpu=gpuArray(im_dft);
for i=1:size(im_s,3)
    moving=im_s(row_nums, col_nums,i);
%     imshowpair(fixed,moving,'montage')
    %imshowpair(fixed,moving)
%     moving = imhistmatch(moving,fixed); int16 vs double
%     [~,movingReg] = imregdemons(moving,fixed,[100 80 50],'DisplayWaitbar',false);
    %use GPU to analyze
    moving_gpu=gpuArray(moving);
    fixed_gpu=gpuArray(fixed);
%     [~,movingReg_gpu] = imregdemons(moving_gpu,fixed_gpu,'DisplayWaitbar',false);
    [~,movingReg_gpu] = imregdemons(moving_gpu,fixed_gpu);
    im_dft_gpu(:,:,i)=movingReg_gpu;
%     movingReg=gather(movingReg_gpu);
%     %imshowpair(fixed,movingReg)
%     im_dft(:,:,i)=movingReg;
end
im_dft=gather(im_dft_gpu);