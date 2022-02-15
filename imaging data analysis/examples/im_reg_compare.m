fixed=imread('C:\Data\pyx397_20211102\im_data_reg\pyx397_20211102_deep20um_920nm_power30_12x_00dftReg_002.tif',181);
moving=imread('C:\Data\pyx397_20211102\im_data_reg\pyx397_20211102_deep20um_920nm_power30_12x_00dftReg_002.tif',184);
imshowpair(fixed.*100, moving.*100,'Scaling','joint');
fixed=imread('C:\Data\pyx397_20211102\im_data_reg\im_data_reg\pyx397_20211102_deep20um_920nm_power30_12x_00dftReg_dftRegDemons_002.tif',115);
moving=imread('C:\Data\pyx397_20211102\im_data_reg\im_data_reg\pyx397_20211102_deep20um_920nm_power30_12x_00dftReg_dftRegDemons_002.tif',116);
imshowpair(fixed.*100, moving.*100,'Scaling','joint');
%next try imregister