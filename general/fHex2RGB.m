function [colorRGB]=fHex2RGB(colorhex)
%FHEX2RGB change color value from hex to RGB(0-255)/255
%Input- colorhex, cell array of color value in hex(string)
%Output- colormat, cell array of 1-by-3 double, indicating the color
colorRGB=cell(size(colorhex));
for i=1:length(colorhex)
    rgb=zeros(1,3);
    for iRGB=1:3
        strRGB=colorhex{i}(2*iRGB-1:2*iRGB);
        rgb(iRGB)=(hex2dec(strRGB))/255;
    end
    colorRGB{i}=rgb;
end
end

