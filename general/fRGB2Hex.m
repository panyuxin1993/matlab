function [colorhex] = fRGB2Hex(colorRGB)
%FHEX2RGB change color value from hex to RGB(0-255)/255
%Input- colorhex, cell array of color value in hex(string)
%Output- colormat, cell array of 1-by-3 double, indicating the color
colorhex=cell(size(colorRGB));
for i=1:length(colorRGB)
    strRGB=[];
    for iRGB=1:3
        strRGB=strcat(strRGB,dec2hex(colorRGB{i}(iRGB)*255));
    end
    colorhex{i}=strRGB;
end
end

