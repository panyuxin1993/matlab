function [rotate_x,rotate_y] = fRotateCurve(x,y,angle,rx0,ry0)
%FROTATECURVE ref https://www.pianshen.com/article/758864846/
%rotate a given curve to specified angle
%   Detailed explanation goes here
% x=-10:0.001:10;
% y=sin(x);
% plot(x,y);hold on;
% middle=fix(length(x)/2);
% rx0=x(middle);
% ry0=y(middle);
% angle=15;
    
angle=pi/180*angle;
for ii=1:length(y)
    x0=(x(ii)-rx0)*cos(angle)-(y(ii)-ry0)*sin(angle)+rx0;
    y0=(x(ii)-rx0)*sin(angle)+(y(ii)-ry0)*cos(angle)+ry0;
    rotate_x(ii)=x0;
    rotate_y(ii)=y0;
end
% plot(rotate_x,rotate_y);
end

