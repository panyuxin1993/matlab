close all;
path='D:\xulab\behavior\pyx292';
% [animal_name,dataChoice]=fFindChoice('Y:\Lab_Members\Duan\mice_performance\PYX179',20);
[animal_name,dataChoice]=fFindChoice(path,10);
fig=fPlotCorrectRate(animal_name,dataChoice,10);      
saveas(fig,[path,filesep, animal_name,'.png']);
% frame= getframe(fig);
% im= frame2im(frame);
% imwrite(im,[path,filesep, animal_name,'\performance of ',animal_name,'.png'],'png'); %alternative, saveas(); 