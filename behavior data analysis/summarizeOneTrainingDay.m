%summarize one training day in all aspects, including performance across
%recent day, performance along delay length, licking raster, and licking
%consistency index, etc. 
%difficult to plot in one figure, so plot separately and arrange in good
%way
close all;
clear;
%specify the path and file name, etc. 
path='H:\xulab\behavior';
ndelaygroup=4;
date='2020_02_03';
animal='pyx355';
fileFolder=[path, filesep,animal];
npast=10;%how many sessions performance will be plotted
%% plot performance in recent trainning days
[animal_name,dataChoice]=fFindChoice(fileFolder,npast);
figPerformance=fPlotCorrectRate(animal_name,dataChoice,npast);      
saveas(figPerformance,[fileFolder,filesep, animal_name,'.png']);
set(figPerformance, 'position', [0 50 900 300]);%¿ØÖÆfig³ß´ç
%% plot licking consistency index
figLCI = fPlotLickConsistency(path,date,animal,npast);
set(figLCI, 'position', [900 50 600 300]);%¿ØÖÆfig³ß´ç
%% plot different delay performance
figDiffDelay = fPlotDiffDelay(ndelaygroup,path,date,animal);
set(figDiffDelay, 'position', [0 400 1200 300]);%¿ØÖÆfig³ß´ç
