%example
%{
videofile='F:\2P\example\pyx290_20200528_920nm_power50_3X_dftReg_015-syn-contra.avi';
% videofile='H:\2P\pyx290_20200528\im_data_reg\pyx290_20200528_920nm_power50_3X_dftReg_015.tif';
behfile='H:\2P\pyx290_20200528\im_data_reg\result_save\pyx290_20200528-imaging.mat';
load('H:\2P\pyx290_20200528\im_data_reg\result_save\CaTrialsSIM_pyx290_20200528_920nm_power50_3X_dftReg_.mat');
framerate=1000/SavedCaTrials.FrameTime;
trialInd=15;
fLabel2PvideoBeh(videofile,behfile,trialInd,framerate)
%}

%for 2020 annual meeting
%
%videofile='F:\2P\example\pyx290_20200528_920nm_power50_3X_dftReg_015-syn-contra.avi';
videofile='H:\2P\pyx285_20200520\im_data_reg\pyx285_20200520_920nm_power50_3x_dftReg_049.tif';
behfile='H:\2P\pyx285_20200520\im_data_reg\result_save\pyx285_20200520_imaging.mat';
load('H:\2P\pyx285_20200520\im_data_reg\result_save\CaTrialsSIM_pyx285_20200520_920nm_power50_3x_dftReg_.mat');
framerate=1000/SavedCaTrials.FrameTime;
trialInd=49;
fLabel2PvideoBeh(videofile,behfile,trialInd,framerate,0.5)
%}
