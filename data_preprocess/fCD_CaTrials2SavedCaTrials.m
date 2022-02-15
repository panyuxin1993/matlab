function [SavedCaTrials] = fCD_CaTrials2SavedCaTrials(folder)
%fCD_CaTrials2SavedCaTrials change data form of imaging extracted calcium
%data from CaTrials (struct array) back to SavedCaTrials (struct) 
%   folder is where CaTrials saved
filedir=dir(folder);
dircell=struct2cell(filedir);
filenames=dircell(1,:);
indCaTrials=cellfun(@(x) contains(x,'CaTrials_'), filenames);
disp(['processing:',folder,filesep,filenames{indCaTrials}]);
if sum(indCaTrials)>0
    load([folder,filesep,filenames{indCaTrials}]);
else
    disp(['no file and return:',folder,filesep,filenames{indCaTrials}]);
    return;
end
indSavedCaTrials=cellfun(@(x) contains(x,'CaTrialsSIM_'), filenames);
load([folder,filesep,filenames{indSavedCaTrials}]);
SavedCaTrials.FileName_prefix= CaTrials(1).FileName_prefix;
SavedCaTrials.nFrames = CaTrials(1).nFrames;
SavedCaTrials.FrameTime = CaTrials(1).FrameTime;
SavedCaTrials.nROIs= CaTrials(1).nROIs;
nTrials=length(CaTrials);
SavedCaTrials.TrialNum=nTrials;
SavedCaTrials.f_raw=cell(nTrials,1);
SavedCaTrials.RingF=cell(nTrials,1);
for i=1:nTrials
    SavedCaTrials.f_raw{i}=CaTrials(i).f_raw;
    SavedCaTrials.RingF{i}=CaTrials(i).RingF;
end
save([folder,filesep,filenames{indSavedCaTrials}],'SavedCaTrials');
end

