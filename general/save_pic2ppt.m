%ref https://itectec.com/matlab/matlab-how-to-add-multiple-png-files-to-an-existing-powerpoint-pptx-file-using-a-for-loop/
%removed the clear all, close all, clc, as they are irrelevant here and can only cause problems
project_dir = pwd();
%new% List all current folder contents ending with .png. Resulting names will
% appear in the order returned by the operating system.
files = dir( fullfile(project_dir, '*.png'));
%changed% Create Common Object Model (COM) server so MATLAB can export data to
% PowerPoint
g = actxserver('powerpoint.application');
% Open PowerPoint and make it visible
g.Visible = 1;
Presentation = g.Presentation;
% Prompt the user for the PowerPoint file to amend
[fn, pn] = uigetfile('*.pptx', 'Select PowerPoint File To Amend');
filename = fullfile(pn, fn);
%changedPresentation = invoke(Presentation, 'open', filename);
% Get current number of slidesslide_count = get(Presentation.Slides, 'Count');
% Export all PNGs in the current directory to the PowerPoint file specified% above. The following slides will be added to the END of the PowerPoint
% file. All slides will have a common title.
for i=1:length(files)
    slide_count = int32(double(slide_count)+1);
    slide = invoke(Presentation.Slides, 'Add', slide_count{i}, 11);
    set(slide.Shapes.Title.Textframe.Textrange, 'Text', 'SomeTitle');
    slidefile = fullfile(project_dir, files(i).name);
    %new
    Image{i} = slide.Shapes.AddPicture(slidefile, 'msoFalse', 'msoTrue', 0, 80, 720, 440);
    %changed
end
% Save the amended PowerPoint presentation to the current directory
outfile = fullfile(project_dir, 'DRAFT.pptx');
%new
Presentation.SaveAs(outfile);
%changed% Close PowerPoint as a COM automation
serverg.Quit;
g.delete;