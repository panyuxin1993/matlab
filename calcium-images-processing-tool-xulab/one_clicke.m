function varargout = one_clicke(varargin)
% ONE_CLICKE MATLAB code for one_clicke.fig
%      ONE_CLICKE, by itself, creates a new ONE_CLICKE or raises the existing
%      singleton*.
%
%      H = ONE_CLICKE returns the handle to a new ONE_CLICKE or the handle to
%      the existing singleton*.
% 
%      ONE_CLICKE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ONE_CLICKE.M with the given input arguments.
%
%      ONE_CLICKE('Property','Value',...) creates a new ONE_CLICKE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before one_clicke_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to one_clicke_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help one_clicke

% Last Modified by GUIDE v2.5 30-Jul-2019 14:34:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @one_clicke_OpeningFcn, ...
                   'gui_OutputFcn',  @one_clicke_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

%**************************************************************************%
function Image_ButtonDownFcn(hObject,eventdata, handles)
global CaSignal
CaSignal = image_buttonDown_fcn(hObject,eventdata, handles, CaSignal);
CaSignal = Update_Image_Fcn(handles, CaSignal, true);

function CaSignal = Update_Image_Fcn(handles, CaSignal, restore_zoom)
CaSignal = update_image_show(handles, CaSignal, restore_zoom);
CaSignal.h_image.ButtonDownFcn = {@Image_ButtonDownFcn, handles};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes just before one_clicke is made visible.
function one_clicke_OpeningFcn(hObject, eventdata, handles, varargin)
global CaSignal
handles.output = hObject;
% initialize CaSignal
% about ROI
CaSignal.ROIs = {};
CaSignal.ROI_num = 0;
CaSignal.CurrentROINo = 1;
CaSignal.ROIDiameter = 12;

% about showing image
CaSignal.imagePathName = pwd;
CaSignal.image_width = 0;
CaSignal.image_height = 0;
CaSignal.imageFilenames = '';
CaSignal.imagePathName = '';
CaSignal.raw_images = [];
CaSignal.mean_images = [];
CaSignal.showing_image = [];
CaSignal.current_trial = 0;
CaSignal.total_trial = 0;
CaSignal.current_frame = 0;
CaSignal.NumOfFrames = 0;
CaSignal.top_percentile = 100.0;
CaSignal.bottom_percentile = 0.0;
% about showing sub_image
CaSignal.TempXY = [64, 64];
% CaSignal.RedrawBasedOnTempROI = true;
CaSignal.SummarizedMask = [];

cd(fileparts(which(mfilename)));
addpath('./');
addpath('./utils');
addpath('./ssh2_v2_m1_v7');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes one_clicke wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = one_clicke_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on button press in ChoseFileButton.
function ChoseFileButton_Callback(hObject, eventdata, handles)


function DataPathEdit_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function DataPathEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function ImageShowAxes_ButtonDownFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function SubimageShowAxes_CreateFcn(hObject, eventdata, handles)

function ImageShowAxes_CreateFcn(hObject, eventdata, handles)




% --- Executes on button press in ReDrawButton.
function ReDrawButton_Callback(hObject, eventdata, handles)
global CaSignal
CaSignal = redraw_fcn(handles, CaSignal);
CaSignal = Update_Image_Fcn(handles, CaSignal, true);



% --- Executes on button press in DeleteButton.
function DeleteButton_Callback(hObject, eventdata, handles)
global CaSignal
if CaSignal.CurrentROINo <= CaSignal.ROI_num
	fprintf('delete ROI %d\n\r', CaSignal.CurrentROINo);
	CaSignal.SummarizedMask(CaSignal.SummarizedMask(:, :) == CaSignal.CurrentROINo) = 0;
	CaSignal.ROIs(CaSignal.CurrentROINo) = [];
	CaSignal.ROI_num = numel(CaSignal.ROIs);
	CaSignal.SummarizedMask(CaSignal.SummarizedMask > CaSignal.CurrentROINo) = ...
	CaSignal.SummarizedMask(CaSignal.SummarizedMask > CaSignal.CurrentROINo) - 1;
	CaSignal = Update_Image_Fcn(handles, CaSignal, true);
	CaSignal = update_subimage_show(handles, CaSignal, true);
	CaSignal = update_signal_show(handles, CaSignal);
	set(handles.ROINumShowText, 'String', num2str(CaSignal.ROI_num));
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
global CaSignal
CaSignal.top_percentile = get(hObject, 'Value');
disp(CaSignal.top_percentile)
if CaSignal.top_percentile <= CaSignal.bottom_percentile
	CaSignal.top_percentile = CaSignal.bottom_percentile + 0.1;
end
set(hObject, 'Value', CaSignal.top_percentile);
CaSignal = Update_Image_Fcn(handles, CaSignal, true);
CaSignal = update_subimage_show(handles, CaSignal, true);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
global CaSignal
CaSignal.bottom_percentile = get(hObject, 'Value');
if CaSignal.top_percentile <= CaSignal.bottom_percentile
	CaSignal.bottom_percentile = CaSignal.top_percentile - 0.1;
end
set(hObject, 'Value', CaSignal.bottom_percentile);

CaSignal = Update_Image_Fcn(handles, CaSignal, true);
CaSignal = update_subimage_show(handles, CaSignal, true);

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function TrialNumEdit_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function TrialNumEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NextTrialButton.
function NextTrialButton_Callback(hObject, eventdata, handles)
global CaSignal
CaSignal.current_trial = CaSignal.current_trial + 1;
if CaSignal.current_trial > CaSignal.total_trial
	CaSignal.current_trial = 1;
end
[CaSignal.raw_images, CaSignal.mean_images, CaSignal.max_images] = load_raw_tiff(CaSignal.imageFilenames{CaSignal.current_trial});
set(handles.TrialNumEdit, 'String', sprintf('%d/%d', CaSignal.current_trial, CaSignal.total_trial));

CaSignal.current_frame = 1;
CaSignal.NumOfFrames = size(CaSignal.raw_images, 3);
CaSignal.showing_image = CaSignal.raw_images(:, :, CaSignal.current_frame);
if get(handles.checkbox_mean, 'Value') == 1
	CaSignal.showing_image = CaSignal.mean_images;
elseif get(handles.checkbox_max, 'Value') == 1
	CaSignal.showing_image = CaSignal.max_images;
end
set(handles.FrameNoText, 'String', sprintf('%d/%d', CaSignal.current_frame, CaSignal.NumOfFrames));

CaSignal = Update_Image_Fcn(handles, CaSignal, false);
CaSignal = update_signal_show(handles, CaSignal);
CaSignal = update_subimage_show(handles, CaSignal, true);


% --- Executes on button press in PreviousTrialButton.
function PreviousTrialButton_Callback(hObject, eventdata, handles)
global CaSignal
CaSignal.current_trial = CaSignal.current_trial - 1;
if CaSignal.current_trial < 1
	CaSignal.current_trial = CaSignal.total_trial;
end
[CaSignal.raw_images, CaSignal.mean_images, CaSignal.max_images] = load_raw_tiff(CaSignal.imageFilenames{CaSignal.current_trial});
set(handles.TrialNumEdit, 'String', sprintf('%d/%d', CaSignal.current_trial, CaSignal.total_trial));

CaSignal.current_frame = 1;
CaSignal.NumOfFrames = size(CaSignal.raw_images, 3);
CaSignal.showing_image = CaSignal.raw_images(:, :, CaSignal.current_frame);
if get(handles.checkbox_mean, 'Value') == 1
	CaSignal.showing_image = CaSignal.mean_images;
elseif get(handles.checkbox_max, 'Value') == 1
	CaSignal.showing_image = CaSignal.max_images;
end
set(handles.FrameNoText, 'String', sprintf('%d/%d', CaSignal.current_frame, CaSignal.NumOfFrames));


CaSignal = Update_Image_Fcn(handles, CaSignal, false);
CaSignal = update_signal_show(handles, CaSignal);
CaSignal = update_subimage_show(handles, CaSignal, true);


function GoToTrialNoEdit_Callback(hObject, eventdata, handles)
global CaSignal
CaSignal.current_trial = int16(str2num(get(handles.GoToTrialNoEdit, 'String')));
if CaSignal.current_trial > CaSignal.total_trial
	CaSignal.current_trial = CaSignal.total_trial;
elseif CaSignal.current_trial < 1
	CaSignal.current_trial = 1;
end
[CaSignal.raw_images, CaSignal.mean_images, CaSignal.max_images] = load_raw_tiff(CaSignal.imageFilenames{CaSignal.current_trial});
set(handles.TrialNumEdit, 'String', sprintf('%d/%d', CaSignal.current_trial, CaSignal.total_trial));

CaSignal.current_frame = 1;
CaSignal.NumOfFrames = size(CaSignal.raw_images, 3);
CaSignal.showing_image = CaSignal.raw_images(:, :, CaSignal.current_frame);
if get(handles.checkbox_mean, 'Value') == 1
	CaSignal.showing_image = CaSignal.mean_images;
elseif get(handles.checkbox_max, 'Value') == 1
	CaSignal.showing_image = CaSignal.max_images;
end
set(handles.FrameNoText, 'String', sprintf('%d/%d', CaSignal.current_frame, CaSignal.NumOfFrames));

CaSignal = Update_Image_Fcn(handles, CaSignal, false);
CaSignal = update_signal_show(handles, CaSignal);
CaSignal = update_subimage_show(handles, CaSignal, true);


% --- Executes during object creation, after setting all properties.
function GoToTrialNoEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function SaveROIInfoTool_ClickedCallback(hObject, eventdata, handles)
global CaSignal
save_roi(CaSignal);


% --------------------------------------------------------------------
function OpenFileTool_ClickedCallback(hObject, eventdata, handles)
global CaSignal
[filename, pathName] = uigetfile(fullfile(CaSignal.imagePathName, '*.tif*'), 'Load Image File');
if isequal(filename,0)
	return;
end
CaSignal.imagePathName = pathName;
d = dir(fullfile(CaSignal.imagePathName, '*.tif*'));
CaSignal.imageFilenames = {};
for i = 1:size(d, 1)
	CaSignal.imageFilenames{i} = fullfile(CaSignal.imagePathName, d(i).name);
end
disp('Loading Image Data')
[CaSignal.raw_images, CaSignal.mean_images, CaSignal.max_images] = load_raw_tiff(CaSignal.imageFilenames{1});
disp('Done')
CaSignal.image_height = size(CaSignal.raw_images, 1);
CaSignal.image_width = size(CaSignal.raw_images, 2);
CaSignal.current_frame = 1;
CaSignal.NumOfFrames = size(CaSignal.raw_images, 3);
set(handles.FrameNoText, 'String', sprintf('%d/%d', CaSignal.current_frame, CaSignal.NumOfFrames));
CaSignal.showing_image = CaSignal.raw_images(:, :, CaSignal.current_frame);
CaSignal.current_trial = 1;
CaSignal.total_trial = numel(CaSignal.imageFilenames);
set(handles.TrialNumEdit, 'String', sprintf('%d/%d', CaSignal.current_trial, CaSignal.total_trial));
CaSignal.SummarizedMask = zeros(CaSignal.image_height, CaSignal.image_width);
CaSignal = Update_Image_Fcn(handles, CaSignal, false);
CaSignal = update_subimage_show(handles, CaSignal, true);



% --- Executes on button press in LoadROIButton.
function LoadROIButton_Callback(hObject, eventdata, handles)
global CaSignal
[filename, pathname] = uigetfile(fullfile(CaSignal.imagePathName, '*.*'), 'Load ROI File');
if isequal(filename,0)
	return;
end
[CaSignal, ROIs] = load_roi(fullfile(pathname, filename), CaSignal);
CaSignal.ROIs = ROIs;
CaSignal.ROI_num = size(CaSignal.ROIs, 2);
if CaSignal.ROI_num > 0
	CaSignal.CurrentROINo = 1;
	set(handles.CurrentROINoEdit, 'String', '1');
	set(handles.ROINumShowText, 'String', num2str(CaSignal.ROI_num));
	CaSignal = update_subimage_show(handles, CaSignal, true);
	CaSignal = update_signal_show(handles, CaSignal);
end
CaSignal = Update_Image_Fcn(handles, CaSignal, true);


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% global CaSignal
% if strcmp(eventdata.Key, 'd')
% 	DeleteButton_Callback(hObject, eventdata, handles);
% elseif strcmp(eventdata.Key, 'space')
% 	ReDrawButton_Callback(hObject, eventdata, handles);
% elseif strcmp(eventdata.Key, 'c')
% 	NextFrame_Callback(hObject, eventdata, handles);
% elseif strcmp(eventdata.Key, 'z')
% 	PreviousFrame_Callback(hObject, eventdata, handles);
% else
% 	return
% end


function figure1_KeyPressFcn(hObject, eventdata, handles)
global CaSignal
if strcmp(eventdata.Key, 'd')
	DeleteButton_Callback(hObject, eventdata, handles);
elseif strcmp(eventdata.Key, 'space')
	ReDrawButton_Callback(hObject, eventdata, handles);
elseif strcmp(eventdata.Key, 'c')
	NextFrame_Callback(hObject, eventdata, handles);
elseif strcmp(eventdata.Key, 'z')
	PreviousFrame_Callback(hObject, eventdata, handles);
else
	return
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
delete(hObject);
answer = questdlg('Save ROI information ?', 'Save query');
switch answer
    case 'Yes'
		SaveROIInfoTool_ClickedCallback(hObject, eventdata, handles)
        delete(hObject);
    case 'No'
       delete(hObject);
    case 'Cancel'
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)


% --- Executes on button press in DeleteAllButton.
function DeleteAllButton_Callback(hObject, eventdata, handles)
global CaSignal
answer = questdlg('Do you want to delete all ROI information ?', 'Alert query');
switch answer
    case 'Yes'
		CaSignal = delete_all_roi(CaSignal);
		set(handles.CurrentROINoEdit, 'String', '0');
		set(handles.ROINumShowText, 'String', '0');
		CaSignal = Update_Image_Fcn(handles, CaSignal, true);
		CaSignal = update_subimage_show(handles, CaSignal, true);
		CaSignal = update_signal_show(handles, CaSignal);
    case 'No'
		return
    case 'Cancel'
		return
end


function CurrentROINoEdit_Callback(hObject, eventdata, handles)
global CaSignal
current_roi_no = str2double(get(hObject,'String'));
current_roi_no = floor(current_roi_no);
if numel(CaSignal.ROIs) > 0
	if current_roi_no > CaSignal.ROI_num
		current_roi_no = CaSignal.ROI_num;
	elseif current_roi_no < 1
		current_roi_no = 1;
	end
	CaSignal.CurrentROINo = current_roi_no;
	set(handles.CurrentROINoEdit, 'String', num2str(CaSignal.CurrentROINo));
	CaSignal = Update_Image_Fcn(handles, CaSignal, true);
	CaSignal = update_subimage_show(handles, CaSignal, true);
	CaSignal = update_signal_show(handles, CaSignal);
end


% --- Executes during object creation, after setting all properties.
function CurrentROINoEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ShowROINoCheckbox.
function ShowROINoCheckbox_Callback(hObject, eventdata, handles)
global CaSignal
CaSignal = Update_Image_Fcn(handles, CaSignal, true);
CaSignal = update_subimage_show(handles, CaSignal, true);


% --- Executes on button press in NextROIButton.
function NextROIButton_Callback(hObject, eventdata, handles)
global CaSignal
CaSignal.CurrentROINo = CaSignal.CurrentROINo + 1;
if CaSignal.CurrentROINo > CaSignal.ROI_num
	CaSignal.CurrentROINo = 1;
end
set(handles.CurrentROINoEdit, 'String', num2str(CaSignal.CurrentROINo));
CaSignal = Update_Image_Fcn(handles, CaSignal, true);
CaSignal = update_subimage_show(handles, CaSignal, true);
CaSignal = update_signal_show(handles, CaSignal);



% --- Executes on button press in PreviousROIButton.
function PreviousROIButton_Callback(hObject, eventdata, handles)
global CaSignal
CaSignal.CurrentROINo = CaSignal.CurrentROINo - 1;
if CaSignal.CurrentROINo < 1
	CaSignal.CurrentROINo = CaSignal.ROI_num;
end
set(handles.CurrentROINoEdit, 'String', num2str(CaSignal.CurrentROINo));
CaSignal = Update_Image_Fcn(handles, CaSignal, true);
CaSignal = update_subimage_show(handles, CaSignal, true);
CaSignal = update_signal_show(handles, CaSignal);


function ROIDiameter_Callback(hObject, eventdata, handles)
global CaSignal
CaSignal.ROIDiameter = str2double(get(handles.ROIDiameter, 'String'));
CaSignal = update_subimage_show(handles, CaSignal, true);


% --- Executes during object creation, after setting all properties.
function ROIDiameter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PreviousFrame.
function PreviousFrame_Callback(hObject, eventdata, handles)
global CaSignal
CaSignal.current_frame = CaSignal.current_frame - 1;
if CaSignal.current_frame < 1
	CaSignal.current_frame  = CaSignal.NumOfFrames;
end
CaSignal.showing_image = CaSignal.raw_images(:, :, CaSignal.current_frame);
set(handles.FrameNoText, 'String', sprintf('%d/%d', CaSignal.current_frame, CaSignal.NumOfFrames));
CaSignal = Update_Image_Fcn(handles, CaSignal, true);
CaSignal = update_subimage_show(handles, CaSignal, true);

% --- Executes on button press in NextFrame.
function NextFrame_Callback(hObject, eventdata, handles)
global CaSignal
CaSignal.current_frame = CaSignal.current_frame + 1;
if CaSignal.current_frame > CaSignal.NumOfFrames
	CaSignal.current_frame  = 1;
end
set(handles.checkbox_mean, 'Value', 0)
set(handles.checkbox_max, 'Value', 0)
CaSignal.showing_image = CaSignal.raw_images(:, :, CaSignal.current_frame);
set(handles.FrameNoText, 'String', sprintf('%d/%d', CaSignal.current_frame, CaSignal.NumOfFrames));
CaSignal = Update_Image_Fcn(handles, CaSignal, true);
CaSignal = update_subimage_show(handles, CaSignal, true);


function FrameNoText_Callback(hObject, eventdata, handles)
global CaSignal
temp = strsplit(get(hObject, 'String'), '/');
CaSignal.current_frame = int16(str2double(temp(1)));
if CaSignal.current_frame > CaSignal.NumOfFrames
	CaSignal.current_frame  = CaSignal.NumOfFrames;
elseif CaSignal.current_frame < 1
	CaSignal.current_frame = 1;
end
if get(handles.Concatenate, 'Value') == 0
	CaSignal.showing_image = CaSignal.mean_images(:, :, CaSignal.current_frame);
else
	CaSignal.showing_image = CaSignal.mean_images(:, :, CaSignal.current_frame:CaSignal.current_frame + 2);
end
set(handles.FrameNoText, 'String', sprintf('%d/%d', CaSignal.current_frame, CaSignal.NumOfFrames));
CaSignal = Update_Image_Fcn(handles, CaSignal, true);
CaSignal = update_subimage_show(handles, CaSignal, true);


% --- Executes on button press in checkbox_mean.
function checkbox_mean_Callback(hObject, eventdata, handles)
global CaSignal
if get(handles.checkbox_mean, 'Value') == 1
	set(handles.checkbox_max, 'Value', 0)
	CaSignal.showing_image = CaSignal.mean_images;
else
	CaSignal.showing_image = CaSignal.raw_images(:, :, CaSignal.current_frame);
end
CaSignal = Update_Image_Fcn(handles, CaSignal, false);
CaSignal = update_subimage_show(handles, CaSignal, true);




% --- Executes on button press in checkbox_max.
function checkbox_max_Callback(hObject, eventdata, handles)
global CaSignal
if get(handles.checkbox_max, 'Value') == 1
	set(handles.checkbox_mean, 'Value', 0)
	CaSignal.showing_image = CaSignal.max_images;
else
	CaSignal.showing_image = CaSignal.raw_images(:, :, CaSignal.current_frame);
end
CaSignal = Update_Image_Fcn(handles, CaSignal, false);
CaSignal = update_subimage_show(handles, CaSignal, true);



% --- Executes on button press in checkbox13.
function checkbox13_Callback(hObject, eventdata, handles)
global CaSignal
a = CaSignal;


% --- Executes on button press in BatchButton.
function BatchButton_Callback(hObject, eventdata, handles)
global CaSignal
disp('Batching ...')
extract_calcium_signal(CaSignal);
disp('Done')

