function varargout = Record(varargin)
% ʹ�õ�����һ����豸����Matlab���ݲɼ�����
% RECORD MATLAB code for Record.fig
%      RECORD, by itself, creates a new RECORD or raises the existing
%      singleton*.
%
%      H = RECORD returns the handle to a new RECORD or the handle to
%      the existing singleton*.
%
%      RECORD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RECORD.M with the given input arguments.
%
%      RECORD('Property','Value',...) creates a new RECORD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Record_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Record_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Record

% Last Modified by GUIDE v2.5 16-May-2016 22:11:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Record_OpeningFcn, ...
                   'gui_OutputFcn',  @Record_OutputFcn, ...
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

% --- Executes just before Record is made visible.
function Record_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Record (see VARARGIN)

% Choose default command line output for Record
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using Record.
if strcmp(get(hObject,'Visible'),'off')
    % plot(rand(5));
end

global dirty
dirty = 0;					% �ĵ�����

% UIWAIT makes Record wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Record_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});



function edit_record_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_record_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_record_time as text
%        str2double(get(hObject,'String')) returns contents of edit_record_time as a double


% --- Executes during object creation, after setting all properties.
function edit_record_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_record_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_filename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filename as text
%        str2double(get(hObject,'String')) returns contents of edit_filename as a double


% --- Executes during object creation, after setting all properties.
function edit_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_start.
function pushbutton_start_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dirty
if dirty == 1
	choice = questdlg('������δ���棬�����ɼ��ᶪʧǰ��ɼ������ݡ��Ƿ���Ҫ�������ݲ���ʼ�ɼ������ݣ�','����', 'Yes','No', 'No');
	if strcmp(choice,'No')
		return
	end
end
save('tempData', 'hObject', 'eventdata', 'handles', 'dirty');
clear all
load('tempData', 'hObject', 'eventdata', 'handles', 'dirty');
global data port ExtAdc freeze windowStartPosition windowEndPosition endPoint ax startMark Running Yscale freqScale h
% 0. ���ó�ֵ
startMark = 1;								% ��������ʶ�����ʼ���
set(handles.pushbutton_stop,'UserData', 0);	% ֹͣ����
set(handles.slider_data_position, 'Value', 1);		% �������ŵ����
if isempty(get(handles.pushbutton_spectrom_hide , 'UserData'))
	set(handles.pushbutton_spectrom_hide , 'UserData', 1);	% ����ʱƵͼ��ʾ
end
Yscale = 100;								% ��ʾ��Χ��100uV
freqScale = 0.12;							% Ƶ����ʾ��31.25Hz��0.25��500Hz��
Running = 0;								% �������еı��
% 1. ���豸
delete(instrfindall);						% �ȹر����д���
com = ComboQuery;							% comΪ�ַ������飬�������еĵ�����һ����豸����
h = ComboOpen(com(1,:));					% �򿪵�һ���豸

% 2. �ɼ����ݡ���ο����ظ���ν��С�
ax(1) = handles.axes_data;		hold off;
ax(2) = handles.axes_ext_adc;	hold off;
ax(3) = handles.axes_port;		hold off;
linkaxes(ax,'x')							% ��������
nDuration = str2double(get(handles.edit_record_time, 'String'));	% ȡ�òɼ�ʱ�䣬��λ��
nLength = (nDuration+1) * 30000;
data = zeros(nLength, 1);
port = zeros(nLength, 1, 'uint8');
ExtAdc = zeros(nLength, 2);
t1 = tic;
startPoint = 1;
Running = 1;								% �������еı��
dirty = 1;									% �������ݵı�־
while 1
	len = ComboGetLength(h);				% ȡ��һ������Ѿ��ɼ������ݳ���
	[A, E, P] = ComboGetData(h, len);		% �����������ݳ��ȣ���������
	endPoint = startPoint + len - 1;		% ��¼����ĩβ��λ��
	% data = [data, double(A).*9.9341e-09];	% 1 LSB = ��2.5V/30/2^24
	%ExtAdc = [ExtAdc; double((2048 - E') .* 10) / 4095];	% ���ȡ�10V������Ŵ�����
	data(startPoint:endPoint) = double(A).*9.9341e-09;	% 1 LSB = ��2.5V/30/2^24
	port(startPoint:endPoint) = P;
	ExtAdc(startPoint:endPoint,:) = double((2048 - E') .* 10) / 4095;	% ���ȡ�10V������Ŵ�����
	t2 = toc(t1);
	% ����̶���ʾ�������¼��㲢����sliderλ��
	if (freeze)
		% ����������Ƶ���ı�������������������Ϣ�õ�����ֵ��Զ��������ı���λ��
	else
		windowStartPosition = max(1, endPoint-3*30000);			% ��ʾ���3�������
		windowEndPosition = endPoint;							% ��ʾ���3�������
		set(handles.text_time, 'String', strcat(num2str(t2), ' s'));
	end
	plotData(handles);
	if (t2 >= nDuration)										% �ɼ�nDuration��
		break;
	end
	if (get(handles.pushbutton_stop, 'UserData') == 1)			% ������ֹͣ��
		break;
	end
	pause(0.001)												% �������ȣ���������������������
	startPoint = endPoint + 1;
end
Running = 0;								% �������еı��
% linkaxes(ax,'x')							% ��������
% �ض����ݣ�����δ�ɼ��Ŀ�����
data = data(1:endPoint);
ExtAdc = ExtAdc(1:endPoint, :);
port = port(1:endPoint);
% �����ݱ�������һ�ؼ����˴�Ϊstart��ť�����£��Է���������֮�����ݶ�ʧ
% ע�⣺���ܱ�����base workspace�У���Ϊ����ֻ����matlab����������m�ű���Ҫ����.exe��ִ���ļ�����base workspace�����ڣ����ݾͶ�ʧ�ˡ�
setappdata(handles.pushbutton_start,'data',data);
setappdata(handles.pushbutton_start,'ExtAdc',ExtAdc);
setappdata(handles.pushbutton_start,'port',port);

% 3. �ر��豸���˳�
ComboClose(h)								% �ر��豸
delete('tempData.mat');


% �������ź�����ͼ�ĺ�������Slider��Start�ص���������
function plotData(handles)
global data port ExtAdc windowStartPosition windowEndPosition ax Yscale freqScale h
try
Ncomp = int32(floor(single(windowEndPosition - windowStartPosition) / 30));	% 30������ѹ�����������
dataRange = windowEndPosition-Ncomp*30+1:windowEndPosition;					% ��ʾ�����ݷ�Χ
viewEEG = mean(reshape(data(dataRange), 30, Ncomp), 1);			% EEG����ѹ��30��
viewExt = squeeze(mean(reshape(ExtAdc(dataRange,:), 30, Ncomp, 2), 1));	% ͬ���ź�ѹ��30��
plot(ax(1), viewEEG)											% ��ʾADC24����
ylim(ax(1), [-Yscale*1e-6 Yscale*1e-6]);
plot(ax(2), viewExt)											% ��ʾͬ���ź�����
% clear DI;
DI(:,1) = bitget(port(dataRange), 1) * 0.8 + 0;					% ȡ��PB9(port�˿ڵ�0ͨ��)
DI(:,2) = bitget(port(dataRange), 2) * 0.8 + 1;					% ȡ��PC6(port�˿ڵ�1ͨ��)
DI(:,3) = bitget(port(dataRange), 3) * 0.8 + 2;					% ȡ��PC7(port�˿ڵ�2ͨ��)
DI(:,4) = bitget(port(dataRange), 4) * 0.8 + 3;					% ȡ��PC8(port�˿ڵ�3ͨ��)
plot(ax(3), squeeze(mean(reshape(DI, [30 Ncomp 4]), 1)));		% ���ƶ˿�����
ylim(ax(3), [0 4])
hold off
% len = (floor(length(viewEEG)/4))*4;								% ��ʱƵͼ����������4����
% specEEG = mean(reshape(viewEEG(1:len), 4, len/4), 1);			%  �ز�����250 sps
% if length(specEEG) >= 64
% 	[S,F,T,P] = spectrogram(specEEG,64,63,1024,250);			% ����ʱƵͼ
% 	freqTop = floor(length(F) * freqScale);
% 	surf(handles.axes_spectrom, T,F(1:freqTop),10*log10(P(1:freqTop,:)),'edgecolor','none'), axis tight, view(0,90)	% ��ʾ
% end
bIsShow = get(handles.pushbutton_spectrom_hide , 'UserData');
if bIsShow && (length(viewEEG) >= 512)
	[S,F,T,P] = spectrogram(viewEEG,512,510,2048,1E3);			% ����ʱƵͼ
	freqTop = floor(length(F) * freqScale);
	surf(handles.axes_spectrom, T,F(1:freqTop),10*log10(P(1:freqTop,:)),'edgecolor','none'), axis tight, view(0,90)	% ��ʾ
end
catch err
	ComboClose(h)											% ��������ʱ�ر��豸
	err
	error(err)
end

% --- Executes on slider movement.
function slider_data_position_Callback(hObject, eventdata, handles)
% hObject    handle to slider_data_position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global data port ExtAdc freeze windowStartPosition windowEndPosition endPoint startMark Running
persistent lastPos lastEndPoint;
if startMark										% ��һ�λع�ʱ�����þ�̬������ֵ
	startMark = 0;
	lastPos = 1;
	lastEndPoint = 1;
end
if isempty(lastEndPoint)
	lastEndPoint = 1;
end
if isempty(freeze)
	freeze = 0;
end
if isempty(lastPos)
	lastPos = 1;
end
% position��Χ��[0..1]�������ͷ�Ĳ���0.001��
pos = get(hObject,'Value');
% ���¹�����λ��
% if endPoint > 2*30000
% 	position = min(1.0, double(windowStartPosition) / double(endPoint-2*30000));	% ÿ�βɼ������µ�ǰ��ʾ���������������е�λ��
% 	set(handles.slider_data_position, 'Value', position);
% end
% set(hObject, 'Value', (minor+maxor)/2);
if pos==1
	freeze = 0;										% ����ʱ��ı�־��=0����������ˢ�´���
	windowStartPosition = max(1, endPoint-2*30000);	% ��ʾ���2�������
	windowEndPosition = endPoint;					% ��ʾ���2�������
else
	if freeze ~= 1							% �ոս��յ�����ָ��
		freeze = 1;								% ����ʱ��ı�־��=1������ˢ�´���
		lastPos = 1;
		lastEndPoint = endPoint;
	end
	%else									% ��ǰ�������̵ģ����ڸ��¡�Ҫ�����Ƿ����˼�ͷ
		lastPoint = lastPos * lastEndPoint;		% �����ϴ�λ�õĵ�������
		thisPoint = pos * endPoint;				% ���㱾��λ�õĵ�������
		posChange = double(abs(thisPoint-lastPoint))/double(endPoint);	% λ�øı�İٷֱ�
		if posChange < 0.002					% ͨ�������ͷ�ı��λ��
			if thisPoint < lastPoint				% ��ǰ��ƪ
				windowStartPosition = max(1, windowStartPosition-0.5*30000);	% һ�η�0.5�������
			else									% ���ƪ
				windowStartPosition = min(endPoint-2*30000+1, windowStartPosition+0.5*30000);	% һ�η�0.5�������
			end
			windowEndPosition = min(endPoint, windowStartPosition+2*30000);
			if endPoint > 2*30000
				pos = min(1, double(windowStartPosition) / double(endPoint-2*30000));	% ���µ�ǰ��ʾ���������������е�λ��
				set(hObject,'Value', pos);
			end
		else									% ͨ���϶������ı��λ��
			if endPoint > 2*30000
				windowStartPosition = (endPoint-2*30000) * pos;
				windowEndPosition = windowStartPosition + 2*30000;
			else
				windowStartPosition = 1;
				windowEndPosition = endPoint;
			end
		end
		lastPos = pos;
	%end
end
set(handles.text_time, 'String', strcat(num2str(double(windowStartPosition)/30000.0), ' s'));
lastEndPoint = endPoint;
plotData(handles);

% --- Executes on button press in pushbutton_dir.
function pushbutton_dir_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder_name = uigetdir;
set(handles.pushbutton_dir,'UserData',folder_name);


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dirty
folder_name = strcat(get(handles.pushbutton_dir, 'UserData'), '\');
file_name = strcat(folder_name, get(handles.edit_filename, 'String'));
data = getappdata(handles.pushbutton_start,'data');
ExtAdc = getappdata(handles.pushbutton_start,'ExtAdc');
port = getappdata(handles.pushbutton_start,'port');
uisave({'data', 'ExtAdc', 'port'}, file_name);
dirty = 0;

% --- Executes on button press in pushbutton_scale_up.
function pushbutton_scale_up_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_scale_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Yscale Running
switch Yscale
	case 1
		Yscale = 1;
	case 2
		Yscale = 1;
	case 5
		Yscale = 2;
	case 10
		Yscale = 5;
	case 20
		Yscale = 10;
	case 50
		Yscale = 20;
	case 100
		Yscale = 50;
	case 200
		Yscale = 100;
	case 500
		Yscale = 200;
	case 1000
		Yscale = 500;
	case 2000
		Yscale = 1000;
	case 5000
		Yscale = 2000;
	case 10000
		Yscale = 5000;
	case 20000
		Yscale = 10000;
	case 50000
		Yscale = 20000;
	case 100000
		Yscale = 50000;
end
if Running == 0
	plotData(handles);
end

% --- Executes on button press in pushbutton_scale_down.
function pushbutton_scale_down_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_scale_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Yscale Running
switch Yscale
	case 1
		Yscale = 2;
	case 2
		Yscale = 5;
	case 5
		Yscale = 10;
	case 10
		Yscale = 20;
	case 20
		Yscale = 50;
	case 50
		Yscale = 100;
	case 100
		Yscale = 200;
	case 200
		Yscale = 500;
	case 500
		Yscale = 1000;
	case 1000
		Yscale = 2000;
	case 2000
		Yscale = 5000;
	case 5000
		Yscale = 10000;
	case 10000
		Yscale = 20000;
	case 20000
		Yscale = 50000;
	case 50000
		Yscale = 100000;
	case 100000
		Yscale = 100000;
end
if Running == 0
	plotData(handles);
end

% --- Executes during object creation, after setting all properties.
function slider_data_position_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_data_position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton_stop.
function pushbutton_stop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton_stop,'UserData',1);


% --- Executes on button press in pushbutton_freq_expand.
function pushbutton_freq_expand_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_freq_expand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global freqScale Running
if freqScale < 1
	freqScale = freqScale / 0.8;
end
if freqScale > 1
	freqScale = 1.0;
end
if Running == 0
	plotData(handles);
end

% --- Executes on button press in pushbutton_freq_shrink.
function pushbutton_freq_shrink_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_freq_shrink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global freqScale Running
if freqScale > 0.01
	freqScale = freqScale * 0.8;
end
if freqScale < 0.01
	freqScale = 0.01;
end
if Running == 0
	plotData(handles);
end


% --- Executes on button press in pushbutton_spectrom_hide.
function pushbutton_spectrom_hide_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_spectrom_hide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
bIsShow = get(hObject, 'UserData');
if (bIsShow == 0)
	bIsShow = 1;
	set(hObject, 'String', 'ͣ��');
else
	bIsShow = 0;
	set(hObject, 'String', 'ˢ��');
end
set(handles.pushbutton_spectrom_hide , 'UserData', bIsShow);
