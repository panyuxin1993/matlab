function varargout = AP(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AP_OpeningFcn, ...
                   'gui_OutputFcn',  @AP_OutputFcn, ...
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

function AP_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
c=get(handles.current,'Value');
d=get(handles.duration,'Value');
tt=get(handles.time,'Value');
int=get(handles.interval,'Value');
set(handles.text1,'String',[num2str(tt),'ms']);
set(handles.text2,'String',[num2str(d),'ms']);
set(handles.text3,'String',[num2str(int),'ms']);
set(handles.text4,'String',[num2str(c),'mA']);
ENa=60;EK=-74.7;EL=-65.8;Erest=-68;vm(1)=Erest;i=0;
gNabar=120;gKbar=30;n(1)=0.3;m(1)=0.0665;h(1)=0.6;
GK(1)=12;GNa(1)=30;gl=0.001;Cm=1;dt=0.01;
for t=0:0.01:tt
    i=i+1;
    v=vm(i)-Erest;
    
    aln=0.01*(10-v)/(exp(10-v)/10-1);
    betn=0.125*exp(-v/80);
    alm=0.1*(25-v)/(exp(25-v)/10-1);
    betm=4*exp(-v/18);
    alh=0.07*exp(-v/20);
    beth=1/(exp(30-v)/10+1);
    
    n=n+(aln*(1-n)-betn*n)*dt;
    m=m+(alm*(1-m)-betm*m)*dt;
    h=h+(alh*(1-h)-beth*h)*dt;
    nn(i)=n;mm(i)=m;hh(i)=h;
    gK(i)=GK*n^4;
    gNa(i)=GNa*m^3*h;
    IK(i)=gK(i)*(vm(i)-EK);
    INa(i)=gNa(i)*(vm(i)-ENa);
    IL(i)=gl*(vm(i)-EL);
    if mod(t,int)<d
        Id(i)=c;
    else
        Id(i)=0;
    end
    vm(i+1)=vm(i)+(Id(i)-(INa(i)+IK(i)+IL(i)))*dt/Cm;
end
z=0:0.01:tt; 

plot(handles.axes1,z,Id,'r',z,vm(1:end-1),'b','linewidth',1)
xlabel(handles.axes1,'Time (ms)')
ylabel(handles.axes1,'V_M , I_{Inj}')
axis(handles.axes1,[0 tt min(vm)-10 c+10])

function varargout = AP_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function time_Callback(hObject, eventdata, handles)
c=get(handles.current,'Value');
d=get(handles.duration,'Value');
tt=get(handles.time,'Value');
int=get(handles.interval,'Value');
set(handles.text1,'String',[num2str(tt),'ms']);
set(handles.text2,'String',[num2str(d),'ms']);
set(handles.text3,'String',[num2str(int),'ms']);
set(handles.text4,'String',[num2str(c),'mA']);
ENa=60;EK=-74.7;EL=-65.8;Erest=-68;vm(1)=Erest;i=0;
gNabar=120;gKbar=30;n(1)=0.3;m(1)=0.0665;h(1)=0.6;
GK(1)=12;GNa(1)=30;gl=0.001;Cm=1;dt=0.01;
for t=0:0.01:tt
    i=i+1;
    v=vm(i)-Erest;
    
    aln=0.01*(10-v)/(exp(10-v)/10-1);
    betn=0.125*exp(-v/80);
    alm=0.1*(25-v)/(exp(25-v)/10-1);
    betm=4*exp(-v/18);
    alh=0.07*exp(-v/20);
    beth=1/(exp(30-v)/10+1);
    
    n=n+(aln*(1-n)-betn*n)*dt;
    m=m+(alm*(1-m)-betm*m)*dt;
    h=h+(alh*(1-h)-beth*h)*dt;
    nn(i)=n;mm(i)=m;hh(i)=h;
    gK(i)=GK*n^4;
    gNa(i)=GNa*m^3*h;
    IK(i)=gK(i)*(vm(i)-EK);
    INa(i)=gNa(i)*(vm(i)-ENa);
    IL(i)=gl*(vm(i)-EL);
    if mod(t,int)<d
        Id(i)=c;
    else
        Id(i)=0;
    end
    vm(i+1)=vm(i)+(Id(i)-(INa(i)+IK(i)+IL(i)))*dt/Cm;
end
z=0:0.01:tt; 

plot(handles.axes1,z,Id,'r',z,vm(1:end-1),'b','linewidth',1)
xlabel(handles.axes1,'Time (ms)')
ylabel(handles.axes1,'V_M , I_{Inj}')
axis(handles.axes1,[0 tt min(vm)-10 c+10])

function time_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function duration_Callback(hObject, eventdata, handles)

c=get(handles.current,'Value');
d=get(handles.duration,'Value');
tt=get(handles.time,'Value');
int=get(handles.interval,'Value');
set(handles.text1,'String',tt);
set(handles.text2,'String',d);
set(handles.text3,'String',int);
set(handles.text4,'String',c);
ENa=60;EK=-74.7;EL=-65.8;Erest=-68;vm(1)=Erest;i=0;
gNabar=120;gKbar=30;n(1)=0.3;m(1)=0.0665;h(1)=0.6;
GK(1)=12;GNa(1)=30;gl=0.001;Cm=1;dt=0.01;
for t=0:0.01:tt
    i=i+1;
    v=vm(i)-Erest;
    
    aln=0.01*(10-v)/(exp(10-v)/10-1);
    betn=0.125*exp(-v/80);
    alm=0.1*(25-v)/(exp(25-v)/10-1);
    betm=4*exp(-v/18);
    alh=0.07*exp(-v/20);
    beth=1/(exp(30-v)/10+1);
    
    n=n+(aln*(1-n)-betn*n)*dt;
    m=m+(alm*(1-m)-betm*m)*dt;
    h=h+(alh*(1-h)-beth*h)*dt;
    nn(i)=n;mm(i)=m;hh(i)=h;
    gK(i)=GK*n^4;
    gNa(i)=GNa*m^3*h;
    IK(i)=gK(i)*(vm(i)-EK);
    INa(i)=gNa(i)*(vm(i)-ENa);
    IL(i)=gl*(vm(i)-EL);
    if mod(t,int)<d
        Id(i)=c;
    else
        Id(i)=0;
    end
    vm(i+1)=vm(i)+(Id(i)-(INa(i)+IK(i)+IL(i)))*dt/Cm;
end
z=0:0.01:tt; 

plot(handles.axes1,z,Id,'r',z,vm(1:end-1),'b','linewidth',1)
xlabel(handles.axes1,'Time (ms)')
ylabel(handles.axes1,'V_M , I_{Inj}')
axis(handles.axes1,[0 tt min(vm)-10 c+10])

function duration_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function interval_Callback(hObject, eventdata, handles)

c=get(handles.current,'Value');
d=get(handles.duration,'Value');
tt=get(handles.time,'Value');
int=get(handles.interval,'Value');
set(handles.text1,'String',tt);
set(handles.text2,'String',d);
set(handles.text3,'String',int);
set(handles.text4,'String',c);
ENa=60;EK=-74.7;EL=-65.8;Erest=-68;vm(1)=Erest;i=0;
gNabar=120;gKbar=30;n(1)=0.3;m(1)=0.0665;h(1)=0.6;
GK(1)=12;GNa(1)=30;gl=0.001;Cm=1;dt=0.01;
for t=0:0.01:tt
    i=i+1;
    v=vm(i)-Erest;
    
    aln=0.01*(10-v)/(exp(10-v)/10-1);
    betn=0.125*exp(-v/80);
    alm=0.1*(25-v)/(exp(25-v)/10-1);
    betm=4*exp(-v/18);
    alh=0.07*exp(-v/20);
    beth=1/(exp(30-v)/10+1);
    
    n=n+(aln*(1-n)-betn*n)*dt;
    m=m+(alm*(1-m)-betm*m)*dt;
    h=h+(alh*(1-h)-beth*h)*dt;
    nn(i)=n;mm(i)=m;hh(i)=h;
    gK(i)=GK*n^4;
    gNa(i)=GNa*m^3*h;
    IK(i)=gK(i)*(vm(i)-EK);
    INa(i)=gNa(i)*(vm(i)-ENa);
    IL(i)=gl*(vm(i)-EL);
    if mod(t,int)<d
        Id(i)=c;
    else
        Id(i)=0;
    end
    vm(i+1)=vm(i)+(Id(i)-(INa(i)+IK(i)+IL(i)))*dt/Cm;
end
z=0:0.01:tt; 

plot(handles.axes1,z,Id,'r',z,vm(1:end-1),'b','linewidth',1)
xlabel(handles.axes1,'Time (ms)')
ylabel(handles.axes1,'V_M , I_{Inj}')
axis(handles.axes1,[0 tt min(vm)-10 c+10])

function interval_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function current_Callback(hObject, eventdata, handles)

c=get(handles.current,'Value');
d=get(handles.duration,'Value');
tt=get(handles.time,'Value');
int=get(handles.interval,'Value');
set(handles.text1,'String',[num2str(tt),'ms']);
set(handles.text2,'String',[num2str(d),'ms']);
set(handles.text3,'String',[num2str(int),'ms']);
set(handles.text4,'String',[num2str(c),'mA']);
ENa=60;EK=-74.7;EL=-65.8;Erest=-68;vm(1)=Erest;i=0;
gNabar=120;gKbar=30;n(1)=0.3;m(1)=0.0665;h(1)=0.6;
GK(1)=12;GNa(1)=30;gl=0.001;Cm=1;dt=0.01;
for t=0:0.01:tt
    i=i+1;
    v=vm(i)-Erest;
    
    aln=0.01*(10-v)/(exp(10-v)/10-1);
    betn=0.125*exp(-v/80);
    alm=0.1*(25-v)/(exp(25-v)/10-1);
    betm=4*exp(-v/18);
    alh=0.07*exp(-v/20);
    beth=1/(exp(30-v)/10+1);
    
    n=n+(aln*(1-n)-betn*n)*dt;
    m=m+(alm*(1-m)-betm*m)*dt;
    h=h+(alh*(1-h)-beth*h)*dt;
    nn(i)=n;mm(i)=m;hh(i)=h;
    gK(i)=GK*n^4;
    gNa(i)=GNa*m^3*h;
    IK(i)=gK(i)*(vm(i)-EK);
    INa(i)=gNa(i)*(vm(i)-ENa);
    IL(i)=gl*(vm(i)-EL);
    if mod(t,int)<d
        Id(i)=c;
    else
        Id(i)=0;
    end
    vm(i+1)=vm(i)+(Id(i)-(INa(i)+IK(i)+IL(i)))*dt/Cm;
end
z=0:0.01:tt; 

plot(handles.axes1,z,Id,'r',z,vm(1:end-1),'b','linewidth',1)
xlabel(handles.axes1,'Time (ms)')
ylabel(handles.axes1,'V_M , I_{Inj}')
axis(handles.axes1,[0 tt min(vm)-10 c+10])

function current_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
