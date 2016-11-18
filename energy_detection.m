function varargout = energy_detection(varargin)
% ENERGY_DETECTION MATLAB code for energy_detection.fig
%      ENERGY_DETECTION, by itself, creates a new ENERGY_DETECTION or raises the existing
%      singleton*.
%
%      H = ENERGY_DETECTION returns the handle to a new ENERGY_DETECTION or the handle to
%      the existing singleton*.
%
%      ENERGY_DETECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ENERGY_DETECTION.M with the given input arguments.
%
%      ENERGY_DETECTION('Property','Value',...) creates a new ENERGY_DETECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before energy_detection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to energy_detection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help energy_detection

% Last Modified by GUIDE v2.5 25-Aug-2015 13:32:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @energy_detection_OpeningFcn, ...
                   'gui_OutputFcn',  @energy_detection_OutputFcn, ...
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


% --- Executes just before energy_detection is made visible.
function energy_detection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to energy_detection (see VARARGIN)

% Choose default command line output for energy_detection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
clc;
% files_dir = uigetdir('matlabroot','Select the measurement files directory');
% save('Files_Directory','files_dir');

set(handles.popupmenu_N,'Value',2);
% UIWAIT makes energy_detection wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = energy_detection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu_txsignal.
function popupmenu_txsignal_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_txsignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_txsignal contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_txsignal


% --- Executes during object creation, after setting all properties.
function popupmenu_txsignal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_txsignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_channel.
function popupmenu_channel_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_channel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_channel


% --- Executes during object creation, after setting all properties.
function popupmenu_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_N.
function popupmenu_N_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_N contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_N


% --- Executes during object creation, after setting all properties.
function popupmenu_N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_snr.
function popupmenu_snr_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_snr contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_snr


% --- Executes during object creation, after setting all properties.
function popupmenu_snr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_awgn.
function checkbox_awgn_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_awgn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_awgn


% --- Executes on button press in checkbox_slc2.
function checkbox_slc2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_slc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_slc2


% --- Executes on button press in checkbox_sls2.
function checkbox_sls2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_sls2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_sls2


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
TXsignal = get(handles.popupmenu_txsignal,'Value');
channel = get(handles.popupmenu_channel,'Value');
RXtech(1) = get(handles.checkbox_awgn,'Value');
RXtech(2) = get(handles.checkbox_slc2,'Value');
RXtech(3) = get(handles.checkbox_slc4,'Value');
RXtech(4) = get(handles.checkbox_sls2,'Value');
RXtech(5) = get(handles.checkbox_sls4,'Value');
N = 10^get(handles.popupmenu_N,'Value');
SNR = 5*(get(handles.popupmenu_snr,'Value')-2);
run_ED(TXsignal, channel, RXtech, N, SNR)


% --- Executes on button press in checkbox_sls4.
function checkbox_sls4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_sls4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_sls4


% --- Executes on button press in checkbox_slc4.
function checkbox_slc4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_slc4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_slc4
