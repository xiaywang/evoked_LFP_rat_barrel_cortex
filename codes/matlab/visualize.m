function varargout = visualize(varargin)
% visualize MATLAB code for visualize.fig
%      visualize, by itself, creates a new visualize or raises the existing
%      singleton*.
%
%      H = visualize returns the handle to a new visualize or the handle to
%      the existing singleton*.
%
%      visualize('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in visualize.M with the given input arguments.
%
%      visualize('Property','Value',...) creates a new visualize or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before visualize_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to visualize_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help visualize

% Last Modified by GUIDE v2.5 03-Apr-2018 11:15:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @visualize_OpeningFcn, ...
                   'gui_OutputFcn',  @visualize_OutputFcn, ...
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


% --- Executes just before visualize is made visible.
function visualize_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to visualize (see VARARGIN)

% Choose default command line output for visualize
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes visualize wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = visualize_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in next_botton.
function next_botton_Callback(hObject, eventdata, handles)
% hObject    handle to next_botton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% contents = get(handles.select_data,'String'); 
% selected_data = contents{get(handles.select_data,'Value')};
% load(selected_data)

curtime = str2num(get(handles.current_t, 'string'));
updated_time = curtime + 1;
handles.updated_time = updated_time;
guidata(hObject, handles);
max_t = getappdata(handles.select_data, 'max_t');
if updated_time <= max_t
    plot_image(hObject, eventdata, handles)
    calculate_mean(hObject, eventdata, handles)
    plot_hist(hObject, eventdata, handles)
    set(handles.current_t, 'string', num2str(updated_time))
else
    disp('no more time points')
end


% --- Executes on button press in previous_botton.
function previous_botton_Callback(hObject, eventdata, handles)
% hObject    handle to previous_botton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% contents = get(handles.select_data,'String'); 
% selected_data = contents{get(handles.select_data,'Value')};
% load(selected_data)

curtime = str2num(get(handles.current_t, 'string'));
updated_time = curtime - 1;
handles.updated_time = updated_time;
guidata(hObject, handles);
if updated_time >= 1
    plot_image(hObject, eventdata, handles)
    calculate_mean(hObject, eventdata, handles)
    plot_hist(hObject, eventdata, handles)
    set(handles.current_t, 'string', num2str(updated_time))
else
    disp('no more time points')
end



function plot_image(hObject, eventdata, handles)

time_point = handles.updated_time;
% data = handles.data;
data = getappdata(handles.select_data , 'data');
I = data(:, :, time_point);
if get(handles.checkbox1, 'Value') == 1.0
    sigma = str2num(get(handles.sigma, 'string'));
    I = imgaussfilt(I, sigma);
end
maxIntensity = max(max(max(data)));
minIntensity = min(min(min(data)));
% BW = edge(I,'sobel', 0.5);
if (get(handles.radiobutton1,'Value') == get(handles.radiobutton1,'Max'))
    set(handles.uipanel2,'visible','on')
    set(handles.uipanel1,'visible','off')
    axes(handles.axes4)
    imshow(I, [], 'Colormap', flipud(parula), 'InitialMagnification','fit') %, 'Parent',handles.axes1
    colorbar
    caxis([minIntensity maxIntensity]);
    axes(handles.axes5)
    time = getappdata(handles.select_data, 'time');
    max_t = getappdata(handles.select_data, 'max_t');
    mean1 = zeros(max_t,1);
    for t = 1:size(data,3)
        mean1(t) = mean2(data(:,:,t));
    end
    plot(time, mean1, 'linewidth', 1.5)
    xlabel('Time (ms)')
    ylabel('Amplitude (mV)')
    hold on
    plot(time(time_point), mean1(time_point), 'or', 'markerfacecolor', 'r')
    hold off
    legend('Average in time', 'Current time point', 'location', 'best')
else
    set(handles.uipanel2,'visible','off')
    set(handles.uipanel1,'visible','on')
    axes(handles.axes1)
    imshow(I, [], 'Colormap', flipud(parula), 'InitialMagnification','fit') %, 'Parent',handles.axes1
    % hold on
    % imshow(BW, [])
    colorbar
    caxis([minIntensity maxIntensity]);
end


function calculate_mean(hObject, eventdata, handles)
    data = getappdata(handles.select_data , 'data');
    time_point = handles.updated_time;
    I = data(:, :, time_point);
    if get(handles.checkbox1, 'Value') == 1.0
        sigma = str2num(get(handles.sigma, 'string'));
        I = imgaussfilt(I, sigma);
    end
    if (get(handles.radiobutton2,'Value') == get(handles.radiobutton2,'Max'))
        text_mean = mean(mean(I)); text_std = std(std(I));
        set(handles.text_mean, 'string', num2str(text_mean))
        set(handles.text_std, 'string', num2str(text_std))
    else
        if (get(handles.radiobutton3,'Value') == get(handles.radiobutton3,'Max'))
            text_mean = mean(mean(I)); text_std = std(std(I));
            set(handles.text_mean, 'string', num2str(text_mean))
            set(handles.text_std, 'string', num2str(text_std))
        else
            set(handles.text_mean2, 'visible', 'on')
            set(handles.text_std2, 'visible', 'on')
            text_mean = mean(mean(I(1:8,:))); text_std = std(std(I(1:8,:)));
            text_mean2 = mean(mean(I(9:16,:))); text_std2 = std(std(I(9:16,:)));
            set(handles.text_mean, 'string', num2str(text_mean))
            set(handles.text_std, 'string', num2str(text_std))
            set(handles.text_mean2, 'string', num2str(text_mean2))
            set(handles.text_std2, 'string', num2str(text_std2))
        end
    end


function plot_hist(hObject, eventdata, handles)

time_point = handles.updated_time;
% data = handles.data;
data = getappdata(handles.select_data , 'data');
I = data(:, :, time_point);
maxIntensity = max(max(max(data)));
minIntensity = min(min(min(data)));
% E = linspace(maxIntensity, minIntensity, 10);
% disp(E)
nbins = str2num(get(handles.n_bins, 'string'));

if (get(handles.radiobutton2,'Value') == get(handles.radiobutton2,'Max'))
    axes(handles.axes2);
    histogram(I, nbins) %, 'Parent',handles.axes2
else
    if (get(handles.radiobutton3,'Value') == get(handles.radiobutton3,'Max'))
        axes(handles.axes2);
        histogram(I, nbins) %, 'Parent',handles.axes2
    else
        set(handles.axes3, 'visible', 'on')
        I1 = I(1:8,:); I2 = I(9:16,:);
        axes(handles.axes2);
        histogram(I1, nbins);
        axes(handles.axes3);
        histogram(I2, nbins);
    end
end


function current_t_Callback(hObject, eventdata, handles)
% hObject    handle to current_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of current_t as text
%        str2double(get(hObject,'String')) returns contents of current_t as a double

% contents = get(handles.select_data,'String'); 
% selected_data = contents{get(handles.select_data,'Value')};
% load(selected_data)

updated_time = str2num(get(hObject, 'String'));
handles.updated_time = updated_time;
guidata(hObject, handles);
max_t = getappdata(handles.select_data, 'max_t');
% axes(handles.axes1);
if updated_time <= max_t || updated_time > 0
    plot_image(hObject, eventdata, handles)
    calculate_mean(hObject, eventdata, handles)
    plot_hist(hObject, eventdata, handles)
else
    disp('invalid time point')
end


% --- Executes during object creation, after setting all properties.
function current_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in select_data.
function select_data_Callback(hObject, eventdata, handles)
% hObject    handle to select_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns select_data contents as cell array
%        contents{get(hObject,'Value')} returns selected item from select_data

contents = get(handles.select_data,'String');
selected_data = contents{get(handles.select_data,'Value')};
folderPath = get(handles.browse_folder, 'String');
% load(selected_data, 'data');
% setappdata(handles.select_data, 'data', data);

if (get(handles.radiobutton1,'Value') == get(handles.radiobutton1,'Max'))
    selected_data = strcat(folderPath, '\', selected_data);
	load(selected_data);
    setappdata(handles.select_data, 'data', data);
else
    folderName = getappdata(handles.browse_botton, 'folderName');
    directory = dir(fullfile(folderName, strcat('*#', num2str(selected_data),'.mat')));
    fileNames = {directory.name}';
    trial_data = [];
    for i = 1:length(fileNames)
        load(strcat(folderPath, '\', fileNames{i}));
        disp(fileNames{i})
        trial_data = cat(1, trial_data, data);
    end
    setappdata(handles.select_data, 'data', trial_data);
end
setappdata(handles.select_data, 'max_t', size(data, 3));
setappdata(handles.select_data, 'time', time);
set(handles.current_t, 'string', '0')

% rmappdata(handles.GUIHandle, 'yourVariable')


% --- Executes during object creation, after setting all properties.
function select_data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function browse_folder_Callback(hObject, eventdata, handles)
% hObject    handle to browse_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of browse_folder as text
%        str2double(get(hObject,'String')) returns contents of browse_folder as a double


% --- Executes during object creation, after setting all properties.
function browse_folder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to browse_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_botton.
function browse_botton_Callback(hObject, eventdata, handles)
% hObject    handle to browse_botton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

startingFolder = pwd;
folderName = uigetdir(startingFolder, 'Select a folder containing the data');
set(handles.browse_folder, 'string', folderName);
setappdata(handles.browse_botton, 'folderName', folderName);
directory = dir(fullfile(folderName,'*.mat'));
files = {directory.name}';
setappdata(handles.browse_botton, 'files', files);
% set(handles.select_data, 'string', files)

if (get(handles.radiobutton1,'Value') == get(handles.radiobutton1,'Max'))
    files(2:end+1) = files(1:end);
    files(1) = {'Please select data'};
	set(handles.select_data, 'string', files)
else
	for i = 1:length(files)
        s = files{i};
        num(i) = str2num(s((strfind(s, '#')+1):(strfind(s, '.mat')-1)));
    end
    files = cellstr(string([min(num):max(num)]'));
%     files = num2cell((min(num):max(num))');
    files(2:end+1) = files(1:end);
    files(1) = {'Please select number of trial'};
    set(handles.select_data, 'string', files)
end


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  structure with the following fields
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
files = getappdata(handles.browse_botton, 'files');
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton1'
        files(2:end+1) = files(1:end);
        files(1) = {'Please select data'};
        set(handles.select_data, 'string', files)
        set(handles.uipanel1,'visible','off')
        set(handles.uipanel2,'visible','on')
    case 'radiobutton2'
        for i = 1:length(files)
            s = files{i};
            num(i) = str2num(s((strfind(s, '#')+1):(strfind(s, '.mat')-1)));
        end
        files = cellstr(string([min(num):max(num)]'));
%     files = num2cell((min(num):max(num))');
        files(2:end+1) = files(1:end);
        files(1) = {'Please select number of trial'};
        set(handles.select_data, 'string', files)
        data = 0;
        setappdata(handles.select_data, 'data', data);
        set(handles.select_data, 'Value', 1)
        set(handles.uipanel2,'visible','off')
        set(handles.uipanel1,'visible','on')
        set(handles.axes3, 'visible', 'off')
        set(handles.text_mean2,'visible','off')
        set(handles.text_std2,'visible','off')
end



function n_bins_Callback(hObject, eventdata, handles)
% hObject    handle to n_bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_bins as text
%        str2double(get(hObject,'String')) returns contents of n_bins as a double


% --- Executes during object creation, after setting all properties.
function n_bins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1



function sigma_Callback(hObject, eventdata, handles)
% hObject    handle to sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma as text
%        str2double(get(hObject,'String')) returns contents of sigma as a double


% --- Executes during object creation, after setting all properties.
function sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2

if get(handles.checkbox2, 'Value') == 1.0
    set(handles.uipanel3, 'visible', 'on')
else
    set(handles.uipanel3, 'visible', 'off')
end


% --- Executes when selected object is changed in uibuttongroup2.
function uibuttongroup2_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup2 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton3'
        calculate_mean(hObject, eventdata, handles)
        set(handles.axes3,'visible','off')
        set(handles.text_mean2,'visible','off')
        set(handles.text_std2,'visible','off')
    case 'radiobutton4'
        calculate_mean(hObject, eventdata, handles)
        if get(handles.checkbox2, 'Value') == 1.0
            set(handles.uipanel3, 'visible', 'on')
            set(handles.axes3,'visible','on')
            plot_hist(hObject, eventdata, handles)
        end
end
