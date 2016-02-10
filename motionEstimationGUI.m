function varargout = motionEstimationGUI(varargin)
% MOTIONESTIMATIONGUI MATLAB code for motionEstimationGUI.fig
%      MOTIONESTIMATIONGUI, by itself, creates a new MOTIONESTIMATIONGUI or raises the existing
%      singleton*.
%
%      H = MOTIONESTIMATIONGUI returns the handle to a new MOTIONESTIMATIONGUI or the handle to
%      the existing singleton*.
%
%      MOTIONESTIMATIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOTIONESTIMATIONGUI.M with the given input arguments.
%
%      MOTIONESTIMATIONGUI('Property','Value',...) creates a new MOTIONESTIMATIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before motionEstimationGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to motionEstimationGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help motionEstimationGUI

% Last Modified by GUIDE v2.5 01-Sep-2015 11:11:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @motionEstimationGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @motionEstimationGUI_OutputFcn, ...
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


% --- Executes just before motionEstimationGUI is made visible.
function motionEstimationGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to motionEstimationGUI (see VARARGIN)

% Choose default command line output for motionEstimationGUI
handles.output = hObject;

[hObject,handles] = hideNorm(hObject, eventdata, handles);
addpath('tools');
addpath(['algorithms',filesep,'matlab']);
addpath(['algorithms',filesep,'mex']);
addpath('operators');

    set(handles.axisRawImage1,'xtick',[]);
    set(handles.axisRawImage1,'ytick',[]);
    set(handles.axisRawImage2,'xtick',[]);
    set(handles.axisRawImage2,'ytick',[]);
    set(handles.axisResultImage,'xtick',[]);
    set(handles.axisResultImage,'ytick',[]);
    set(handles.axisFlowMagnitude,'xtick',[]);
    set(handles.axisFlowMagnitude,'ytick',[]);
    

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes motionEstimationGUI wait for user response (see UIRESUME)
% uiwait(handles.mainFigure);


% --- Outputs from this function are returned to the command line.
function varargout = motionEstimationGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in buttonLoadImages.
function buttonLoadImages_Callback(hObject, eventdata, handles)
% hObject    handle to buttonLoadImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if (exist('handles.pfad','var'))
        [FileName,PathName] = uigetfile([handles.pfad,'*.*'],'Load First Image');
    else
        [FileName,PathName] = uigetfile('*.*','Load First Image');
    end
    
    if (~FileName)
        return;
    end
    
    handles.pfad = PathName;
    
    [FileName2,PathName2] = uigetfile([handles.pfad,'*.*'],'Load Second Image');
    
    if (~FileName2)
        return;
    end
    
    image1 = imread([PathName,FileName]);
    image2 = imread([PathName2,FileName2]);
    
    if (ndims(image1)==3)
        image1 = rgb2gray(image1);
        image2 = rgb2gray(image2);
    end
    
    handles.rawImage1 = im2double(image1);
    handles.rawImage2 = im2double(image2);
    
 
    axes(handles.axisRawImage1);
    imagesc(handles.rawImage1);colormap(gray);
    set(handles.axisRawImage1,'xtick',[]);set(handles.axisRawImage1,'ytick',[]);
    axis image;
    
    axes(handles.axisRawImage2);
    imagesc(handles.rawImage2);colormap(gray);
    set(handles.axisRawImage2,'xtick',[]);set(handles.axisRawImage2,'ytick',[]);
    axis image;
guidata(hObject, handles);


% --- Executes on selection change in selectAlgorithm.
function selectAlgorithm_Callback(hObject, eventdata, handles)
    algorithmNumber = get(handles.selectAlgorithm,'Value');
    switch algorithmNumber
        case {1,2,7,8}
            [hObject,handles] = hideNorm(hObject, eventdata, handles);
        otherwise
            [hObject,handles] = showNorm(hObject, eventdata, handles);
    end
    switch algorithmNumber
        case {1,5,6,7,8,9}
            set(handles.parameterAlpha,'String','0.1');
        case {2}
            set(handles.parameterAlpha,'String','1');
        case {3,4}
            set(handles.parameterAlpha,'String','0.01');
    end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function selectAlgorithm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectAlgorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function parameterAlpha_Callback(hObject, eventdata, handles)
% hObject    handle to parameterAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of parameterAlpha as text
%        str2double(get(hObject,'String')) returns contents of parameterAlpha as a double


% --- Executes during object creation, after setting all properties.
function parameterAlpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameterAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function parameterTol_Callback(hObject, eventdata, handles)
% hObject    handle to parameterTol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of parameterTol as text
%        str2double(get(hObject,'String')) returns contents of parameterTol as a double


% --- Executes during object creation, after setting all properties.
function parameterTol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameterTol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function parameterMaxIterations_Callback(hObject, eventdata, handles)
% hObject    handle to parameterMaxIterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of parameterMaxIterations as text
%        str2double(get(hObject,'String')) returns contents of parameterMaxIterations as a double


% --- Executes during object creation, after setting all properties.
function parameterMaxIterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameterMaxIterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in parameterNorm.
function parameterNorm_Callback(hObject, eventdata, handles)
% hObject    handle to parameterNorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns parameterNorm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from parameterNorm


% --- Executes during object creation, after setting all properties.
function parameterNorm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameterNorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [hObject,handles] = hideNorm(hObject, eventdata, handles)
    set(handles.textNorm,'Visible','Off');
    set(handles.parameterNorm,'Visible','Off');
guidata(hObject, handles);

function [hObject,handles] = showNorm(hObject, eventdata, handles)
    set(handles.textNorm,'Visible','On');
    set(handles.parameterNorm,'Visible','On');
guidata(hObject, handles);


% --- Executes on button press in buttonRun.
function buttonRun_Callback(hObject, eventdata, handles)
% hObject    handle to buttonRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

alpha = str2double(get(handles.parameterAlpha,'String'));
tol = str2double(get(handles.parameterTol,'String'));
maxIterations = str2double(get(handles.parameterMaxIterations,'String'));
numSteps = str2double(get(handles.parameterSteps,'String'));
norm = get(handles.parameterNorm,'Value');
if (norm == 1)
    norm = 3;
else
    norm = 1;
end

if (isfield(handles, 'rawImage1') && isfield(handles, 'rawImage2'))
    set(handles.textStatus,'Visible','On','String','Calculating, please wait...');drawnow
    
    algorithmNumber = get(handles.selectAlgorithm,'Value');
    
    u = cat(3,handles.rawImage1,handles.rawImage2);
    dimsU = size(u);
    
    switch algorithmNumber
        case 1
            [x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L2L2OpticalFlow',4,'maxIt',maxIterations,'numsteps',numSteps);
            v1 = x(:,:,1,1);
            v2 = x(:,:,1,2);
        case 2
            [x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L2L2MassPreservation',5,'maxIt',maxIterations,'numsteps',numSteps);
            v1 = x(:,:,1,1);
            v2 = x(:,:,1,2);
        case 3
            [x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L2TVOpticalFlow',4,'maxIt',maxIterations,'numsteps',numSteps,'typeNorm',norm);
            v1 = x(:,:,1,1);
            v2 = x(:,:,1,2);
        case 4
            [x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L2TVMassPreservation',5,'maxIt',maxIterations,'numsteps',numSteps,'typeNorm',norm);
            v1 = x(:,:,1,1);
            v2 = x(:,:,1,2);
        case 5
            [x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L1TVOpticalFlow',4,'maxIt',maxIterations,'numsteps',numSteps,'typeNorm',norm);
            v1 = x(:,:,1,1);
            v2 = x(:,:,1,2);
        case 6
            [x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L1TVMassPreservation',5,'maxIt',maxIterations,'numsteps',numSteps,'typeNorm',norm);
            v1 = x(:,:,1);
            v2 = x(:,:,2);
        case 7
            [x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L1L2OpticalFlow',4,'maxIt',maxIterations,'numsteps',numSteps);
            v1 = x(:,:,1,1);
            v2 = x(:,:,1,2);
        case 8
            [x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L1L2MassPreservation',5,'maxIt',maxIterations,'numsteps',numSteps);
            v1 = x(:,:,1);
            v2 = x(:,:,2);
        case 9
            [x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L1TVOpticalFlowNonlinear',4,'maxIt',maxIterations,'numsteps',numSteps,'typeNorm',norm);
            v1 = x(:,:,1,1);
            v2 = x(:,:,1,2);
    end
    
    axes(handles.axisResultImage);
    imagesc(flowToColorV2(cat(3,v1,v2)));
    set(handles.axisResultImage,'xtick',[]);
    set(handles.axisResultImage,'ytick',[]);
    axis image;
    
    axes(handles.axisFlowMagnitude);
    imagesc(sqrt(v1.^2 + v2.^2));
    set(handles.axisFlowMagnitude,'xtick',[]);
    set(handles.axisFlowMagnitude,'ytick',[]);
    axis image;colorbar;
    
    handles.resultV1 = v1;
    handles.resultV2 = v2;
    
    set(handles.textStatus,'Visible','On','String','Finished');drawnow
else
    errordlg('Please load images first')
end
guidata(hObject, handles);

% --- Executes on selection change in exportFormat.
function exportFormat_Callback(hObject, eventdata, handles)
% hObject    handle to exportFormat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns exportFormat contents as cell array
%        contents{get(hObject,'Value')} returns selected item from exportFormat


% --- Executes during object creation, after setting all properties.
function exportFormat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to exportFormat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in buttonExport.
function buttonExport_Callback(hObject, eventdata, handles)

if (isfield(handles, 'resultV1'))
    exportFormat = get(handles.exportFormat,'Value');
    
    switch exportFormat
        case 1
            [FileName,PathName] = uiputfile('*.flo','Please specify export folder');

            if (~PathName)
                return;
            end
            writeFlowFile(cat(3,handles.resultV1,handles.resultV2), [PathName,FileName]);
        case 2
            [FileName,PathName, filterindex] = uiputfile('*.txt','Please specify export folder');
            if (~PathName)
                return;
            end
            
            FileName = FileName(1:end-4);
            
            v1 = handles.resultV1;
            v2 = handles.resultV2;
            save([PathName,FileName,'v1.txt'],'v1','-ASCII');
            save([PathName,FileName,'v2.txt'],'v2','-ASCII');
    end

else
    errordlg('Please calculate flow field first');
end



function parameterSteps_Callback(hObject, eventdata, handles)
% hObject    handle to parameterSteps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of parameterSteps as text
%        str2double(get(hObject,'String')) returns contents of parameterSteps as a double


% --- Executes during object creation, after setting all properties.
function parameterSteps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameterSteps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
