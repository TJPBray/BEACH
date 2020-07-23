function varargout = BEACH_FF_R2(varargin)
%% Dr Tim Bray
% t.bray@ucl.ac.uk
%
% This code implements the tool for ADC measurement described in:
% Histographic analysis of oedema and fat in inflamed bone marrow based on qMRI.
% Bray et al. European Radiology April 2020
%
% The code is based on the earlier 'roianal' software with additional
% functionality to enable 
% (1) semi-automated ROI placement once the observer has defined the joint
% (2) histographic analysis of voxel values in the ROI
% (3) in BEACH_FF, the ROIs are automatically propagated onto R2* maps  
% (acquired together with the PDFF maps) to enable 2D histogram / density
% analysis - this version therefore requires two inputs
%
%% 
% BEACH_FF_R2(FFvol,R2vol)

% BEACH_FF_R2(FFvol,R2vol, matp, 'parameter',value,...)
%
% [handles, roi] = BEACH_FF_R2(...) 
%
% outp is a structure from d2mat with outp.geom.XData and YData, or
% outp.XData and YData
%
% Parameter value pairs allowed:
%  'ismodal'  true | {false}
%  'Name' 
%
%      BEACH_FF_R2, by itself, creates a new BEACH_FF_R2 or raises the existing
%      singleton*.
%
%      H = BEACH_FF_R2 returns the handle to a new BEACH_FF_R2 or the handle to
%      the existing singleton*.
%
%      BEACH_FF_R2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEACH_FF_R2 with the given input arguments.
%
%      BEACH_FF_R2('Property','Value',...) creates a new BEACH_FF_R2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before roianal_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to roianal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES, profile_anal

% Edit the above text to modify the response to help roianal

% Last Modified by GUIDE v2.5 17-Jul-2020 16:22:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @roianal_OpeningFcn, ...
                   'gui_OutputFcn',  @roianal_OutputFcn, ...
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


% --- Executes just before roianal is made visible.
function roianal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to roianal (see VARARGIN)

% Choose default command line output for roianal
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes roianal wait for user response (see UIRESUME)
% uiwait(handles.figure1);

hf = figure('Name','roianalfig') ;
set(hf,'Tag','roianalfig')
disp(['Opened fig: ',num2str(fignum(hf))])


vol = varargin{1} ;
vol2 = varargin{2};
if size(vol,3) == 1
    vol = squeeze(vol) ;
    vol2=squeeze(vol2) ;
end
nslice = size(vol,3) ;

%initial slice
slice1 = floor(nslice/2) + 1;

%set slider and edit box
if nslice > 1
  set(handles.slider_slice,'Min',1,'Max',nslice,'SliderStep',[1/(nslice-1) max(0.1,1/(nslice-1)) ],...
      'Value', slice1) ;
else
    set(handles.slider_slice,'Visible','off')
end
set(handles.edit_slice,'String',num2str(slice1)) ;

if length(varargin) < 2 || ~isstruct(varargin{2})
    XData(1) = 1 ; YData(1) = 1 ;
    XData(2) = XData(1) + size(vol,2)-1 ;
    YData(2) = YData(1) + size(vol,1)-1 ;
else
    if isfield(varargin{2},'geom')
      geom = varargin{2}.geom ;
    else
        geom = varargin{2} ;
    end
    XData = geom.XData ;
    YData = geom.YData ;
    
    if isfield(varargin{2},'TinSeriesVec')
        UDfig.TinSeries = varargin{2}.TinSeriesVec/1e6 ;
    end
end

%warm-up axes
himshow = imshow(vol(:,:,slice1),[],'XData',XData,'YData',YData) ;
ha = gca ;
set(himshow,'Visible','off')

%apply all slices to axes, but make in visible
ip = {'XData',XData,'YData',YData,'CDataMapping','scaled','Visible','off','CData'} ;

for islice = 1:size(vol,3)
    him(islice) = image(ip{:},vol(:,:,islice)) ;
end

ismodal= false ;

if nargin > 2
    for ip =3 : 2 : length(varargin)
        switch varargin{ip}
            case 'ismodal'
                if varargin{ip+1} ~= false 
                    ismodal = true ;
                    set(handles.pushbutton_finish,'Visible','on')
                end
            case {'Name','name'}
                Name = varargin{ip+1} ;
                set(hf,'Name',Name)
            otherwise
                warning(['Unknown parameter: ',varargin{ip}])
        end
    end
end

if ~exist('buildroianal.m','file')
    set(handles.pushbutton_report,'Enable','off')
end

UDcont.ismodal = ismodal ;
UDfig.nslice = nslice ;
UDfig.him = him ;
UDcont.hf = hf ; % control panel records active figure
UDfig.ha = ha ;

UDfig.vol = vol ;
UDfig.vol2= vol2 ;


set(handles.figure1,'UserData',UDcont) ;
set(hf,'UserData',UDfig) ;

draw_slice(handles) ;
update_activefig(handles) ;
turn_off_link(handles) ;
update_listbox(handles)

if ismodal
    uiwait(handles.figure1);
end


%-------------------------------------

function draw_slice(handles)
slicestr = get(handles.edit_slice,'String') ;
islice = str2num(slicestr) ;

UDcont = get(handles.figure1,'UserData') ;
hf = UDcont.hf ;
UD = get(hf,'UserData') ;
draw_slice_one(UD,islice,hf) ;
islinked = get(handles.checkbox1,'Value');
if islinked
    val = get(handles.popupmenu_linkfig,'Value') ;
    str = get(handles.popupmenu_linkfig,'String') ;

    fig = sscanf(str{val},'Fig: %d%s') ;
    hflink = fig(1) ;
    
    UDlink = get(hflink,'UserData') ;
    draw_slice_one(UDlink,islice,hflink) ;
end

function roi_on_slice = draw_slice_one(UD,islice,hf)
set(UD.him,'Visible','off') 
set(UD.him(islice),'Visible','on')
UD.disp_slice = islice ;
set(hf,'UserData',UD) ;
roi_on_slice = false ;

if isfield(UD,'hroi')
    hrois = [UD.hroi] ;
    for iroi = 1: length(hrois)
        if ~isempty(hrois{iroi})
            sl = UD.roislice(iroi) ;
            if islice == sl
                set(hrois{iroi},'Visible','on')
                roi_on_slice = true ;
            else
                set(hrois{iroi},'Visible','off')
            end
        end
    end
end

function draw_montage(handles)
% DRAW_MONTAGE Compiles montage from slices with ROI drawn
slicestr = get(handles.edit_slice,'String') ;
curr_slice = str2num(slicestr) ;

UDcont = get(handles.figure1,'UserData') ;
hf = UDcont.hf ;
UD = get(hf,'UserData') ;

nslice = size(UD.vol , 3) ;

iframe = 0 ;
M = [] ;
for islice = 1:nslice
    roi_on_slice = draw_slice_one(UD,islice,hf) ;
    if roi_on_slice
        iframe = iframe + 1 ;
        f = getframe(UD.ha) ;
        [im,map] = frame2im(f);    %Return associated image data 
        if isempty(map)            %Truecolor system
          rgb = im;
        else                       %Indexed system
          rgb = ind2rgb(im,map);   %Convert image data
        end
         M = cat(4,M,rgb) ;
    end
end

draw_slice_one(UD,curr_slice,hf) ;


if iframe > 0
    figure('Name','roianal_montage','Tag','roianal_montage')
    montage(M)
end



% --- Outputs from this function are returned to the command line.
function varargout = roianal_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


varargout{1} = handles.output;
if nargout > 1
    UDcont = get(handles.figure1,'UserData') ;
    UD = get(UDcont.hf,'Userdata') ;
    
    lstr = get(handles.listbox_rois,'String') ;
    lval = get(handles.listbox_rois,'Value') ;
    
    nsel = length(lval) ;
    rois = zeros([1 nsel]) ;
    
    for isel = 1:nsel
        sstr = lstr{lval(isel)} ;
        rois(isel) = strread(sstr(1:4),'%u') ;
    end
    
    posc = cell([1 nsel]) ;
    slice = zeros([1 nsel]) ;
    BWc = cell([1 nsel]) ;
    labelc = cell([1 nsel]) ;
    
    for isel = 1:nsel
        pos = getPosition(UD.hroi{rois(isel)}) ;
        posc{isel} = pos ;
        slice(isel) = UD.roislice(rois(isel)) ;
        BW = createMask(UD.hroi{rois(isel)}, UD.him(slice(isel))) ;
        BWc{isel} = BW ;
        labelc{isel} = UD.roilabel(rois(isel));
    end
    
    roi.posc = posc ;
    roi.slice = slice ;
    roi.BWc = BWc ;
    roi.labelc = labelc ;
    
    
    varargout{2} = roi ;
end


% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(handles.figure1)


% --- Executes on slider movement.
function slider_slice_Callback(hObject, eventdata, handles)
% hObject    handle to slider_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

islice = get(hObject,'Value') ;
islice = round(islice) ;
set(hObject,'Value',islice) ;

set(handles.edit_slice,'String',num2str(islice) ) ;

draw_slice(handles) ;


% --- Executes during object creation, after setting all properties.
function slider_slice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_slice_Callback(hObject, eventdata, handles)
% hObject    handle to edit_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_slice as text
%        str2double(get(hObject,'String')) returns contents of edit_slice as a double

islice = str2double(get(hObject,'String')) ;
set(handles.slider_slice,'Value',islice)

draw_slice(handles) ;


% --- Executes during object creation, after setting all properties.
function edit_slice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_imcontrast.
function pushbutton_imcontrast_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_imcontrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


UDcont = get(handles.figure1,'UserData') ;
islice = str2double(get(handles.edit_slice,'String')) ;
UD = get(UDcont.hf,'UserData') ;
imcontrast(UD.him(islice)) 

% --- Executes on selection change in listbox_rois.
function listbox_rois_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_rois contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_rois


% --- Executes during object creation, after setting all properties.
function listbox_rois_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'Max',2) ; % allow multiple selection


% --- Executes on button press in pushbutton_newroi.
function pushbutton_newroi_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_newroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UDcont = get(handles.figure1,'UserData') ;
hf = UDcont.hf ;
UD = get(hf,'UserData') ; %UD contains stack of images, number of slices, current slice, and some scaling info

[hroi_1,hroi_2] = draw_imroi(handles,UD.ha) ; %hroi_this is a polygon

if isfield(UD,'hroi')
    hrois = UD.hroi ;
    nroi = length(hrois) ;
else
    nroi = 0 ;
end

islice = str2double(get(handles.edit_slice,'String')) ; %islice is the slice for the roi?
labstr = get(handles.edit_roilabel,'String') ; %labstr is label

UD.roislice(nroi+1) = islice ;
UD.roislice(nroi+2) = islice ;
UD.hroi{nroi+1} = hroi_1 ; %UD.hroi{} is a list of polygons
UD.hroi{nroi+2} = hroi_2 ;
UD.roilabel{nroi+1} = strcat(labstr,' Left');
UD.roilabel{nroi+2} = strcat(labstr, ' Right') ;

set(UDcont.hf,'UserData',UD)
update_listbox(handles, nroi+1)

%------------------------------------------
function [hroi,hroi2] = draw_imroi(handles, varargin)

if length(varargin)==3
    roitype = varargin{3} ;
    vin{1}  = varargin{1};
    vin{2}  = varargin{2};
else
    lval = get(handles.popupmenu_roitype,'Value') ;
    lstr = get(handles.popupmenu_roitype,'String') ;
    roitype = lstr{lval} ;
    vin = varargin ;
end
switch roitype
    
    
    case 'impoly'
        
        hroi = impoly(vin{:},'Closed',false) ;
        setColor(hroi,'red');
        posn = getPosition(hroi) ;
        if size(posn,1) < 3
           warning(['Less than 3 points in impoly'])
           delete(hroi)
           hroi = draw_imroi(handles, varargin{:}) ;
        end
        
      
    get(handles.edit5)    
  width=handles.edit5.Value
  get(handles.edit8)
   widthinner=handles.edit8.Value    
   
im=imgcf
im=double(im)
polypoints=hroi.getPosition
delete(hroi)

%Get wide and narrow polygons
[xy1wide,xy2wide]=create2poly(polypoints,width,im);
[xy1narrow,xy2narrow]=create2poly(polypoints,widthinner,im);

%Create left and right polygons
xyL=[xy1wide; flipud(xy1narrow)];
xyR=[xy2wide; flipud(xy2narrow)];

hroi=impoly(gca,xyL,'Closed',true);
hroi2=impoly(gca,xyR,'Closed',true);
setColor(hroi,'red');
setColor(hroi2,'red');

        set(hroi,'Tag','impoly')
        set(hroi2,'Tag','impoly')
   
end


function hroi = draw_normroi(handles, varargin)

if length(varargin)==3
    roitype = varargin{3} ;
    vin{1}  = varargin{1};
    vin{2}  = varargin{2};
else
    lval = get(handles.popupmenu_roitype,'Value') ;
    lstr = get(handles.popupmenu_roitype,'String') ;
    roitype = lstr{lval} ;
    vin = varargin ;
end

        hroi = impoly(vin{:}) ;
        posn = getPosition(hroi) ;
        if size(posn,1) < 3
           warning(['Less than 3 points in impoly'])
           delete(hroi)
           hroi = draw_normroi(handles, varargin{:}) ;
        end
        set(hroi,'Tag','normroi')



%------------------------------------------
function update_listbox(handles, varargin)

if length(varargin) > 0
    roi_selected = varargin{1} ;
else
    roi_selected = 0 ;
end

UDcont = get(handles.figure1,'UserData') ;
UD = get(UDcont.hf,'UserData') ;
if ~isfield(UD,'hroi')
    lstr{1} = 'ROI information';
    lbval = 1 ;
else
    hrois = [UD.hroi] ;
    ilstr = 0 ;
    lbline = 0 ;
    lbval = [] ;
    lstr = {'ROI information'} ;
    nroi = length(hrois);
    
    for iroi = 1: nroi
        if ~isempty(hrois{iroi})
            ilstr = ilstr+1;
            sl = UD.roislice(iroi) ;
            labstr = UD.roilabel{iroi};
            roitype = get(hrois{iroi},'Tag') ;
            switch roitype
                 case 'impoly'
                    BW = createMask(hrois{iroi},UD.him(sl));
                  stats = regionprops(BW,UD.vol(:,:,sl), ...
                 'Area','MeanIntensity','MaxIntensity','MinIntensity') ;
             
             tempim=double(UD.vol(:,:,sl));
             masked=BW.*tempim;
             masked(masked==0)=NaN;
             values=(masked(~isnan(masked)));
             p95 = prctile(values,95);
                                                     
                  dstr = [ '. IMPOLY. Area ',num2str(stats.Area), ...
                '. Mean ',num2str(stats.MeanIntensity), ...
                ', Min ',num2str(stats.MinIntensity), ...
                ', Max ',num2str(stats.MaxIntensity), ...
                ', P95 ',num2str(p95)];
           
               
            
                case 'normroi'
                    BW = createMask(hrois{iroi},UD.him(sl));
                  stats = regionprops(BW,UD.vol(:,:,sl), ...
                 'Area','MeanIntensity','MaxIntensity','MinIntensity') ;
             
             tempim=double(UD.vol(:,:,sl));
             masked=BW.*tempim;
             masked(masked==0)=NaN;
             values=(masked(~isnan(masked)));
             p95 = prctile(values,95);
                                                     
                  dstr = [ '. IMPOLY. Area ',num2str(stats.Area), ...
                '. Mean ',num2str(stats.MeanIntensity), ...
                ', Min ',num2str(stats.MinIntensity), ...
                ', Max ',num2str(stats.MaxIntensity), ...
                ', P95 ',num2str(p95)];
                    
                    
            %Set p95 for use by histogram button
            get(handles.edit9);
            handles.edit9.Value=p95;
            
                otherwise
                    warning(['Unknown imroi type'])
            end
            
            
            lstr{ilstr} = [num2str(iroi,'%04u'), ...
                '. Slice ',num2str(sl,'%04u'), ...
                '. <',labstr,'>', ...
                dstr ] ;
            
            if ismember(iroi,roi_selected)
                lbval = [ lbval ilstr] ;
            end
        end
        
    end
end
set(handles.listbox_rois,'String',lstr)
set(handles.listbox_rois,'Value',lbval)

% --- Executes on button press in pushbutton_roidelete.
function pushbutton_roidelete_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_roidelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% removes "deleted" from list (doesn't actually delete object)
rois = get_selected_rois(handles) ;
UDcont = get(handles.figure1,'UserData') ;
UD = get(UDcont.hf,'UserData') ;
hrois = UD.hroi ;
for iroi = 1:length(rois)
  set(hrois{rois(iroi)},'Visible','off')
  hrois{rois(iroi)} = [] ;
end

UD.hroi = hrois ;
set(UDcont.hf,'UserData',UD) ;
update_listbox(handles)


% --- Executes on button press in pushbutton_roisave.
function pushbutton_roisave_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_roisave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%loop over all rois in listbox, save pos, XData, YData, vol
%size
UDcont = get(handles.figure1,'UserData') ;
UD = get(UDcont.hf,'UserData') ;

% update_listbox(handles)
lstr = get(handles.listbox_rois,'String') ;
%  lval = [1:length(lstr)]; % fakes all selected in listbox
lval = get(handles.listbox_rois,'Value') ;
if length(lval)<length(lstr)
    warndlg(['Selected ',num2str(length(lval)),' out of ',...
        num2str(length(lstr)),' ROIs.'],'Warning on Save','replace')
end

nsel = length(lval) ; 
rois = zeros([1 nsel]) ;

for isel = 1:nsel
    sstr = lstr{lval(isel)} ;
    rois(isel) = strread(sstr(1:4),'%u') ;
end

posc = cell([1 nsel]) ;
labstrc = cell([1 nsel]);
roitype = cell([1 nsel]);
slice = zeros([1 nsel]) ;
roistats = cell([1 nsel]);
BW = cell([1 nsel]) ;

hrois = [UD.hroi];

for isel = 1:nsel
    pos = getPosition(UD.hroi{rois(isel)}) ;
    posc{isel} = pos ;
    labstrc{isel} = UD.roilabel(rois(isel)) ;
    slice(isel) = UD.roislice(rois(isel)) ;
    roitype{isel} = get(hrois{rois(isel)},'Tag') ;
    roistats{isel} = calc_roistats(handles,roitype{isel}, pos, slice(isel)) ;
    BW{isel} = createMask(hrois{rois(isel)}, UD.him(slice(isel)) ) ;
end

if ispref('roianal','save_dir')
    defdir = getpref('roianal','save_dir');
else
    defdir = [] ;
end

% in future could also save in OsiriX format
[fn,pn,FilterIndex] = uiputfile({'*.mat','MATLAB file (*.mat)'; ...
                                 '*.txt','Text file (*.txt)' },'Enter filename',defdir) ;

if fn == 0
    warning(['No file specified'])
    return
end

if FilterIndex == 1  % *.mat
    sz_vol = size(UD.vol) ;
    XData = get(UD.him(1),'XData') ;
    YData = get(UD.him(1),'YData') ;
    
    save(fullfile(pn,fn),'posc','slice','sz_vol','XData','YData',...
        'labstrc','roitype','roistats','BW') ;
elseif FilterIndex==2  % *.txt
    [fid, msg] = fopen(fullfile(pn,fn),'w') ;
    if fid == -1 
        error(['File opening error: ',msg])
    end
    fprintf(fid,'%s\n',datestr(now)) ;
    
    for isel = 1:nsel
      fprintf(fid,'%s Label: %s,  Slice: %s Profile\n',roitype{isel}, ...
                               labstrc{isel}{1}, num2str(slice(isel)) ) ;
      fprintf(fid,'Distance , Profile \n')                     
      if ~isempty(roistats{isel})
        prof = roistats{isel}.profile ;
        dist = roistats{isel}.dist ;
        for ip = 1:length(prof)
            fprintf(fid,'%f ,  %f\n',dist(ip), prof(ip)) ;
        end
      end
    
    fprintf(fid,'\n') ;
    
    fprintf(fid,'%s Label: %s,  Slice: %s Area\n',roitype{isel}, ...
                               labstrc{isel}{1}, num2str(slice(isel)) ) ;
    fprintf(fid,'Threshold , Area \n')  
      if ~isempty(roistats{isel})
        A = roistats{isel}.areas ;
        ths = roistats{isel}.threshs ;
        
        for ip = 1:length(A)
            fprintf(fid,'%f ,  %f\n',ths(ip),A(ip)) ;
        end
      end
    
    fprintf(fid,'\n') ;
    
    end
    
    fclose(fid) ;
    disp(['Written to file: ',fullfile(pn,fn)])
else
    warning('FilterIndex ??')
end
setpref('roianal','save_dir',pn)


% --- Executes on button press in pushbutton_roiload.
function pushbutton_roiload_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_roiload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ispref('roianal','load_dir')
    defdir = getpref('roianal','load_dir');
else
    defdir = [] ;
end

% in future could also load in OsiriX format
[fn,pn,FilterIndex] = uigetfile({'*.mat','MATLAB mat file'},'Select ROI filename',defdir) ;

if fn == 0
    warning(['No roi file specified'])
    return
end

setpref('roianal','load_dir',pn) ;

UDcont = get(handles.figure1,'UserData') ;
hf = UDcont.hf ;
UD = get(hf,'UserData') ;

sz_vol = size(UD.vol) ;
XData = get(UD.him(1),'XData') ;
YData = get(UD.him(1),'YData') ;

S = load(fullfile(pn,fn)) ;

if ~isequal(S.sz_vol,sz_vol)
    warning(['Input ROIs created on different sized volume'])
    if S.sz_vol(3) ~= sz_vol(3)
        return
    end
end
if ~isequal(S.XData,XData)
    warning(['Input ROIs created on differently scaled image in X'])
end
if ~isequal(S.YData,YData)
    warning(['Input ROIs created on differently scaled image in Y'])
end

currsl = str2num(get(handles.edit_slice,'String')) ;

if isfield(UD,'hroi')
    hrois = [UD.hroi] ;
    nrois = length(hrois) ;
else
    nrois = 0 ;
end

nroi_load = length(S.posc) ;

for iroi = 1:nroi_load
    h_this = draw_imroi(handles,UD.ha, S.posc{iroi}, S.roitype{iroi}) ;
    if S.slice(iroi) == currsl
        set(h_this,'Visible','on')
    else
        set(h_this,'Visible','off')
    end
    
    nrois = nrois + 1 ;
    UD.hroi{nrois} = h_this ;
    UD.roislice(nrois) = S.slice(iroi) ;
    UD.roilabel(nrois) = S.labstrc{iroi} ;
    
end
set(hf,'UserData',UD) ;

draw_slice(handles)
update_listbox(handles)
        

% --- Executes on button press in pushbutton_histogram.
function pushbutton_histogram_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_histogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rois = get_selected_rois(handles) ;

nroi = length(rois) ;
if nroi == 0
    disp(['No ROIs selected'])
    return
end

UDcont = get(handles.figure1,'UserData') ;
UD = get(UDcont.hf,'UserData') ;

hrois = [UD.hroi] ;

pv = [] ;
pv2 = [];
for iroi = 1: nroi
    roi_this = rois(iroi) ;
    if ~isempty(hrois{roi_this})
        slice_this = UD.roislice(roi_this) ;
        BW = createMask(hrois{roi_this},UD.him(slice_this)) ;
        %Get stats for vol 1
        stats = regionprops(BW,UD.vol(:,:,slice_this), ... 
            'PixelValues') ;
        statinf=stats.PixelValues;
        pv = [pv ; statinf] ;
        %Get stats for vol 2
                stats2 = regionprops(BW,UD.vol2(:,:,slice_this), ... 
            'PixelValues') ;
        pv2 = [pv2 ; stats2.PixelValues] ;
    else
        warning(['Empty ROI: ', num2str(roi_this)])
    end
end

%Calculate key FF values
FFmean=round(nanmean(pv),2);
FF10=round(prctile(pv,10),2);
FF25=round(prctile(pv,25),2);
FFmedian=round(prctile(pv,50),2);
FF75=round(prctile(pv,75),2);
FF90=round(prctile(pv,90),2);
pvmax=round(max(pv),2);


%Percentage above thresholds 
thresh1=69.7;
totalpixels=numel(pv);
highFFpixels=sum(pv>thresh1); 
highFFperc=highFFpixels/totalpixels;

thresh2=28.2;
lowFFpixels=sum(pv<thresh2); 
lowFFperc=lowFFpixels/totalpixels;

%Calculate key R2* values
R2mean=round(nanmean(pv2),3);

R2_10=round(prctile(pv2,10),3);
R2_25=round(prctile(pv2,25),3);
R2median=round(prctile(pv2,50),3);
R2_75=round(prctile(pv2,75),3);
R2_90=round(prctile(pv2,90),3);
R2max=round(max(pv2),3);

%Create matrix with data for easy export
HistData(1,1)=FFmean;
HistData(1,2)=FFmedian;
HistData(1,3)=FF10;
HistData(1,4)=FF25;
HistData(1,5)=FF75;
HistData(1,6)=FF90;
HistData(1,7)=highFFperc;
HistData(1,8)=lowFFperc;
HistData(1,9)=R2mean;
HistData(1,10)=R2median;
HistData(1,11)=R2_10;
HistData(1,12)=R2_25;
HistData(1,13)=R2_75;
HistData(1,14)=R2_90;

HistData=HistData.'

%Create column names
Param={'FFmean';'FFmedian';'FF10';'FF25';'FF75';'FF90';'HighFF';'LowFF';'R2mean';'R2median';'R2_10';'R2_25';'R2_75';'R2_90'};

T=table(HistData,'RowNames',Param)

%Calculate fixed free diffusion fraction
%diffthresh=957;
%totalpixels=numel(pv);
%abovethreshpixels=sum(pv>diffthresh); %957 is defined threshold for free diffusion fraction
%FDFfixed=abovethreshpixels/totalpixels;

%Calculate variable free diffusion fraction (above P95 normal)
   %get(handles.edit9);
   %N95=handles.edit9.Value;

%abovethreshpixels2=sum(pv>N95);
%FDFvar=abovethreshpixels2/totalpixels;


%%%% CREATE FIGURES
figure ('Name',['Histogram ROI # ',num2str(rois)], 'Tag', 'roianal_histogram')

% (1) FF histogram
subplot(2,2,1)

    %Specific bin widths
    xbins=-10:5:110;

hist(double(pv),xbins, 'FontSize', 14)
xlim([0 100]);
hhist = gca ;
hold on
title(['Fat Fraction Histogram'])
xlabel(['FF'])
ylabel(['Number of voxels'])
set(gca,'FontSize',16);

%Plot key values
yl=ylim;
plot([FFmedian FFmedian],[yl(1) yl(2)],'LineWidth',2,'color','red','Linestyle','--') %median
plot([FF10 FF10],[yl(1) yl(2)],'LineWidth',2,'color','red','Linestyle','--') %FF10
plot([FF25 FF25],[yl(1) yl(2)],'LineWidth',2,'color','red','Linestyle','--') %FF25
plot([FF75 FF75],[yl(1) yl(2)],'LineWidth',2,'color','red','Linestyle','--') %FF75
plot([FF90 FF90],[yl(1) yl(2)],'LineWidth',2,'color','red','Linestyle','--') %FF90

%title([' Mean= ',num2str(meanpv), ' Median = ', num2str(pvmedian), ' IQR = ', num2str(IQR), ' P95 = ',num2str(pv95th), ' Max = ',num2str(pvmax)])

%loc = find(pv > level) ;
%xlabel(['Fixedthresh = ',num2str(diffthresh), 'Varthresh = ',num2str(N95), '    FDFfixed = ',num2str(FDFfixed), '    FDFvar = ',num2str(FDFvar)])

% (2) R2* histogram
subplot(2,2,2)

    %Specific bin widths
    xbins2=0:0.5:10;

hist(double(pv2),xbins2)
xlim([0 10]);
hhist = gca ;
hold on
title(['R2star Histogram'])
xlabel(['R2star'])
ylabel(['Number of voxels'])
set(gca,'FontSize',16);

%Plot key values
yl=ylim;
plot([R2median R2median],[yl(1) yl(2)],'LineWidth',2,'color','red','Linestyle','--') %median
plot([R2_10 R2_10],[yl(1) yl(2)],'LineWidth',2,'color','red','Linestyle','--') %FF10
plot([R2_25 R2_25],[yl(1) yl(2)],'LineWidth',2,'color','red','Linestyle','--') %FF25
plot([R2_75 R2_75],[yl(1) yl(2)],'LineWidth',2,'color','red','Linestyle','--') %FF75
plot([R2_90 R2_90],[yl(1) yl(2)],'LineWidth',2,'color','red','Linestyle','--') %FF90

% (3) Display data

t=uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames);
%uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
%t.Position= [70 30 190 220]

% (4) FF/R2 scatterplot
subplot(2,2,4)
dscatter_mod(pv,pv2,'PLOTTYPE','CONTOUR','SMOOTHING',30)

xlim([0 100]);
ylim([0 10]);
set(gca,'FontSize',16);
set(gca,'LineWidth',2);
xlabel(['Fat fraction'])
ylabel(['R2*'])
title(['FF / R2* Scatterplot'])
hold off


function rois = get_selected_rois(handles)
lstr = get(handles.listbox_rois,'String') ;
lval = get(handles.listbox_rois,'Value') ;

nsel = length(lval) ;
rois = zeros([1 nsel]) ;

for isel = 1:nsel
    sstr = lstr{lval(isel)} ;
    if isstrprop(sstr(1:4),'digit')
      rois(isel) = strread(sstr(1:4),'%u') ;
    else
        rois = [] ; % 'ROI information'
    end
end
   
% --- Executes on selection change in popupmenu_activefig.
function popupmenu_activefig_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_activefig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_activefig contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_activefig
val = get(hObject,'Value') ;
str = get(hObject,'String') ;

fig = sscanf(str{val},'Fig: %d%s') ;
hf_num = fig(1) ;

if verLessThan('matlab','8.4.0')
    hf = hf_num ;
else
    r = groot ;
    hf = findobj(r.Children,'Type','Figure','Number',hf_num);
end

UDcont = get(handles.figure1,'UserData') ;
UDcont.hf = hf ;

UD = get(hf,'UserData') ;
islice = UD.disp_slice ;

if UD.nslice > 1
  set(handles.slider_slice,'Min',1,'Max',UD.nslice,'SliderStep',[1/(UD.nslice-1) max(0.1,1/(UD.nslice-1)) ],...
      'Value', islice) ;
else
    set(handles.slider_slice,'Visible','off')
end

set(handles.slider_slice,'Value',islice) ;
set(handles.edit_slice,'String',num2str(islice) ) ;


set(handles.figure1,'UserData',UDcont) ;

update_listbox(handles) ;

% set possible link figures
hfigs = findobj(0,'Tag','roianalfig') ;

if isempty(hfigs)
    turn_off_link(handles) ;
	return
end

lfig = 0 ;
for ifig = 1:length(hfigs)
    UDlink = get(hfigs(ifig),'UserData') ;
    if verLessThan('matlab','8.4.0')
        hfigs_num = hfigs(ifig) ;
    else
        hfigs_num = hfigs(ifig).Number ;
    end
    
    % error here if previous unwanted roianal fig is still open
    if UDlink.nslice == UD.nslice && hfigs_num ~= hf_num
        lfig = lfig + 1 ;
        name = get(hfigs(ifig),'Name') ;
        
        st{lfig} = ['Fig: ',num2str(hfigs_num),' ',name];
        
    end
end

if lfig > 0
    set(handles.popupmenu_linkfig,'String',st) ;
    set(handles.popupmenu_linkfig,'value',1) ;
    turn_on_link(handles)
else
    turn_off_link(handles)
end



% --- Executes during object creation, after setting all properties.
function popupmenu_activefig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_activefig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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
state = get(hObject,'Value') ;
if state == 1
    set(handles.popupmenu_linkfig,'Enable','on')
    set(handles.pushbutton_newlinkedroi,'Enable','on')
    draw_slice(handles) % forces linked to same slice as active
else
    set(handles.popupmenu_linkfig,'Enable','off')
    set(handles.pushbutton_newlinkedroi,'Enable','off')
end

% --- Executes on selection change in popupmenu_linkfig.
function popupmenu_linkfig_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_linkfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_linkfig contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_linkfig


% --- Executes during object creation, after setting all properties.
function popupmenu_linkfig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_linkfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_newlinkedroi.
function pushbutton_newlinkedroi_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_newlinkedroi (see GCBO)
% eventdata  reserved - to be defined in, a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% first add roi as per New ROI button
UDcont = get(handles.figure1,'UserData') ;
hf = UDcont.hf ;
UD = get(hf,'UserData') ;


hroi_this = draw_imroi(handles, UD.ha) ;
setColor(hroi_this,'yellow')

if isfield(UD,'hroi')
    hrois = UD.hroi ;
    nroi = length(hrois) ;
else
    nroi = 0 ;
end

islice = str2double(get(handles.edit_slice,'String')) ;
labstr = get(handles.edit_roilabel,'String') ;

UD.roislice(nroi+1) = islice ;
UD.hroi{nroi+1} = hroi_this ;
UD.roilabel{nroi+1} = labstr ;

set(hf,'UserData',UD)
update_listbox(handles, nroi+1)

% add copy of this ROI to linked figure and link each
pos_main = getPosition(hroi_this) ;

val = get(handles.popupmenu_linkfig,'Value') ;
str = get(handles.popupmenu_linkfig,'String') ;
fig = sscanf(str{val},'Fig: %d%s') ;
hflink = fig(1) ;
UDlink = get(hflink,'UserData') ;

% Add copy to linked figure
hroi_linked = draw_imroi(handles, UDlink.ha,pos_main) ;
setColor(hroi_linked,'yellow')
%update linked figure's list of ROIs but not list box as this is for
%current active figure
if isfield(UDlink,'hroi')
    hrois_link = UDlink.hroi ;
    nroi_link = length(hrois_link) ;
else
    nroi_link = 0 ;
end
UDlink.roislice(nroi_link+1) = islice ;
UDlink.hroi{nroi_link+1} = hroi_linked ;
UDlink.roilabel{nroi_link+1} = labstr ;
set(hflink,'UserData',UDlink)

% each ROI contains a UserData with its linked ROI's handle in UD.hlinkedROI
UDroilink.hlinkedROI = hroi_this ;
UDroi.hlinkedROI = hroi_linked ;

set(hroi_linked,'UserData',UDroilink)
set(hroi_this,'UserData',UDroi)

id1 = addNewPositionCallback(hroi_this,@ (pos) roicallb(pos, hroi_this)) ;
id2 = addNewPositionCallback(hroi_linked,@ (pos) roicallb(pos, hroi_linked)) ;

function roicallb(pos, hpoly)

% hr2 = impoly ;
% hr = impoly ;
% UD.hlinkedROI = hr2 ;
% set(hr,'UserData',UD)
%
% id2 = addNewPositionCallback(hr,@ (pos) roicallb(pos, hr)) ;

UD = get(hpoly,'UserData');

if isfield(UD,'hlinkedROI')
    setPosition(UD.hlinkedROI,pos) ;
end


function turn_off_link(handles)
set(handles.popupmenu_linkfig,'Enable','off')
set(handles.pushbutton_newlinkedroi,'Enable','off')
set(handles.checkbox1,'Value',0)
set(handles.checkbox1,'Enable','off')

function turn_on_link(handles)
set(handles.checkbox1,'Enable','on')

function update_activefig(handles)
% updates active figure popupmenu
UDcont = get(handles.figure1,'UserData') ;
hf = UDcont.hf ;
hfigs = findobj(0,'Tag','roianalfig') ;

if length(hfigs)==0
	warning(' There are no roianal figs open')
end

for ifig = 1:length(hfigs)
  name = get(hfigs(ifig),'Name') ;
  if verLessThan('matlab','8.4.0')
      st{ifig} = ['Fig: ',num2str(hfigs(ifig)),' ',name];
  else
      st{ifig} = ['Fig: ',num2str(hfigs(ifig).Number),' ',name];
  end
end

set(handles.popupmenu_activefig,'String',st) ;
val = find(hfigs== hf) ;
set(handles.popupmenu_activefig,'value',val) ;

roianal('popupmenu_activefig_Callback',handles.popupmenu_activefig, [], handles) ;


% --- Executes on button press in pushbutton_update.
function pushbutton_update_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_listbox(handles)

% --- Executes on button press in radiobutton_contour.
function radiobutton_contour_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_contour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_contour
state = get(hObject,'Value') ;
UDcont = get(handles.figure1,'UserData') ;
UD = get(UDcont.hf,'UserData') ;

if state == 1
    % calculate and display contours
    currsl = str2num(get(handles.edit_slice,'String')) ;
    level = str2num(get(handles.edit_threshold,'String')) ;
    XData = get(UD.him(1),'XData') ;
    YData = get(UD.him(1),'YData') ;
    X = linspace(XData(1),XData(2),size(UD.vol,2)) ;
    Y = linspace(YData(1),YData(2),size(UD.vol,1)) ;
    C = contourc(X,Y,double(UD.vol(:,:,currsl)),double([level level]));
    
    ncC = size(C,2) ;
    ic = 1 ;
    hold(UD.ha,'on')
    while(ic <=ncC )
        np = C(2,ic) ;
        Xp = C(1,ic+1:ic+np) ;
        Yp = C(2,ic+1:ic+np) ;
        
        plot(UD.ha,Xp,Yp,'Tag','contour')
        ic = ic+1+np ;
    end
    hold(UD.ha,'off')
else
    hcont = findobj(UDcont.hf, 'Tag','contour') ;
    delete(hcont)
end


function edit_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_threshold as text
%        str2double(get(hObject,'String')) returns contents of edit_threshold as a double


% --- Executes during object creation, after setting all properties.
function edit_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_finish.
function pushbutton_finish_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_finish (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UD = get(handles.figure1,'UserData') ;
if UD.ismodal == true
    uiresume(handles.figure1) ;
end

set(hObject,'Visible','off')

function roistats = calc_roistats(handles, roitype, pos, islice)
UDcont = get(handles.figure1,'UserData') ;
hf = UDcont.hf ;
UD = get(hf,'UserData') ;

switch roitype
    case 'imline'
        x1 = pos(1,1); x2 = pos(2,1) ; y1=pos(1,2); y2 = pos(2,2) ;
        llen = sqrt(( (x1-x2).^2 + (y1-y2).^2)) ;
        np = round(4*llen) ; % 1/4 mm
        dx = llen/np;
        
        hims = UD.him ;
        him = hims(islice) ;
        XData = get(him,'XData') ;  YData = get(him,'YData') ;
        
        [cx,cy,c] = improfile(XData,YData,UD.vol(:,:,islice),[x1 x2],[y1 y2],np) ;
        dist = sqrt((cx-x1).^2+(cy-y1).^2)-llen/2 ;
        
        % Area above threshold
        % threshold analysis
        mint = min(c) ;
        maxt = max(c) ;
        nt = 100 ;
        A = zeros([1 nt]) ;
        ts = linspace(mint,maxt,nt) ;
        for it =1:nt
            ct = c-ts(it) ;
            ct(ct<0)=0;
            A(it) = sum(ct)*dx;
        end
        
        roistats.threshs = ts ;
        roistats.areas = A ;
        roistats.dist = dist ;
        roistats.profile = c ;
        roistats.llen = llen ;
    otherwise
        roistats = [] ;

end

        
        

% --- Executes on button press in pushbutton_linestats.
function pushbutton_linestats_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_linestats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UDcont = get(handles.figure1,'UserData') ;
hf = UDcont.hf ;
UD = get(hf,'UserData') ;

hrois = [UD.hroi];

rois = get_selected_rois(handles) ;
nsel = length(rois) ;

if verLessThan('matlab','8.4.0')
    hf_num = hf ;
else
    hf_num = hf.Number ;
end

descrip = ['Profiles on fig ', num2str(hf_num),': ',get(hf,'Name')] ;
    
hprofile = figure('Name',descrip, 'Tag', 'roianal_lineprofile') ;
haprofile = gca ;
hold on
    
hthresh = figure('Name',['thresh ',descrip], 'Tag', 'roianal_athresh') ;
hathresh = gca ;
hold on

labstrc = cell([1 nsel]);
roitype = cell([1 nsel]);
slice = zeros([1 nsel]) ;

    
for isel = 1:nsel
    roitype = get(hrois{rois(isel)},'Tag') ;
    switch roitype
        case 'imline'
            pos = getPosition(UD.hroi{rois(isel)}) ;
            slice(isel) = UD.roislice(rois(isel)) ;
            
            lab = UD.roilabel(rois(isel)) ;
            labstrc{isel} = ['sl ',num2str(slice(isel)),': ',lab{1}];
            
            roistats = calc_roistats(handles,roitype, pos, slice(isel)) ;
            
            hp = plot(haprofile,roistats.dist,roistats.profile) ;
            col = getColor(UD.hroi{rois(isel)}) ;
            set(hp,'Color',col) % set to same color as line on image
            
            hpt = plot(hathresh,roistats.threshs,roistats.areas);
            set(hpt,'Color',col)
        otherwise
    end
end

grid(hathresh,'on')
title(hathresh,['Areas vs threshold. ',descrip])
xlabel(hathresh,'Threshold')
ylabel(hathresh,'Area')
legend(hathresh,labstrc{:})

grid(haprofile,'on')
title(haprofile,descrip)
xlabel(haprofile,'Distance')
ylabel(haprofile,'Profile')
legend(haprofile,labstrc{:})

return
if nlines > 0
    islice = str2num(get(handles.edit_slice,'String') ) ;
    
    descrip = ['Profiles on fig ',...
        num2str(hf),': ',get(hf,'Name')] ;
    
    hprofile = figure('Name',descrip) ;
    haprofile = gca ;
    hold on
    
    hthresh = figure('Name',['thresh ',descrip]) ;
    hathresh = gca ;
    hold on
    
    for iline = 1:nlines
        pos = getPosition(hlines{iline}) ;
        x1 = pos(1,1); x2 = pos(2,1) ; y1=pos(1,2); y2 = pos(2,2) ;
        col = getColor(hlines{iline}) ;
        
        llen = sqrt(( (x1-x2).^2 + (y1-y2).^2)) ;
        np = round(4*llen) ; % 1/4 mm
        dx = llen/np;
        
        hims = UD.him ;
        him = hims(islice) ;
        XData = get(him,'XData') ;
        YData = get(him,'YData') ;
        
        [cx,cy,c] = improfile(XData,YData,UD.vol(:,:,islice),[x1 x2],[y1 y2],np) ;
        
        hp = plot(haprofile,sqrt((cx-x1).^2+(cy-y1).^2)-llen/2,c) ;
        set(hp,'Color',col) % set to same color as line on image
        
        % threshold analysis
        mint = min(c) ;
        maxt = max(c) ;
        nt = 100 ;
        A = zeros([1 nt]) ;
        ts = linspace(mint,maxt,nt) ;
        for it =1:nt
            ct = c-ts(it) ;
            ct(ct<0)=0;
            A(it) = sum(ct)*dx;
        end
        hpt = plot(hathresh,ts,A);
        set(hpt,'Color',col)
            
    end
    
    grid(hathresh,'on')
    title(hathresh,['Areas vs threshold. ',descrip])
    xlabel(hathresh,'Threshold')
    ylabel(hathresh,'Area')
    
    grid(haprofile,'on')
    title(haprofile,descrip)
    xlabel(haprofile,'Distance')
    ylabel(haprofile,'Profile')
end




% --- Executes on selection change in popupmenu_roitype.
function popupmenu_roitype_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_roitype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_roitype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_roitype


% --- Executes during object creation, after setting all properties.
function popupmenu_roitype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_roitype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_roilabel_Callback(hObject, eventdata, handles)
% hObject    handle to edit_roilabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_roilabel as text
%        str2double(get(hObject,'String')) returns contents of edit_roilabel as a double


% --- Executes during object creation, after setting all properties.
function edit_roilabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_roilabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_applylabel.
function pushbutton_applylabel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_applylabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rois = get_selected_rois(handles) ;
UDcont = get(handles.figure1,'UserData') ;
UD = get(UDcont.hf,'UserData') ;

labstr = get(handles.edit_roilabel,'String') ;
for iroi = 1:length(rois)
   UD.roilabel{rois(iroi)} = labstr ;
end

set(UDcont.hf,'UserData',UD) 
update_listbox(handles) 


% --- Executes on button press in pushbutton_report.
function pushbutton_report_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_report (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ispref('roianal','save_dir')
    defdir = getpref('roianal','save_dir');
else
    defdir = [] ;
end

% in future could offer choice of report format
[fn,pn,FilterIndex] = uiputfile({'*.pdf','PDF report'},'Enter filename',defdir) ;

if fn == 0
    warning(['No file specified'])
    return
end

ffn = fullfile(pn,fn) ;

roianal_list = get(handles.listbox_rois,'String') ;
assignin('base','roianal_list',roianal_list)

%generate montage figure that will go into report
draw_montage(handles)

if exist('RptgenML.CReport','class')
    [RptgenML_CReport1] = buildroianal ;
    disp(['Writing report to: ',ffn])
    report(RptgenML_CReport1,['-o',ffn],'-noview','-debug') ;
else
    disp(['MATLAB Report toolbox not present. No report created'])
end


% --- Executes on button press in pushbutton_CopyToAll.
function pushbutton_CopyToAll_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_CopyToAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

iroi = get_selected_rois(handles) ;
if length(iroi) == 1
    UDcont = get(handles.figure1,'UserData') ;
    UD = get(UDcont.hf,'UserData') ;
    hrois = UD.hroi ;
    nrois = length(hrois) ;
    
    posn = getPosition(hrois{iroi}) ;
    col = getColor(hrois{iroi}) ;
    roi_slice = UD.roislice(iroi) ;
    labstr = UD.roilabel{iroi} ; 
    
    
    this_slice_display = str2double(get(handles.edit_slice,'String')) ;
    
    nslice = size(UD.vol , 3) ;
    
    for islice = 1: nslice
        if islice ~= roi_slice
            
            draw_slice_one(UD, islice, UDcont.hf) ;
            
            hroi_new = draw_imroi(handles, UD.ha, posn) ;
            setColor(hroi_new,col)
            
            nrois = nrois+1 ;
            UD.roislice(nrois) = islice ;
            UD.hroi{nrois} = hroi_new ;
            UD.roilabel{nrois} = labstr ;
        end
    end
    
    draw_slice_one(UD, this_slice_display, UDcont.hf) ;
    
    set(UDcont.hf,'UserData',UD) ;
    update_listbox(handles)        
else
    disp(['Copy To All requires exactly one ROI selected.'])
end


% --- Executes on button press in pushbuttonMeanVsT.
function pushbuttonMeanVsT_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonMeanVsT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rois = get_selected_rois(handles) ;
if length(rois) == 0 
    disp(['Need at least one ROI selected'])
else
    UDcont = get(handles.figure1,'UserData') ;
    UD = get(UDcont.hf,'UserData') ;
    hrois = UD.hroi ;
    
    labs = [UD.roilabel]  ;
   
    ulabs = unique(labs(rois)) ;
    
    hfig = figure('Name','roianal: Mean Area vs Time','Tag','roianal_time') ;
    hold on
    for ilab = 1:length(ulabs)
       loc_in_all = find(strcmp(ulabs{ilab}, labs)) ;
       loc = intersect(loc_in_all, rois) ;
       
       sls = [UD.roislice(loc)] ;
       [sls_sort, idx] = sort(sls) ;
       
       for isls = 1:length(sls_sort)
         iroi = loc(idx(isls)) ;
         if ~isempty(hrois{iroi})
            sl = UD.roislice(iroi) ;
            roitype = get(hrois{iroi},'Tag') ;
            switch roitype
                case 'impoly'
                  BW = createMask(hrois{iroi},UD.him(sl)) ;
                  stats = regionprops(BW,UD.vol(:,:,sl),'MeanIntensity') ;
             
                  mns(isls) = stats.MeanIntensity ;
            end
         end
       end
        
       if isfield(UD,'TinSeries') && (length(UD.TinSeries) == length(mns))
          hp = plot(UD.TinSeries, mns,'o-') ;
       else
          hp = plot(mns,'o-')
       end
       col = getColor(hrois{iroi}) ;
       set(hp,'Color',col)
       
    end
    legend(ulabs{:})
    grid
    
end


% --- Executes on button press in pushbutton_snapshot.
function pushbutton_snapshot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_snapshot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


hfig = figure('Name','Snapshot','Tag','roianal_snapshot')
hax = gca ;

slicestr = get(handles.edit_slice,'String') ;
set(hfig,'Name',['Snapshot slice: ',num2str(slicestr)])

UDcont = get(handles.figure1,'UserData') ;
hf = UDcont.hf ;
UD = get(hf,'UserData') ;

frame = getframe(UD.ha) ;
[X,map] = frame2im(frame) ;
imshow(X,map,'Parent',hax,'Border','tight')

if verLessThan('matlab','8.4.0')
    hfig_num = hfig ;
else
    hfig_num = hfig.Number ;
end

disp(['Snapshot created figure: ',num2str(hfig_num)])





function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

width = str2num(get(hObject,'String')); %Gets width from textbox
handles.edit5.Value = width

%set(handles.edit5,'Value',width) ;

disp(width)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
    

end




% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
widthINNER = str2num(get(hObject,'String')); %Gets width from textbox
handles.edit8.Value = widthINNER

disp(widthINNER)
% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

normP95 = str2num(get(hObject,'String')); %Gets width from textbox
handles.edit9.Value = normP95;

%set(handles.edit5,'Value',width) ;

disp(perc95)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



