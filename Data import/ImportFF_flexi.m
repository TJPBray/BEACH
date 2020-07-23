%Import mDixon Quant FF map and convert to usable FF and R2* volumes
%Colourmap as normal

clear all

%Select folder
foldername=uigetdir;
folderinfo=dir(foldername);
n=numel(folderinfo)-2;

%Acquire images from DICOM files
for k=1:n
    filename=folderinfo(k+2).name;
    filename=fullfile(foldername,filename);
    image=dicomread(filename);
    stack(:,:,k)=image;
    
    try
    info(:,:,k)=dicominfo(filename);
    catch
        ;
    end
    
end

stack=double(stack);

%Sort into component images
dim=size(stack);
depth=dim(3)/4;


%Use pos to assign the position of the FF map in the stack (1,2,3,4)    
pos=2;

%Generate FF maps using rescale info
    FFinfo=info(1,1,(pos-1)*depth+1);
    
    %Get rescale info if possible, otherwise assign hard value
    try
        RescaleSlope=FFinfo.RescaleSlope;
    catch
        warning('Problem interpreting DICOM info.  Assigning a fixed rescale slope of 0.0733. Usually accurate for mDixon Quant.');
        RescaleSlope=0.07326;
    end
    
    %If rescale slope field is present but empty
    if isempty(RescaleSlope)==1
        RescaleSlope=0.07326;
    else;
    end
    
%find the rescale intercept 
    %Get rescale info if possible, otherwise assign hard value
    try
        RescaleIntercept=FFinfo.RescaleIntercept;
    catch
        warning('Problem interpreting DICOM info.  Assigning a fixed rescale intercept of -100. Usually accurate for mDixon Quant.');
        RescaleIntercept=-100;
    end

    %If rescale slope field is present but empty
    if isempty(RescaleIntercept)==1
        RescaleIntercept=-100;
    else;
    end
    
%Generate FF maps with rescale slope value
      FF=stack(:,:,((pos-1)*depth+1):(pos*depth));
%     FF=double(FF+RescaleIntercept);
%     FF=FF*RescaleSlope;
      FF=double(FF);
      FF = FF*RescaleSlope;
      FF=FF+RescaleIntercept;
    
    %Remove extreme outliers
    FF(FF>110)=110; 
    FF(FF<-10)=-10;
 
%Generate water images
Water=stack(:,:,(1:depth));

%Generate R2* maps
    T2star=double(stack(:,:,(3*depth+1):(4*depth)));
    
    %Replace 0s with NaN (prevents bright background in final image)
    T2star(T2star==0)=NaN;
    R2star=1000./T2star;
    
    %Remove extreme values
    R2star(R2star>10)=10;
    
% newanal2(FF) change this to use freehand ROI
newanal2(FF)